#xcms wrapper function to generate
gc_xcms_wrapper<- function( xcms_params= xcms_params,
                            merge_peaks_rt= 25,
                            mapfile= mapfile,
                            subset_alignment= TRUE,
                            sample_class_subsets= subset_samples,
                            std_alignment= TRUE,
                            istds= istds,
                            xcms_outloc= xcms_outloc,
                            study_id= study_id
) }

    print("Stage 1: Peak detection and alignment using xcms")

    print("Performing Step 1: Peak Detection")
    #Step 1, peak detection
    #Define CentWave parameterds
    xcms::cwp <- CentWaveParam(
      ppm= xcms_params$cwp_ppm,
      peakwidth= xcms_params$cwp_peakwidth,
      snthr= xcms_params$cwp_snthr,
      mzdiff= xcms_params$cwp_mzdiff,
      noise= xcms_params$cwp_noise,
      prefilter= xcms_params$cwp_prefilter,
      mzCenterFun= xcms_params$cwp_mzCenterFun,
      integrate= xcms_params$cwp_integrate,
      fitgauss= xcms_params$cwp_fitgauss,
      extendLengthMSW= xcms_params$cwp_extendLengthMSW)


    #Detect peaks using cwp
    step_1_res <- xcms::findChromPeaks(raw_data, param = cwp)


    #Step 1a, merge split peaks
    print("Performing Step 1a: Merge Split Peaks")

    #Merge split peaks using parameters in mpp
    mpp<- xcms::MergeNeighboringPeaksParam(
      expandRt= merge_rt,
      expandMz= xcms_params$mrg_expandMz,
      ppm= xcms_params$mrg_ppm,
      minProp= xcms_params$mrg_minProp)


    step_1a_res <- xcms::refineChromPeaks(step_1_res, param= mpp)

    #Step 1b, First grouping before retention time correction
    print("Performing Step 1b: First Grouping")

    pdp.1 <- xcms::PeakDensityParam(sampleGroups= mapfile[,3],
      minFraction= xcms_params$grp1_minFraction,
      bw= xcms_params$grp1_bw,
      ppm= xcms_params$grp1_ppm,
      binSize= xcms_params$grp1_binSize)

    step_1b_res <- xcms::groupChromPeaks(step_1a_res, param = pdp.1)

    #saveRDS(step_1b_res, paste(xcms_outloc, study_id, "-Step01_XCMSnExp.rds", sep=""))


    #######
    #Step 2; Retention time correction
    print("Performing Step 2: Retention time Correction")


    #Select samples used for alignment (should be bio samples only)
    if(subset_alignment == TRUE) {
      p_final<- sort(which(mapfile[, 4] %in% sample_class_subsets))

      }	else	{
      p_final<- 1:nrow(mapfile)

      }


    #Perform standard based retention time correction (if std_alignment= TRUE)
    if(std_alignment == TRUE) {
      istd_grp_matrix<- retcor_group_matrix(istd_df= istds,
        xcms_grp1= step_1b_res,
        align_index= p_final)

      #Retention time correction parameters using subset and std alignment
      retcor_params<- PeakGroupsParam(
        minFraction = xcms_params$rtcor_minFraction,
        extraPeaks = xcms_params$rtcor_extraPeaks,
        peakGroupsMatrix= as.matrix(istd_grp_matrix$rt_mat[,-c(1:17)]),
        smooth = xcms_params$rtcor_smooth,
        span = xcms_params$rtcor_span,
        subset = istd_grp_matrix$align_index,
        subsetAdjust = xcms_params$rtcor_subsetAdjust,
        family = xcms_params$rtcor_family)


    } else  {
      retcor_params<- xcms::PeakGroupsParam(
        minFraction = xcms_params$rtcor_minFraction,
        extraPeaks = xcms_params$rtcor_extraPeaks,
        smooth = xcms_params$rtcor_smooth,
        span = xcms_params$rtcor_span,
        subset = p_final,
        subsetAdjust = xcms_params$rtcor_subsetAdjust,
        family = xcms_params$rtcor_family)

    }

    #Run retention time correction
    step_2_res<- xcms::adjustRtime(object= step_1b_res, param = retcor_params)


    #Create adjustment plot
    res_c<-colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nrow(mapfile))
    names(res_c) <- mapfile[, 1]

    pdf(paste(xcms_outloc, study_id, "-Retention_time_deviation.pdf", sep=""))
    xcms::plotAdjustedRtime(step_2_res, col = res_c)

    dev.off()


  ########

    #Step 3: Peak grouping
    print("Performing Step 3: Second Grouping")


    #Grouping 2 parametere
    pdp <- xcms::PeakDensityParam(sampleGroups =mapfile[,3],
                            minFraction = xcms_params$grp2_minFraction,
                            bw = xcms_params$grp2_bw,
                            ppm =xcms_params$grp2_ppm,
                            binSize=xcms_params$grp2_binSize)


    #Secong peak grouping after rt adjustment
    step_3_res <- xcms::groupChromPeaks(step_2_res, param = pdp)


    #######
    #Step 4: Fill in missing values
    print("Performing Step 4: Missing Value Filling")


    step_4_res <- fillChromPeaks(step_3_res, param = ChromPeakAreaParam())


    #Calculate percent decrease in zeroes after missing value filling
    print(paste("Percent decrease in zeroes from missing value filling: ",
                round((sum(is.na(featureValues(step_3_res))) - sum(is.na(featureValues(step_4_res)))) / sum(is.na(featureValues(step_3_res))) * 100, 2),
                "%", sep=""))


    #######
    ### Step 5: Create and output final feature table
    print("Performing Step 5: Exporting Aligned Feature Table")

    ## Get feature definitions and intensities
    featuresDef <- xcms::featureDefinitions(step_4_res)
    featuresIntensities <- xcms::featureValues(step_4_res, value = "into")

    ## generate data table
    feature_table <- merge(featuresDef, featuresIntensities, by = 0, all = TRUE)
    feature_table <- feature_table[, !(colnames(feature_table) %in% c("peakidx"))]


    #Replace NA's with zeroes
    feature_table[is.na(feature_table)]<-0


    #Write results table
    file_ii<-match(paste(mzML_files_BASE, ".mzML", sep=""), colnames(feature_table))
    feature_table<-feature_table[, c(1:12, file_ii)]


    #Calculate feature summaries
    #Number detects
    feature_table[, 9]<- rowSums(feature_table[, -c(1:12)] != 0)
    #Percent detects
    feature_table[, 10]<- round(feature_table[, 9] / ncol(feature_table[,- c(1:12)])*100,2)
    #Mass range ppm
    feature_table[, 11]<- round((feature_table$mzmax - feature_table$mzmin) / feature_table$mzmax * 10^6, 2)
    #Time range secs
    feature_table[, 12]<- round((feature_table$rtmax -feature_table$rtmin), 2)

    #Create final table and export
    colnames(feature_table)[c(1:2, 5, 9:12)]<- c("Peak_ID", "mz", "time", "Num_Detected", "Percent_Detected", "mz-group_range_ppm", "time-group_range_sec")
    colnames(feature_table)<- gsub(colnames(feature_table), pattern = ".mzML", replacement = "")
    feature_table<- dplyr::arrange(feature_table, mz, time)

    feature_table<- feature_table %>%
      relocate("mz-group_range_ppm", .after = "mzmax") %>%
      relocate("time", .after = "mz") %>%
      relocate("time-group_range_sec", .after = "rtmax")

    out.name<-paste(xcms_outloc, study_id, "-XCMSv", packageVersion("xcms"), "_Untargeted_RAW_Feature_Table.txt", sep="")
    write.table(feature_table, out.name, sep= "\t", row.names=FALSE)

    print("XCMS processing complete")
    print(paste("Number of peaks detected: ", nrow(feature_table), sep= ""))

    return(feature_table)

}

