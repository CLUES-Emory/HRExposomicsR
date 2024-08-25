#1. Function to create peak groups matrix for standard based alignment
retcor_group_matrix<- function( istd_df= istds,
                                xcms_grp1= step_1b_res,
                                sample_align= p_final,
                                mz_error = 5,
                                time_error = 5
                  ) {

    #Internal standard database
    istds<- istd_df

    #Create matrix of retention times for peaks in all samples
    rt_mat<- featureValues(xcms_grp1, value="rt") %>%
      cbind(featureDefinitions(xcms_grp1)[, c("mzmed", "mzmin", "mzmax", "rtmed","rtmin", "rtmax")],
            featureDefinitions(xcms_grp1)$rtmax - featureDefinitions(xcms_grp1)$rtmin, .) %>%
      add_column(., d = (featureDefinitions(xcms_grp1)$mzmax - featureDefinitions(xcms_grp1)$mzmin)/featureDefinitions(xcms_grp1)$mzmin*10^6, .after = "mzmax")

    colnames(rt_mat)[4]<- "group_mz_range_ppm"
    colnames(rt_mat)[8]<- "group_time_range"


    #Match internal standards to detected peaks using xMSanalyzer
    overlap_res<- xMSanalyzer::find.Overlapping.mzs(istds[,c("mz", "Alkane_RI")], rt_mat[,c("mzmed", "rtmed")],
                                       mz.thresh = mz_error, time.thresh = time_error, alignment.tool = NA)

    overlap_res<-cbind(istds[overlap_res$index.A,1], overlap_res)

    overlap_res$mass_error_ppm<- round((overlap_res$mz.data.A - overlap_res$mz.data.B) / overlap_res$mz.data.A * 10^6, 3)
    overlap_res$time.difference<- round(overlap_res$time.difference, 2)

    colnames(overlap_res)[1]<- "ISTD_name"

    #Combine to create final table
    rt_mat<-cbind(overlap_res, rt_mat[overlap_res$index.B, ])
    rt_mat<-rt_mat[,c(1:17, 17 + c(sample_align))]

    #Count non-detects in samples and remove any with greater than
    na_rm<-colSums(is.na(rt_mat[, -c(1:17)])) / nrow(rt_mat)
    na_rm<-which(na_rm > 0.5)

    #Remove samples with greater than 50% missing
    if(length(na_rm) != 0){
      rt_mat<- rt_mat[, -c(na_rm + 17)]
      sample_align<- sample_align[-na_rm]
    }

    return(list(rt_df= rt_mat[,-c(1:17)], align_index= sample_align))

} #End of function 1








