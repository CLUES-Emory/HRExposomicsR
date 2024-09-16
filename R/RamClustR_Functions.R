##1. Function to create feature table for RamClustR

xcms_to_ramclustr_table<-function(feature_table= feature_table,
                                  ramclustr_out= ramclustr_outloc,
                                  cols_meta= 12,
                                  study_id= study_id


) {

#Create table for RamClustR.
  mz_time<-paste(feature_table$mz, feature_table$time, sep="_")                   # Create single variabler for mz and rt
  sample_id<-colnames(feature_table)[-c(1:cols_meta)]                                    # Sample names

#Create RamClustR feature table
  rcr_input<-t(feature_table[,-c(1:cols_meta)]) %>%          # Remove non-intensity columns from feature table
    cbind(sample_id, .) %>%                 #Add sample ID for each row
    as.data.frame %>%                       #Conver to dataframe
    rename_with(~ c("sample_id", mz_time))  #Rename columns with sample ID and mz_time

#Input file for RamClustR
  rcr_file_name<-paste(out_dir, "1-", study_id, "_RamClustR_Formatted_InputTable.csv", sep="")    # CSV file name
  write.csv(rcr_input, rcr_file_name, row.names=FALSE)


  return(rcr_file_name)


}


##2. Function to run RamClustR and generate clustered feature tables

ramclustR_wrapper<-function(ramclustr_params= ramclustr_params,
                            ramclustr_file= ramclustr_file_path,
                            ramclustr_out= ramclustr_outloc,
                            feature_table= feature_table,
                            cols_meta= 12,
                            study_id= study_id
) {

      #Run ramclustr clustering to identify specral clusters
        t1<-Sys.time()
        RC1 <- ramclustR::ramclustR(ms = ramclustr_file,
                         featdelim = ramclustr_params$featdelim,
                         timepos= ramclustr_params$timepos,
                         st = ramclustr_params$st,
                         sr= ramclustr_params$sr,
                         maxt= ramclustr_params$maxt,
                         blocksize= ramclustr_params$blocksize,
                         deepSplit = ramclustr_params$deepSplit,
                         mult= ramclustr_params$mult,
                         ExpDes= ramclustr_params$ExpDes,
                         sampNameCol= ramclustr_params$sampNameCol,
                         mspout= ramclustr_params$mspout,
                         minModuleSize= ramclustr_params$minModuleSize,
                         mzdec= ramclustr_params$mzdec,
                         linkage=  ramclustr_params$linkage,
                         rt.only.low.n= ramclustr_params$rt.only.low.n,
                         normalize= ramclustr_params$normalize,
                         replace.zeros= ramclustr_params$replace.zeros)


        t2<-difftime(Sys.time(), t1, unit= "hours")
        print(paste("Total time to run RamClustR:", round(t2, 2), "hrs.", sep= " "))
        print("RamClustR finished successfully. Prepararing output tables.")

        #RC1<-readRDS("/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/10-Update_240608/4-Data_Extraction_240807/Stage3-GCHRMS_CLU0102_Spectra_Deconvolution/RamClustR_Results.rds")

        #Create table assigning cluster name to each mz. Output this file to results folder.
        Cluster_name<-paste("C", formatC(RC1$featclus, width = nchar(length(unique(RC1$featclus))), format = "d", flag = "0"), sep = "")

        mz_time<-paste(feature_table$mz, feature_table$time, sep="_")

        ii<- na.omit(match(colnames(RC1$MSdata), mz_time))

        #Create cluster feature table and write output
        cluster_ft<-cbind(feature_table[ii, 1:cols_meta], Cluster_name, RC1$msint, feature_table[-c(1:cols_meta)]) %>%
          arrange(Cluster_name, mz) %>%
          dplyr::rename(Average_int = 14)

        #Write
        write.table(cluster_ft,
                    paste(ramclustr_out, "2-", study_id, "-GCHRMS_ClusterMembership_Feature_Table.txt", sep=""),
                    sep="\t", row.names=FALSE)


        ###Create RamClustR spectral table
        num_per_cluster<- table(cluster_ft$Cluster_name) %>%
          as.data.frame() %>%
          dplyr::rename(Cluster_name = 1, Cluster_size=2)

        #Create spectral clusters feature table with averaged abundance
        SpecAbund <- lapply(RC1$SpecAbund, as.numeric)
        mutate_all(RC1$SpecAbund, as.numeric)

        sprectra_ft<- data.frame(row.names(t(RC1$SpecAbund)), round(RC1$clrt, 2), round(RC1$clrtsd, 2), t(RC1$SpecAbund)) %>%
          dplyr::rename(Cluster_name = 1, Cluster_rt = 2, Cluster_rtsd = 3) %>%
          left_join(., num_per_cluster, by = "Cluster_name") %>%
          relocate(Cluster_size, .after = Cluster_rtsd)

        write.table(sprectra_ft, paste(ramclustr_out, "3-", study_id, "-GCHRMS_PreBatchCorr_SpectraClustered_Feature_Table.txt", sep=""), sep="\t", row.names=FALSE)

        return(sprectra_ft)

}
