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

        res_out<-list()
        res_out$sprectra_ft<- sprectra_ft
        res_out$cluster_ft<- cluster_ft

        return(res_out)

} #end function 2



#3. Function to create msp and mgf files

ramclustr_to_msp_mgf<-function(ramclustr_out= ramclustr_outloc,
                               ramclustr_ft= ramclustr_res$sprectra_ft,
                               ramclustr_clusters= ramclustr_res$cluster_ft,
                               study_id= study_id){



#Create msp and mgf output files
        msp_output<- paste(ramclustr_out, "4a-", study_id, "-GCHRMS_MSP-Spectra.msp", sep="")
        mgf_output<- paste(ramclustr_out, "4b-", study_id, "-GCHRMS_MGF-Spectra.mgf", sep="")


        for(ii in 1:nrow(ramclustr_ft)){
          clust_ii<- ramclustr_ft[ii,1] #Select cluster #

          ms_spectra<- ramclustr_clusters[ramclustr_clusters$Cluster_name==clust_ii,]  %>%              #Select mz for cluster
                        .[,c(2, 14 + which(colSums(.[, -c(1:14)]) == max(colSums(.[,-c(1:14)]))))] %>%     #Select sample with max column sum
                        rename_at(2, ~"intensity") %>%              #Rename intensity column
                        arrange(desc(intensity))  %>%               #Arrange in order
                        mutate(intensity = round(intensity / max(intensity) * 100, 2)) %>%  #Normalize by max
                        filter(intensity >= 1)         #Remove features with intensity %< 1



          if(nrow(ms_spectra) > 1){
            ms_entropy<- msentropy::calculate_spectral_entropy(as.matrix(ms_spectra))/log(nrow(ms_spectra))
          } else  {
            ms_entropy<-NA
          }


          #Code for writing to MSP file

            cat(
            paste("NAME", clust_ii, sep=": "),
              "msLevel: 1",
            paste("RETENTIONTIME", round(as.numeric(ramclustr_ft= ramclustr_res$sprectra_ft[ii,2]),2), sep=": "),
            paste("RETENTIONINDEX", round(as.numeric(ramclustr_ft= ramclustr_res$sprectra_ft[ii,2]),2), sep=": "),
              "INSTRUMENTTYPE: GC-HRMS",
              "INSTRUMENT: Thermo Exploris240",
              "IONMODE: Positive",
              "IONIZATION: EI",
              "Spectrum_type: in-source",
            paste("Notes: Date processed", Sys.Date(), sep= " "),
            paste("COMMENT: ",
                  "Study ID: ", study_id, "; ",
                  "Normalized Spectral Entropy: ", round(ms_entropy,3), sep=""),
            paste("Num Peaks: ", nrow(ms_spectra), sep=""),

            sep = "\n", file = msp_output, append = TRUE)

          #Add spectra
          data.table::fwrite(x =ms_spectra,
                             file = msp_output,
                             sep="\t",
                             col.names=F,
                             append=T)

          #Add empty line after spectra
          cat("",
              sep = "\n", file = msp_output, append = TRUE)




          ###Create mgf file
          # Add all meta data

          cat("BEGIN IONS",
              paste("NAME", clust_ii, sep= "="),
                "CHARGE=1+",
              paste("RTINSECONDS", round(as.numeric(ramclustr_ft= ramclustr_res$sprectra_ft[ii, 2]), 2), sep= "="),
              paste("TITLE=",
                    "Study ID: ", study_id, "; ",
                    "Normalized Spectral Entropy: ", round(ms_entropy,3), sep=""),
              sep = "\n", file = mgf_output, append = TRUE)

          #Add spectra
          data.table::fwrite(x = ms_spectra,
                             file = mgf_output,
                             sep = "\t",
                             col.names = F,
                             append = T)

          #Add empty line after spectra
          cat("END IONS",
              "",
              sep = "\n", file = mgf_output, append = TRUE)





        }

        }




