
#Function to evaluate batch effects
batch_eval_plots<-function(feature_table= rc_ft, #Feature table from data extraction
                           mapfile= mapfile,             #Mapfile
                           istd_df= istds,                #Internal standard mzs and rt
                           pdf_loc= rcr_outloc,           #Where to save pdfs
                           corr_mode= "Pre",             #Batch correction mode
                           istd_plots= FALSE,              #Whether to plot ISTDs
                           study_id=study_id,
                           cols_meta= 4
){

  #Calculate PCs and table for plotting
  pca_nc<- prcomp(t(feature_table[,-c(1:cols_meta)]), scale=TRUE)$x %>%
    data.frame(cbind(mapfile))  %>%
    mutate(Batch=as.factor(Batch))  %>%
    mutate(, sum_sample_int = rowSums(t(feature_table[,-c(1:cols_meta)])))

  #Statistics for pre-correction batch effects
  pc1_anova_nc<- summary(aov(PC1 ~ Batch, data = pca_nc))
  pc2_anova_nc<- summary(aov(PC2 ~ Batch, data = pca_nc))
  int_anova_nc<- summary(aov(sum_sample_int ~ Batch, data = pca_nc))

  #Create plots for batch effects and sum intensities

  #Colors for plots
  if(batch_num_all<9){
    color_index<- brewer.pal(8, "Dark2")[1:batch_num_all]
  } else {
    color_index<-rep(brewer.pal(8, "Dark2"), 100)[1:batch_num_all]
  }

  p<-c()
  p[[1]]<-ggplot(pca_nc, aes(y=sum_sample_int, x= 1:length(sum_sample_int), fill= as.factor(Batch))) +
    geom_bar(stat = "identity", colour="black", linewidth=0.2) +
    scale_fill_manual(values = color_index) +
    labs(x="Sample number", y= "Sum Intensity", size=20) +
    guides(fill=guide_legend(title="Batch Number")) +
    ggtitle(label= paste("Sum intensity for all samples; p=", round(int_anova_nc[[1]]$`Pr(>F)`[1],4), sep=" ")) +
    theme(plot.title = element_text(hjust = 0.5))

  p[[2]]<-ggplot(pca_nc,aes(x=PC1,y=PC2, color=Batch)) +
    geom_point() +
    scale_color_manual(values = color_index) +
    ggtitle("PC1 vs PC2") +
    theme(plot.title = element_text(hjust = 0.5))

  p[[3]]<-ggplot(pca_nc,aes(x=Batch, y=PC1,  color=Batch)) +
    geom_jitter(position=position_dodge(0.8)) +
    scale_color_manual(values = color_index) +
    labs(x="Batch Number", y= "PC1", size=20) +
    guides(fill=guide_legend(title="Batch Number")) +
    ggtitle(label= paste("PC1 Scores; p=", round(pc1_anova_nc[[1]]$`Pr(>F)`[1],4), sep=" ")) +
    theme(plot.title = element_text(hjust = 0.5))

  p[[4]]<-ggplot(pca_nc,aes(x=Batch, y=PC2,  color=Batch)) +
    geom_jitter(position=position_dodge(0.8)) +
    scale_color_manual(values = color_index) +
    labs(x="Batch Number", y= "PC2", size=20) +
    guides(fill=guide_legend(title="Batch Number")) +
    ggtitle(label= paste("PC2 Scores; p=", round(pc2_anova_nc[[1]]$`Pr(>F)`[1],4), sep=" ")) +
    theme(plot.title = element_text(hjust = 0.5))


  pdf_name<-paste(pdf_loc, "1-", study_id, "_", corr_mode, "-Batch_Correction_ResultsPlots.pdf", sep="")
  ggsave(
    filename = pdf_name,
    plot = marrangeGrob(p, nrow=2, ncol=2),
    width = 10, height = 10
  )
  dev.off()


  ###########
  if(istd_plots){
    ###Plot internal standard intensities in each batch
    #Match internal standards to detected peaks using xMSanalyzer and create table of intensities
    istd_res<- find.Overlapping.mzs(istd_df[,c("mz", "Alkane_RI")], feature_table[,-1],
                                    mz.thresh = 5, time.thresh = 5, alignment.tool=NA) %>%
      cbind(istd_df[.$index.A,1], .)  %>%
      mutate(mass_error_ppm = round((.$mz.data.A - .$mz.data.B)/.$mz.data.A*10^6,3)) %>%
      mutate(time.difference = round(.$time.difference,2))  %>%
      rename_with(.cols = 1, ~"ISTD_name") %>%
      cbind(., feature_table[.$index.B,])

    #ISTD intensities only
    istd_int<-t(istd_res[,-c(1:21)])

    #Loop to plot istd intensities in all samples
    pdf_name<-paste(pdf_loc, "1-", study_id, "_", corr_mode,"-Batch_Correction_ISTD-Intensities.pdf", sep="")
    pdf(pdf_name)
    for(kk in 1:nrow(istd_res)){
      print(
        ggplot(istd_int, aes(y=istd_int[,kk], x= 1:nrow(istd_int), fill= as.factor(mapfile$Batch))) +
          geom_bar(stat = "identity", colour="black", size=0.2) +
          scale_fill_manual(values = color_index) +
          labs(x="Sample number", y= "Intensity", size=20) +
          guides(fill=guide_legend(title="Batch Number")) +
          ggtitle(label= paste("Peak area for ",istd_res[kk,1] , sep=" ")) +
          theme(plot.title = element_text(hjust = 0.5))
      )
    }
    dev.off()
  }

}



#Function to perform batch correction using multiple methods
gc_batch_correct<-  function(feature_table= rc_ft,     #feature table with mzs in rows and samples in columns
                             cols_meta= 4,                    #number of columns meta data
                             corr_modes=c("MetaComBat", "WaveICA", "LIMMA"),  #Which batch correction modes
                             outloc=rcr_outloc,                 #output location
                             istds=istds,
                             mapfile=mapfile,
                             study_id,
                             istd_plots= FALSE
){

  print(paste( "Starting batch correction with ", paste(corr_modes, collapse=", "), sep=""))

  bc_res<-c()
  for(jj in 1:length(corr_modes)){


    ###################
    #####MetaCombat####
    if(corr_modes[jj]=="MetaComBat"){

      print("Starting batch correction with MetaCombat")
      #Prepare tables for combat
      combat_dat<-feature_table[,-(1:cols_meta)] %>% #Select intensities only
        `[<-`(. == 0, value = NA) %>%  #Replace zeroes with NA
        log(.,2) %>%
        data.frame(feature_table[,c(2,2:4)], .) %>%
        as.matrix()

      mapfile$Batch<-as.numeric(mapfile$Batch)

      #Run MetaCombat from xMSanalyzer
      mc_res <- xMSanalyzer::MetabComBat(dat=combat_dat, saminfo=mapfile[,1:3], type = "txt",
                                         write = F, covariates = "all", par.prior = T,
                                         filter = F, skip = 0, prior.plots = F)

      #Select intensities only; #Replace NAs and <0 values with zero
      mc_int<-mc_res[,-c(1:4)]


      #Replace NAs and <0 values with zero
      mc_int <- replace(as.matrix(mc_int),
                        which(is.na(mc_int) == TRUE), 0)

      mc_ft_out<-cbind(as.data.frame(feature_table[,1:cols_meta]), as.data.frame(mc_int))
      write.table(mc_ft_out, paste(outloc, "2-", study_id, "_Untargeted_MetaCombat_BatchCorrected_Feature_Table.txt", sep=""),
                  sep="\t", row.names=FALSE)

      bc_res$MetaComBat<- mc_ft_out


      print("Summarizing batches after MetaComBat")
      batch_eval_plots(feature_table= mc_ft_out,     #Feature table from data extraction
                       mapfile= mapfile,             #Mapfile
                       istd_df= istds,               #Internal standard mzs and rt
                       pdf_loc= outloc,           #Where to save pdfs
                       corr_mode= "MetaComBat",      #Batch correction mode
                       istd_plots= istd_plots,
                       study_id=study_id,
                       cols_meta= cols_meta
      )

      #####WaveICA1.0####
    } else if(corr_modes[jj]=="WaveICA") {

      print("Starting batch correction with WaveICA")

      #Prepare tables for WaveICA1.0
      #Transpose and log2 transform feature intensities
      wi_feature_table<-t(feature_table[,-c(1:cols_meta)]) %>%
        `[<-`(. == 0, value = 1) %>%  #Replace zeroes with 1
        log(.,2)
      colnames(wi_feature_table)<-feature_table[,1]

      #Run WaveICA1 using batch and sample type variables
      gc_Wave1<-WaveICA::WaveICA(data=wi_feature_table , wf="haar", batch= mapfile$Batch, group=mapfile$Sample_type, K=10, t=0.05, t2=0.05, alpha=0)
      waveICA_res<-gc_Wave1$data_wave


      #Write Wave ICA corrected table to results
      wICA_ft_out<-cbind(feature_table[,1:cols_meta], as.data.frame(t(waveICA_res)))
      write.table(wICA_ft_out, paste(outloc, "3-", study_id, "_Untargeted_WaveICA1.0_BatchCorrected_Feature_Table.txt", sep=""),
                  sep="\t", row.names=FALSE)

      bc_res$WaveICA1<- wICA_ft_out

      batch_eval_plots(feature_table= wICA_ft_out,     #Feature table from data extraction
                       mapfile= mapfile,               #Mapfile
                       istd_df= istds,                 #Internal standard mzs and rt
                       pdf_loc= outloc,                #Where to save pdfs
                       corr_mode= "WaveICA1.0",        #Batch correction mode
                       istd_plots= istd_plots,
                       study_id=study_id,
                       cols_meta= cols_meta
      )

      #####LIMMA####
    }   else if(corr_modes[jj]=="LIMMA") {

      print("Starting batch correction with LIMMA")

      limma_ft<-feature_table[,-c(1:cols_meta)] %>%
        `[<-`(. == 0, value = NA) %>%  #Replace zeroes with NA
        log(.,2)

      limma_res<-limma::removeBatchEffect(limma_ft, batch=mapfile$Batch, batch2=mapfile$Sample_type) %>%
        `[<-`(is.na(.), value = 0)


      #Write Limma corrected table to results
      limma_ft_out<-cbind(feature_table[,1:cols_meta], limma_res)
      write.table(limma_ft_out, paste(outloc, "3-", study_id, "_Untargeted_LIMMA_BatchCorrected_Feature_Table.txt", sep=""),
                  sep="\t", row.names=FALSE)

      bc_res$LIMMA<- limma_ft_out

      batch_eval_plots(feature_table= limma_ft_out,    #Feature table from data extraction
                       mapfile= mapfile,               #Mapfile
                       istd_df= istds,                 #Internal standard mzs and rt
                       pdf_loc= outloc,             #Where to save pdfs
                       corr_mode= "LIMMA",             #Batch correction mode
                       istd_plots= istd_plots,
                       study_id=study_id,
                       cols_meta= cols_meta
      )


    } else  {
      print("No batch correction method selected. Please select at least one of MetaComBat, WaveICA or LIMMA")
    }
  }

  return(bc_res)

} #Function end


