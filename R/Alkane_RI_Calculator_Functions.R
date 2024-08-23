#Functions related to calculating RI from retention time indices

#Variables include:
#ms_data: MSExperiment object from XCMS
#alkanes: sheet containing alkane retention times for all batches. See example
#average_rt: If TRUE, duplicate alkane runs in a batch will be averaged. If FALSE, will use alkane rts provides
#filter_rt: If TRUE, remove rt before earliest alkane and after latest alkane. If FALSE, no filtering occurs.
#mapfile: Mapfile containing Filename, Sample_ID, Batch_Number, Sample_class. Filename must be the first column.

#1. Function to calculate RI from alkane RTs
xcms_RT_to_RI_calc<-function(	ms_data= raw_data,    #MSExperiment object from XCMS
                              alkanes= "/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/10-Update_240608/3-Extraction_inputs/240806_Test_Batches_Alkanes.xlsx",
                              average_rt= TRUE,
                              filter_rt= TRUE,
                              mapfile= mapfile) {

  #Loop to adjust retention time with retention index
  #Results vector
  RI_rt_FINAL<- c() 	#Save adjusted RTs
  batch_n<- 0			#Batch column adjuster for alkanes_RT


  alkane_rt<-read_xlsx(alkanes)	#Read in alkane retention times

  #Remove any NA rows from alkane rt
  alkane_rt[alkane_rt == "NA"]<- NA
  alkane_rt<-alkane_rt[complete.cases(alkane_rt[, -c(1:5)]), ]

  alkane_rt<- alkane_rt %>%
    mutate_at(c(1, 6:ncol(alkane_rt)), ~ as.numeric(.))


  #If RTs average, first calculate average for each batch
  if(average_rt) {
    avg_rts<- c()
    for(kk in seq(6, ncol(alkane_rt), 2)) {
      avg_rts<- cbind(avg_rts, rowMeans(alkane_rt[, c(kk,kk+1)]))

    }
    alkane_rt<-cbind(alkane_rt[, 1:5], avg_rts)
  }

  #Determine min and max retention time
  min_time<- max(alkane_rt[1, -c(1:5)]) + 0.005
  max_time<- min(alkane_rt[nrow(alkane_rt), -c(1:5)]) - 0.005

  #Filter retention times to only include features with specified range
  if(filter_rt){
    ms_data<- filterRt(ms_data, rt = c(min_time * 60, max_time * 60))
  }

  #for loop to update retention times with retention index
  for(ii in unique(mapfile$Batch)) {

    #List of FileID numbers from sequence file
    batch_rownames<- as.numeric(rownames(mapfile[mapfile$Batch == ii, ]))

    #Select features detected in each file from corresponding batch
    raw_rt<- rtime(ms_data[batch_rownames])
    RI_rt_ii<- rep(0, length(raw_rt))

    #Loop through each pair of alkanes and calculate RI for each alkane pair
    for(jj in 2:nrow(alkane_rt)) {

      #Alkane retention times
      alkane_1<- as.numeric(alkane_rt[jj - 1, 6 + batch_n] * 60) #Tz
      alkane_2<- as.numeric(alkane_rt[jj, 6 + batch_n] * 60)    #Tz+1

      #Select all retention times between alkane_1 and alkane_2
      rt_jj<- which(raw_rt > alkane_1 & raw_rt <= alkane_2)

      #Calculate Retention time index
      RI_rt_ii[rt_jj]<- 100 * ((raw_rt[rt_jj] - alkane_1) / (alkane_2 - alkane_1) + as.numeric(alkane_rt[jj - 1, 1]))
    }
    batch_n<- batch_n + 1					#Advance to next batch RIs
    RI_rt_FINAL<-c(RI_rt_FINAL, RI_rt_ii)	#Final adjusted RTs
  }

  #Update raw retention time with alkane indices
  ms_data@spectra$rtime<- RI_rt_FINAL
  print("Retention time indices calculation complete")

  if(length(which(rtime(ms_data) == 0)) != 0){
    print("Error: Retention index of zero detected. Make sure at least one alkane retenton time exceeds max rt.")
  }

  return(ms_data)

} #End of function 1




########
#Variables include:
#input_file: File containing compounds of interest in each row and detected retention time (in minutes).
#time_col1: Time column in input file
#alkane_rt: File containing alkane retention times (in minutes)
#time_col2:  column to use in alkanes_rt
#outloc: Where to save files

#2. Function to calculate RI from list of chemicals with RT in minutes and alkane file
cmp_rt_to_ri<- function(input_file= "/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/~RI Converter/13C_ISTD_MSSM-Exploris.xlsx",
                        time_col1= "time",
                        alkane_rt= "/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/Alkane_RT_MSSM_Exploris.xlsx",
                        time_col2= "RT_mins_B1", #Time column to use in alkanes_rt
                        outloc= "/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/~RI Converter/") {


  #Loop to adjust retention time with retention index
  cmpds<- as.data.frame(readxl::read_xlsx(input_file))
  alkane_rt<- readxl::read_xlsx(alkane_rt)

  #Results vector
  RI_rt_FINAL<- c() 	#Save adjusted RTs

  #for loop to update retention times with retention index
  #Select raw rt
  raw_rt<- cmpds[, time_col1]
  RI_rt_ii<- rep(0, length(raw_rt))

  #Loop through each pair of alkanes and calculate RI for each alkane pair
  for(jj in 2:nrow(alkane_rt)) {

    #Alkane retention times
    alkane_1<- as.numeric(alkane_rt[jj - 1, time_col2]) #Tz
    alkane_2<- as.numeric(alkane_rt[jj, time_col2])    #Tz+1

    #Select all retention times between alkane_1 and alkane_2
    rt_jj<- which(raw_rt > alkane_1 & raw_rt <= alkane_2)

    #Calculate Retention time index
    RI_rt_ii[rt_jj]<- 100 * ((raw_rt[rt_jj] - alkane_1)/(alkane_2 - alkane_1) + as.numeric(alkane_rt[jj - 1, "Alkane_N"]))
  }


  #Update raw retention time with alkane indices
  cmpds$Alkane_RI<- RI_rt_ii
  out.names<-paste(outloc, "/", "RI-Updated_", basename(input_file), sep="")
  writexl::write_xlsx(cmpds, out.names)

} #End of function 2



########
#Variables include:
#input_file: File containing compounds of interest in each row and detected retention time (in minutes).
#time_col1: Time column in input file
#alkane_rt: File containing alkane retention times (in minutes)
#time_col2:  column to use in alkanes_rt
#outloc: Where to save files

#3. Function to convert RI back to RT in minutes


