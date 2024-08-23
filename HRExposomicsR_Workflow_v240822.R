library(xcms)
library(MsExperiment)
library(BiocParallel)
library(dplyr)
library(RColorBrewer)
library(microbenchmark)
library(readxl)
library(ggplot2)
library(gridExtra)
library(xMSanalyzer)
library(writexl)
library(RAMClustR)
library(Spectra)
library(MsBackendMsp)
library(msentropy)
library(data.table)
library(WaveICA)
library(tibble)
library(limma)


#1. Output directory
#2. Study ID
#3. Study mzML location
#4. Mapfile
#5. Alkanes file
#6. Number of cores
#7. Perform RI based retention time correction prior to alignment (TRUE or FALSE)
#8. Missing value percent filtering prior to batch correction
#9. Remove RTs based on max and min RTs in alkane sheet
#10. Average alkane RTs
#11. Samples for alignment
#12. Standard based retention time correction
#13. 13C internal standard database
#14. Samples to use for alignment

args<- c("/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/10-Update_240608/4-Data_Extraction_240807/",
        "CLU0102",
        "/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/10-Update_240608/2-mzML/",
        "/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/10-Update_240608/3-Extraction_inputs/EXPANSE_Testing_Mapfile_240806.txt",
        "/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/10-Update_240608/3-Extraction_inputs/240806_Test_Batches_Alkanes.xlsx",
        6,
        TRUE,
        10,
        TRUE,
        TRUE,
        TRUE,
        TRUE,
        "/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/10-Update_240608/1-13C_Database/RI_Calculations/RI-Updated_Batch14_13C_ISTD_avg.xlsx",
        c("QAQC","NIST1958","Spike"))


#Read in variables from bash script
args<- commandArgs(trailingOnly = TRUE)
print(sessionInfo())
print(args)


#If using a windows computer, set to TRUE
win_comp<- FALSE

######
#Study ID
study_id<- args[2]


#mzML file location (must be centroided)
mzML_loc<- args[3]


#Sample mapfile . Mapfile columns are ordered as follows: Filename, Sample ID, Batch Number, Sample class. Filename must be the first column.
mapfile<- read.table(args[4], sep= "\t", header= TRUE)


#Variable for subset retention correction. Should be in fourth column in mapfile and only include biological samples
subset_samples<- args[14:length(args)]
print(subset_samples)


#Alkane mapfile
alkane_rt<- read_xlsx(args[5])


#Specify number of cores
c1<- as.numeric(args[6])


#Read in ISTD file. Set ISTD_plots to FALSE if no directory present
if(!is.na(args[13])){
  istd_plots<- TRUE
  istds<- as.data.frame(read_xlsx(args[13]))

}	  else	{
  istd_plots<- FALSE

}


#Register multicores. If using windows computer, use Snow
if(win_comp){
  register(bpstart(SnowParam(c1)))

}	else	{
  register(bpstart(MulticoreParam(c1)))

  }


#######---#######
#1. Create output directory
xcms_outloc<- paste(args[1], "/Stage1-GCHRMS_", study_id, "_XCMSv", packageVersion("xcms"), "/", sep="")
dir.create(xcms_outloc, showWarnings= FALSE)


#2. lists all mzML files in directory and check number of files in folder compared to sequence file. Stop running if there is a difference
file.pattern = ".mzML"

mzML_files_FULL<- list.files(mzML_loc, pattern= file.pattern, full.names= TRUE)
mzML_files_BASE<- mzML_files_FULL %>%
  basename() %>%
  gsub(, pattern=file.pattern, replacement= "")


print(paste("Number of files: ", length(mzML_files_BASE), sep=""))

#Error checking to determine if number files matches sequence file
if(nrow(mapfile) != length(mzML_files_BASE)){
  print("STOP! NUMBER OF FILES DOES NOT EQUAL SEQUENCE FILE")

  if(nrow(mapfile) > length(mzML_files_BASE)){
    print("DATAFILES MISSING FROM FOLDER")
    setdiff(mapfile[,1], mzML_files_BASE)

  }	else	{
    print("DATAFILES MISSING FROM SEQUENCE LIST")
    setdiff(mzML_files_BASE, mapfile[,1])

  }	} else {
    mapfile_arrange<- arrange(mapfile, File.Name)

    if(any(mapfile_arrange[,1] != sort(mzML_files_BASE))){
      print("STOP! MZML FILES AND MAPFILE DO NOT MATCH")
      misnamed_files<- data.frame(cbind(mapfile[which(mapfile_arrange[,1] != sort(mzML_files_BASE)),1], mzML_files_BASE[which(mapfile_arrange[,1] != sort(mzML_files_BASE))]))
      colnames(misnamed_files)<- c("Mapfile", "RAW_Filenames")
      print(misnamed_files)

    }	else {


      #Create filename list from sequence file
      seq_filenames<- paste(mzML_loc, "/", mapfile[,1], ".mzML", sep="")


      #Read raw data into MSExperiment file
      raw_data<- readMsExperiment(spectraFiles = seq_filenames)


      #Perform retenton time indix calculation
      if(args[7] == TRUE){
        print("Performing retention time indices correction")

        raw_data<-xcms_RT_to_RI_calc(ms_data = raw_data,
                                      alkanes = args[5],
                                      average_rt = args[10],
                                      filter_rt = args[9],
                                      mapfile = mapfile)

        merge_rt<-25

      }	else	{
        print("Retention time indices adjustment not performed. Using raw retention times.")
        merge_rt<-10

      }


        #XCMS data processing steps
        #Step 1 XCMS peak detection parameters
        xcms_params<-c()
        xcms_params$cwp_ppm= 5
        xcms_params$cwp_peakwidth= c(5,25)
        xcms_params$cwp_snthr= 50
        xcms_params$cwp_mzdiff= -0.001
        xcms_params$cwp_noise= 20000
        xcms_params$cwp_prefilter= c(5,20000)
        xcms_params$cwp_mzCenterFun= "wMean"
        xcms_params$cwp_integrate= 1
        xcms_params$cwp_fitgauss= FALSE
        xcms_params$cwp_extendLengthMSW=TRUE

        #Step 1a Merge neighboring peaks parameters
        xcms_params$mrg_expandMz = 0
        xcms_params$mrg_ppm = 3
        xcms_params$mrg_minProp = 0.9

        #Step 1b XCMS grouping 1 parameters
        xcms_params$grp1_minFraction = 0.05
        xcms_params$grp1_bw = 2
        xcms_params$grp1_ppm=10
        xcms_params$grp1_binSize=0.001

        #Step 2 XCMS retention time correction
        xcms_params$rtcor_minFraction = 0.75
        xcms_params$rtcor_extraPeaks = 1
        xcms_params$rtcor_smooth = "linear"
        xcms_params$rtcor_span = 0.2
        xcms_params$rtcor_subsetAdjust = "average"
        xcms_params$rtcor_family = "gaussian"

        #Step 3 XCMS grouping 2 parameters
        xcms_params$grp2_minFraction = 0.05
        xcms_params$grp2_bw =2
        xcms_params$grp2_ppm=10
        xcms_params$grp2_binSize=0.001


























    }  #Bracket at end of script
