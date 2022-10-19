#R script from S8Assignment 

# Load packages
library(xcms)
library(RColorBrewer)
library(ggplot2)
library(magrittr)
library(plotly)
library(DT)

# XCMS Analysis 
## Load data
### List .mzML files contained in the working directory

path.mzML <- "C:/Users/Bernat/Desktop/UNI/3r 2010-2021/Omicas/Practicas/S8-20201119-practicum-20201119/Task" 

mzML.files <- list.files(path.mzML, pattern= ".mzML", 
                         recursive=T, full.names = T)


### Phenodata dataframe from the working directory structure

pheno <- phenoDataFromPaths(mzML.files)
pheno$sample_name <- rownames(pheno) 
pheno$sample_group <- pheno$class

#The phenoDataFromPaths function builds a data.frame representing the experimental design from the folder structure in which the files of the experiment are located.
#Each row is a serum sample.

### Reading data .mzML files with readMSData method from the MSnbase package.
raw_data <- readMSData(files = mzML.files,
                       pdata = new("NAnnotatedDataFrame", pheno),
                       mode = "onDisk")

#mode = "onDisk": reads only spectrum header from files, but no data which is retrieved on demand allowing to handle very large experiments with high memory demand.
#It does not upload it to RAM, only when I tell it to (on demand). It reads only the metadata, not the mass specs, which is what it weighs. 

#### **How much time did the entire chromatographic run span?**
#In order to see the entire chromatographic run time we plot TIC: 
  
## We define colors according to experimental groups
group_colors <- paste0(brewer.pal(length(levels(pData(raw_data)$sample_group)),
                                  "Dark2"), "60")
names(group_colors) <- levels(pData(raw_data)$sample_group)

## We plot TIC using chromatogram() function 
## raw_data is the s4 object
TICs <- chromatogram(raw_data, aggregationFun = "sum") 

# Plot of TIC according to experimental groups
plot(TICs, col = group_colors[raw_data$sample_group], main = "TIC")
#We note the entire chromatographic run takes 600 seconds, that is to say, 10 minutes. 



#### **Number of experimental groups and number of samples per group: ** 
table(pheno$sample_group)
#As it is shown, there are 6 experimental groups: CTR, DIA_0, DIA_18, PIO_0, PIO_18 and QC, with 14, 5, 6, 6, 6, and 8 samples respectively. 

## Step1: Chromatographic peak detection through centwave
CentWaveParam()
#it tell us the default parameters (already defined)


#We must modify these parameters according what it is described in the on-line xcms parameters
#We define xcms online parameters 
cw_onlinexcms_default <- CentWaveParam(ppm = 15, peakwidth = c(5, 20), 
                                       snthresh = 6, prefilter = c(3, 100),
                                       mzCenterFun = "wMean", integrate = 1L,
                                       mzdiff = 0.01, fitgauss = FALSE, 
                                       noise = 0, verboseColumns = FALSE, 
                                       roiList = list(), firstBaselineCheck = TRUE, 
                                       roiScales = numeric())


xdata <- findChromPeaks(raw_data, param = cw_onlinexcms_default) 


#findChromPeaks(): to detect chromatographic peaks. We will have an s4 object that stores the peaks it has found for each sample.
#The findChromPeaks() function will give us a result that will be a list of all the peaks it finds with the centWave function in the 45 samples we have in raw_data. 
#xData is the s4 object that contains all the xcms processing results. 
#When the findChromPeaks() finishes running, the results will be stored in the variable xData


#### **Number of detected peaks:**
df2DT <- chromPeaks(xdata)
dim (df2DT)

#We get the peaks it has found in the s4 object. 
#We have detected 819946 peaks.
xdata  

#Chromatographic peak detection:
  #Method we have used: centWave
  #We have identified 819946 peaks in 45 samples

### In order to get the Extracted Ion Chromatogram for L-methionine (retention time ~ 293 s):

## rt and m/z range of the L-methionine peak area 
M <- 149.051049291 #monoisotopic weight of L-methionine 
H <- 1.007276 #proton weight
MH <- M + H #theoretical weight = 150.0583

## L-methionine peak mz range width according to XCMS parameters
error <- cw_onlinexcms_default@ppm #errors in ppm 
mmu_min <- MH-(error*MH/1e6) #proton weight of L-methionine +- a ppm error: 15 ppm 
mmu_max <- MH+(error*MH/1e6)
mzr <- c(mmu_min, mmu_max) #m/z range

## L-methionine peak rt range width according to XCMS parameters
apex.Lmethionine <- 293 #retention time ~ 293 s
rt.window <- cw_onlinexcms_default@peakwidth[2] #max rt width in seconds
rtr <- c(apex.Lmethionine-rt.window, apex.Lmethionine+rt.window)


#### **Extracted Ion Chromatograph (EIC) of L-methionine from raw data** 
chr_Lmethionine  <- chromatogram(raw_data, mz = mzr, rt = rtr)

plot(chr_Lmethionine, col = group_colors[chr_Lmethionine$sample_group], lwd = 2)

#This is the Extracted Ion Chromatogram of L-methionine.
#It plots the chromatogram in a mass range and in a retention time range.
#We note that we have 45 peaks, one for each sample, and colored according to their experimental group. 
#We are in a mass range that goes from a minimum mass range (mmu_min) of 150,0561 to a maximum mass range (mmu_max) of 150,0606. 
#The retention time is 293 s, because when we know that a metabolite is present in the samples, we know which is its the retention time.
#We observe that the maximum intensity at the apex (maxo) is 80000 approximately. 

### **Is there methionine in the raw data?**

#### Plot of Lmethionine MS1 spectrum
spMS1.Lmethionine<- raw_data %>%
  filterRt(rt = c(apex.Lmethionine, apex.Lmethionine+1)) %>%
  filterFile(1) %>%
  spectra

ggplotly(plot(spMS1.Lmethionine[[1]]))

#Yes, we do have L-methionine in the raw data, because we see a peak with an m/z of 150,05732 which is the L-methionine weight we have calculated.  
#We had MH = 150,0583 which was the theoretical monoisotopic weight of L-methionine plus the proton weight, so we have a difference of 0,00098 ppm between the he theoretical weight of L-methionine and the weight of L-methionine we have calculated. 
#We can also see peaks with higher weights than L-methionine which means that when we are doing MS1 scan, we measure in a particular retention time all the fragments we have with their adducts, in in-source fragmentations, etc. We have an endless amount of molecules represented in a mass spectrum. 

## Step 2: Alignment 
#There can be a little bit of movement of the RT from sample to sample.
#These minimums alignments can be corrected and we do use adjustRtime() wrapper function. 
#We have different algorithms to adjust RT. We are using obiwarp. 
#ObiwarpParam() warps the (full) data to a reference sample.

#Settings for the alignment:
  #We are going to use the online xcms parameters: profStep= 1 --> binSize = 1 
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 1))
#It uses sample number 23 as a sample against which you align the remaining ones. 
#It is done for all the samples. 


## Step 3: Correspondence 
#Correspondence group peaks or signals from the same ion across samples. 

### Defining peak density parameters
#Setting on-line xcms default parameters for peak density using PeakDensityParam().
#We define the online xcms peak density parameters: onlinexcms_pdp (go to online xcms and we use the parameters we see in the alignment section).
onlinexcms_pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                                   bw = 5,
                                   minFraction = 0.5,
                                   binSize = 0.015,
                                   minSamples = 1,
                                   maxFeatures = 100)


#Performance of the correspondence analysis using optimized peak density settings on on-line xcms.
#In order to see how correspondence works we use groupChromPeaks().
xdata <- groupChromPeaks(xdata, param = onlinexcms_pdp) 
## xdata contains my features. 


### Correspondence Results 
#Extracting the feature definitions as data frame.
#We can see the correspondence results with featureDefinitions().
feature.info <- as.data.frame(featureDefinitions(xdata))

#### **Number of detected features:**
dim(feature.info)
#We have found 25417 features.

xdata  

# Correspondence:
  # 25417 features have been identified 


### Checking L-methionine features
#I look for what is the index in this feature matrix, where is the feature that represents the L-methionine:
ix_Lmethionine_features <- which(feature.info$mzmed >= mzr[1] &
                                   feature.info$mzmed <= mzr[2] & 
                                   feature.info$rtmed >= rtr[1] &
                                   feature.info$rtmed <= rtr[2])


ix_Lmethionine_features
#It returns 803, so we can say that the feature that represents the L-methionine is the feature 803.


### In order to see how many peaks are associated with L-methionine and if L-methionine was detected in all samples:
Lmethionine.features <- feature.info[ix_Lmethionine_features,]

datatable(Lmethionine.features[,-ncol(Lmethionine.features)]) %>%
  formatRound(c("mzmed", "mzmin", "mzmax"), 4) %>%
  formatRound(c("rtmed", "rtmax", "rtmin"),0) 

#Lmethionine.features finds the feature with the index 803: the feature FT00803. This feature has 45 peaks (npeaks): 14 peaks in the CTR group, 5 peaks in the DIA_0 group, 6 peaks in the DIA_18 group, 6 peaks in the PIO_0 group, 6 peaks in the PIO_18 group, and 8 peaks in the QC group. 
#So, xcms has found a feature that is L-methionine and is well aligned. 

#Is there more than one peak in features associated with methionine?
  #We can see there are 45 peaks associated with L-methionine.

#Was methionine detected in all samples?**
  #We can see that L-methionine was detected in all samples because this feature has 45 peaks (npeaks): 14 peaks in the CTR group, 5 peaks in the DIA_0 group, 6 peaks in the DIA_18 group, 6 peaks in the PIO_0 group, 6 peaks in the PIO_18 group, and 8 peaks in the QC group. This means this feature is in all the samples and therefore L-methionine. 

## Step 4: Missing Values 

#The feature matrix we will create later (fmat) contains NA values for samples in which no chromatographic peak was detected in the feature's m/z-rt region. We must fill in missing peaks in order to fill in the empty areas, and reduce the number of NA values in our matrix. 
## We get missing values before filling in peaks
apply(featureValues(xdata, filled = FALSE), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))


#We use the fillChromPeaks method to fill in intensity data for such missing values from the original files. 
xdata <- fillChromPeaks(xdata)


### **Plot methionine differences (in terms of feature intensities) according to experimental groups (boxplot or scatter plot)**

## We use featureValues to extract the integrated peak intensity per feature/sample
## We get feature intensity matrix and its dimension
fmat <- featureValues(xdata, value = "maxo", method = "maxint") 

dim(fmat)

#The featureValues method returns a matrix with rows being features and columns samples.
#We have 25417 features and 45 samples 


#### In order to make a boxplot with ggplot2, we must convert fmat (matrix) to data frame
df.fmat <- data.frame(fmat)


#### We select the row where our feature (FT00803) is located
myfeature.fmat <- df.fmat[803,]


#### We create onother dataframe that contains the samples, its intensity and its experimental groups of the feature FT00803
df.myInfo <- data.frame(FT00803 = t(myfeature.fmat), GROUP = pheno$sample_group)


#### Finally, we create the boxplot: **L-methionine differences (in terms of feature intensities) according to experimental groups**
ggplot(data = df.myInfo, mapping = aes(x=df.myInfo$GROUP, y= df.myInfo$FT00803, fill = group_colors[raw_data$sample_group] )) + 
  geom_boxplot(alpha = 0.5) +
  theme(legend.position="none")+
  scale_y_continuous(name = "Intensity") +
  scale_x_discrete(name = "Experimental group") +
  ggtitle("L-methionine differences (in terms of feature intensities) according to experimental groups")




#save.image("TaskSession8.Rdata")


























