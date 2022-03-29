
### Raw MS data pipeline for distinguishing bacterial subspecies - Biotyper

## Motivation
## To distinguish between bacterial sub types, Mass Spectrum (MS) data can be used. It allows identification of differences 
## in protein expression between bacterial samples treated and untreated with antibiotic. This pipeline allows processing 
## of raw MS data required for analysis and visualization of the results, which in turn can help distinguish between 
## the bacterial sub types based on their proteome.

## Libraries

library("MALDIquant")
library("MALDIrppa")
library("readxl")
library("plyr")
library("parallel")
library("tidyverse")
library("doParallel")
library("data.table")

## Data import

# Loading spectra

# Importing the spectrum data from .dat files produced by BRUKER MALDI software and the associated metadata, by loading 
## content of entire folder. Immediate conversion into MassSpectrum object. Creating a back up object. Printing summary 
## spectra (first 10).

spectra <- importSpectra(where = "PATH")
spectra2 <- spectra
summarySpectra(spectra[1:10])

# Loading metadata
# Reading in the associated metadata, directly from excel. Selecting for columns of interest only.

metadata <- read_excel("PATH.xlsx")

## Data checks
## Checking for compatibility between metadata and spectra.

# Name matching
# Checking if names in metadata and spectra match.

setdiff(metadata$"File name",names(spectra))
setdiff(names(spectra),metadata$"File name")

# Duplications
# Checking if there is any duplication in both metadata and spectra.

any(duplicated(metadata$"File name")) 
any(duplicated(names(spectra)))

# Correct format
# Checking if file names format from spectra matches file name format in metadata. 

all(names(spectra) == metadata$"File name")


# Setting all File names in spectra and metadata to be of character format. Performing additional check if format matches.

spectra <- spectra[as.character(metadata$"File name")]
p1 <- as.character(metadata$"File name")
p2 <- as.character(names(spectra))
all(p1 == p2) 

# Setting all remaining data in metadata object to factor format

metadata$Species <- as.factor(metadata$Species)
metadata$Number <- as.factor(metadata$Number)
metadata$Strain <- as.factor(metadata$Strain)
metadata$Antibiotic <- as.factor(metadata$Antibiotic)
metadata$Hrs <- as.factor(metadata$Hrs)
metadata$"Bio Rep" <- as.factor(metadata$"Bio Rep")
metadata$"Tech Rep" <- as.factor(metadata$"Tech Rep")
metadata$"Spot Rep" <- as.factor(metadata$"Spot Rep")

## Processing parameter settings
## Setting the optimal combination of parameters for peak profiles, to increase sensibility and decrease false discovery rate.

thScale <- 2.5 # Smoothing
ite <- 105 # Baseline correction
SigNoi <- 3 # Peak extraction
hws <- 20 # Peak extraction
tol <- 0.01 # Peak binning

## Initial screening
## Detection and filtering of low quality spectra. Creating a list of detected low quality spectra.

system.time(sc.results <- screenSpectra(spectra, meta = metadata))
summary(sc.results)
failure <- sc.results$est.table
failure <- failure[failure$Class == "failure",]

## Entire spectrum data plotting
## Plotting of the entire data and spectra detected by screenSpectra(), which will be retained (Intensity of base peak > 1000)

plot(sc.results, labels = TRUE)

## Updating results
## Updating files with filtered data.

spectra <- sc.results$fspectra 
metadata <- sc.results$fmeta 

## Denoising, baseline correction and normalization
## Performing all steps ensuring retention of spectrum ID's.

spectra.ID <- names(spectra)
spectra <- transfIntensity(spectra, fun = sqrt)
spectra <- wavSmoothing(spectra, method = "Wavelet", thresh.scale = thScale)
spectra <- removeBaseline(spectra, method = "SNIP", iterations = ite)
names(spectra) <- spectra.ID
spectra <- calibrateIntensity(spectra, method = "PQN")

## Range trimming
## Trimming spectra according to the data.

spectra <- trim(spectra, range = c(2500, 13000))

## Peak extraction and alignment 
## Detecting peaks, extracting and aligning them according to present parameters.

peaks <- detectPeaks(spectra, SNR = SigNoi, halfWindowSize = hws)
peaks <- alignPeaks(peaks, minFreq = 0.4, tolerance = tol)
summaryPeaks(peaks[1:10])

## Peak count and plotting
## Counting peaks for plotting of their pattern.

cP <- countPeaks(peaks)
plot(cP, type = "n")
text(cP, label = 1:length(cP))
peakPatterns(peaks)

## Elimination of rare peaks across the replicates
## Creating additional column containing information on Strain, Bio Rep, Hrs(time with antibiotic treatment) and Antibiotic Concentration. 

metadata$IsoBio <- as.factor(paste(metadata$Strain, metadata$"Bio Rep", sep = " "))
metadata$IsoBio <- as.factor(paste(metadata$IsoBio, metadata$Antibiotic, sep = " "))
metadata$IsoBio <- as.factor(paste(metadata$IsoBio, metadata$Hrs, sep = " "))
peaks.clean.f <- filterPeaks(peaks, minFreq = 0.4, labels = metadata$Strain)
cP <- countPeaks(peaks.clean.f)
plot(cP)
text(cP, label = 1:length(cP))
peakPatterns(peaks.clean.f)

## Outliers detection
## Detecting outliers and checking ratio to retained peaks. (optional step).

out <- detectOutliers(peaks.clean.f, metadata$"Strain", binary = TRUE)
table(out[,2])

## Discarding outlying peak profiles
## Discarding of the peaks detected in previous step.

peaks.clean.f <- peaks.clean.f[out[,2] == FALSE] 
metadata <- metadata[out[,2] == FALSE,] 

## Saving data 
## Storing intensity data as matrix (int) with NA's replaced by 0. Merging Intensity matrix with metadata.
## Saving data before aggregating at bio rep level (for entire data analysis).

int1 <- intensityMatrix(peaks.clean.f) 
int1[is.na(int1)] <- 0
colnames(int1) <- as.character(round(as.numeric(colnames(int1)),2))
dat1 <- cbind(metadata,int1)
save(dat1, file = "PATH.Rdata")

## Replicates merging for single composite representative peak profile by isolate
## Merge reps by the median at bio rep level for both intensity and metadata.

peaks.clean.BioRep <- mergeMassPeaks(peaks.clean.f, labels = metadata$IsoBio, method="median")
metadata.clean.BioRep <- aggregate(metadata, list(IsoBio = metadata$IsoBio), FUN = function(x) unique(x)[1])[,1:10]

## Conversion of massPeaks object into a matrix of intensity
## Creating intensity matrix merged at Bio Rep level, replacing NA's with 0s.

int2 <- intensityMatrix(peaks.clean.BioRep)
int2[is.na(int2)] <- 0
colnames(int2) <- as.character(round(as.numeric(colnames(int2)),2))
rownames(int2) <- as.character(metadata.clean.BioRep$IsoBio)

## Combining into a single data set
## Merging intensities with metadata after aggregation according to bio rep level.

dat2 <- cbind(metadata.clean.BioRep, int2)

## Saving processed data

save(dat2, file = "PATH.Rdata")
writeIntensity(int2, file = "PATH", format = "R", labels = rownames(int2))
save(metadata.clean.BioRep, file = "PATH.Rdata")





