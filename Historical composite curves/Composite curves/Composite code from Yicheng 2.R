

#########################
detach("package:plyr", unload=TRUE)
###BINNING Z SCORES:
rm(list = ls())

setwd("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/ZTrans_regions/")

dat <- read.csv('MIDE_dataset.csv')

length(unique(dat$ID_ENTITY)) 
#Write entities to comma seaprated vector for mySQL
vec<- as.character(unique(dat$ID_ENTITY))
vec <- paste0(vec, collapse = ",")
vec
# Function that bins data:
library(dplyr)
library(tidyr)

Zbinning <- function(dt, bins, binhw){
  #Reshape data
  dat_res <- dt[,c("ID_ENTITY", "zt")] %>% group_by(ID_ENTITY) %>% mutate(number = row_number()) %>% spread(ID_ENTITY, zt) ##ive modified sarahs original code that worked on her data. if you want to see this look at the original PCA script.
  dat_res <- dat_res[,-1] ##takes away a 'number' column
  dat_res <- as.matrix(dat_res)
  
  Age_res <- dt[,c("ID_ENTITY", "EST_AGE")] %>% group_by(ID_ENTITY) %>% mutate(number = row_number()) %>% spread(ID_ENTITY, EST_AGE)
  Age_res <- Age_res[,-1]
  Age_res <- as.matrix(Age_res)
  
  #binning
  result <- matrix(ncol = length(Age_res[1,]), nrow = length(bins))
  colnames(result) <- colnames(dat_res) ##added this in so that column names of the binned matrix are entity IDs
  for (k in 1:length(dat_res[1,])){
    if(length(dat_res[is.na(dat_res[,k]) == F, k]) != 0){
      for (i in 1:length(bins)){
        t <- na.omit(cbind(as.numeric(Age_res[,k]), as.numeric(dat_res[,k])))
        result[i,k] <- mean(t[t[,1] > bins[i] - binhw & t[,1] < bins[i] + binhw, 2])
      }
    }
  }
  BinnedData <- structure(result, row.names = as.character(bins), col.names = colnames(dat_res), class = "matrix") #originally class = 'matrix'
  return(BinnedData)
}

Zbin_dat <- Zbinning(dt = dat, bins = seq(-70, 210,10), binhw = 10) #for now will go with  bin width of 250
colnames(Zbin_dat)

#give ages of each bin as a column

bin_age <- seq(-70,210,10)


##assume it starts at 0.ive used the sequence from the function so it must be.

Zbin_dat <- cbind(bin_age, Zbin_dat)
Zbin_dat <- as.data.frame(Zbin_dat) 

###################################################################################################

###################################################################################################


###Now make a composite curve 



# names
basename <- "zt60_200"
binname <- "bw10" 

# paths for output .csv files -- modify as appropriate
datapath <- "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Composite curves/"

#Load the locfit library and set parameters. The window width hw is, following convention, the half-window width (of the tricube weight function). The number of bootstrap resampling replications is given by nreps
library(locfit)


# locfit (half) window-width parameter
hw <- 10 # bandwidth (smoothing parameter)

# number of bootstrap samples/replications
nreps <- 1000

#Set the target ages for fitted values.

# target ages for fitted values
targbeg <- -60
targend <- 210
targstep <- 10 ## this is based on the previous binning procedure

# array sizes
maxrecs <- 2000
maxreps <- 1000

# plot output 
plotout <- "screen" # "screen" # "pdf"



##############
##CURVE GENERATION 
 
#########################

## Z-score run

curvecsvpath <- paste(datapath,"Regional_analyses/",sep="")
curvecsvpath

# if output folder does not exist, create it
#dir.create(file.path(datapath, paste(queryname,"_curves/",sep="")), showWarnings=FALSE)
curvefile <- paste(curvecsvpath,"MIDE_",basename,"_",binname,"_",as.character(hw),"_",as.character(nreps),'_',as.character(length(Zbin_dat)-1),'n', ".csv", sep="") #adjust to particular data you use now. .
print(curvefile)


#Code block to implement writing .pdf to file: 

# .pdf plot of bootstrap iterations
if (plotout == "pdf") { ## when you define PDFfile over here, need to make sure you name it with the nsites and what baseperiod overlap level you insisted on.
  pdffile <- paste(curvecsvpath,"MIDE_",basename,"_",binname,"_",as.character(hw),"_",as.character(nreps),'_',as.character(length(Zbin_dat)-1),'n','27_09_2021',".pdf", sep="")
  print(pdffile)
}

#Read the list of sites to be processed. Note that this is the site list may contain sites that MIDEter transforming and presampling may have no useful data. Those sites are ignored when the data are read in


# read the list of sites


ns <- as.numeric(colnames(Zbin_dat))
ns <- na.omit(ns) ##Z score data. note that this excludes the age column and only takes entity ids. so be careful when using the length of this in loop below.
length(ns)

maxrecs
# arrays for data and fitted values
age <- matrix(NA, ncol=length(ns), nrow=maxrecs)
influx <- matrix(NA, ncol=length(ns), nrow=maxrecs)
nsamples <- rep(0, maxrecs)
targage <- seq(targbeg,targend,targstep) #I assume targage is the set intervals which the function uses to make the curve
targage.df <- data.frame(x=targage)
lowage <- targage - hw; highage <- targage + hw
ntarg <- length(targage)
yfit <- matrix(NA, nrow=length(targage.df$x), ncol=maxreps) #the actual curve values (fitted)

# arrays for sample number and effective window span tracking
ndec <- matrix(0, ncol=ntarg, nrow=length(ns))
ndec_tot <- rep(0, ntarg)
xspan <- rep(0, ntarg)
ninwin <- matrix(0, ncol=ntarg, nrow=length(ns))
ninwin_tot <- rep(0, ntarg)


###Read data
# read and store the presample (binned) files as matrices of ages and influx values
ii <- 0

length(ns) #important to note that this vector is one less than the number of columns of the binned data datMIDErame. careful when looping that you dont miss any entities. I use the loop counter as the length of the binned datMIDErame instead of the length of ns.

for(i in 2:length(Zbin_dat)){ ##important that its length. not the actual character. NB that it skips the AGE COLUMN!!!!!!!!!!!!! which is why the ages might be better as last column of dframe.
  indata <- Zbin_dat[,c(i,1)] #1 is the age column. might be better last, and the loop starts at 1?
  indata <- na.omit(indata)  ##remove NaNs.
  nsamp <-  length(indata$bin_age) 
  if (nsamp > 0) {
    ii <- ii+1
    age[1:nsamp,ii] <- indata$bin_age # bins in >45n subset, age in original globe curve
    influx[1:nsamp,ii] <- indata[,1] # 
    nsamples[ii] <- nsamp
  }
  print(i)
}

nsites <- ii


# number of sites with data
nsites

#Trim the input data to an appropriate range given the target ages, and censor (set to NA) any tranformed influx values greater than 10

## trim samples to age range
influx[age >= targend+hw] <- NA
age[age >= targend+hw] <- NA #dont do it for targbeg because beg. is essentially the present and in previous scripts you trimmed down to less than 2020. I think

# censor abs(influx) values > 10
influx[abs(influx) >= 10] <- NA
age[abs(influx) >= 10] <- NA # the ages where influx is greater than 10. not sure how influx is in age matrix. I guess it just works element wise because they are the same lengths. 

#Find number of sites and samples contributing to fitted values

#The number of sites with samples (ndec_tot) and the number of samples (ninwin_tot) that contribute to each fitted value are calculated, along with the effective window width or span (xspan). This code is clunky, but parallels that in the Fortran versions.

# count number of sites that contributed to each fitted value
ptm <- proc.time()
for (i in 1:ntarg) {
  agemax <- -1e32; agemin <- 1e32
  for (j in 1:nsites) {
    for (k in 1:nsamples[j]) {
      if (!is.na(age[k,j])) {
        ii <- (age[k,j]-targage[1])/targstep + 1
        #print (c(i,j,k,ii))
        if (ii > 0 && ii <= ntarg) {ndec[j,ii] = 1}
        if (age[k,j] >= targage[i]-hw && age[k,j] <= targage[i]+hw) {
          ninwin[j,i] = 1
          if (agemax < age[k,j]) {agemax <- age[k,j]}
          if (agemin > age[k,j]) {agemin <- age[k,j]}
        }
      }
    }
  }
  ndec_tot[i] <- sum(ndec[,i])
  ninwin_tot[i] <- sum(ninwin[,i])
  xspan[i] <- agemax - agemin
}
proc.time() - ptm


head(cbind(targage,ndec_tot,xspan,ninwin_tot)) 



###########################################################################################################################

#Curve-fitting and bootstrapping
#First, the overall composite curve, i.e. without bootstrapping, determined and saved for plotting over the individual bootstrap curves later. Second, the nreps individual bootstrap samples are selected, and a composit curve fitted to each and saved.

#4.1 Composite curve
#The steps in getting this curve include:

#  Reshaping the data matrices (age and influx) into vectors
##Running locfit(), and
#Getting fitted values at the target ages
#The composite curve using all data is calculated as follows:

ptm <- proc.time()

# 1. reshape matrices into vectors 
x <- as.vector(age)
y <- as.vector(influx)
lfdata <- data.frame(x,y)
lfdata <- na.omit(lfdata)
x <- lfdata$x; y <- lfdata$y
print(y)
# 2. locfit
# initial fit, unresampled (i.e. all) data
loc01 <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, family="qrgauss")
summary(loc01)
#?locfit

# 3. get  fitted values
pred01 <- predict(loc01, newdata=targage.df, se.fit=TRUE)
loc01_fit <- data.frame(targage.df$x, pred01$fit)
fitname <- paste("locfit_",as.character(hw), sep="")
colnames(loc01_fit) <- c("age", fitname)
head(loc01_fit)


proc.time() - ptm


#Bootstrap-by-site confidence intervals
#The bootstrap confidence intervals are obtained by sampling with replacement sites (as opposed to individual samples), fitting a composite curve using locfit(), and saving these. 

#The 95-percent confidence intervals are then determined for each target age using the quantile() function.

# bootstrap samples
ptm <- proc.time()
# Step 1 -- Set up to plot individual replications - grey lines I assume
curvecsvpath


##check y limits by scaling max and min of fitted vals

max(loc01_fit$locfit_10)
min(loc01_fit$locfit_10)
nsites

if (plotout == "pdf") {dev.out()}
plot(x, y, xlab="Age (BP 1950)", ylab=fitname, xlim=c(-60, 200), xaxp  = c(0, 200, 10), ylim=c(-2,2), type="n") ##make sure xlim is right.
title(main = 'MIDE_scores, bp = -60-200BP: min 50% overlap, bins = 27, hw = 10, nboot = 1000', cex.main = 0.6)

#?set.seed
# Step 2 -- Do the bootstrap iterations, and plot each composite curve
set.seed(10) # do this to get the same sequence of random samples for each run

for (i in 1:nreps) { #nreps is 200
  print(i)
  randsitenum <- sample(seq(1:nsites), nsites, replace=TRUE)
 # print(head(randsitenum))
  x <- as.vector(age[,randsitenum])
  y <- as.vector(influx[,randsitenum])
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  locboot <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, maxit=20, family="qrgauss")
  predboot <- predict(locboot, newdata=targage.df, se.fit=TRUE)
  yfit[,i] <- predboot$fit
  # note plotting lines is slowww
  lines(targage.df$x, yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
  if (i %% 10 == 0) {print(i)}
}
warnings()
# Step 3 -- Plot the unresampled (initial) fit
fitname <- paste("Zlocfit_45N_",as.character(hw), sep="")
colnames(loc01_fit) <- c("age", fitname)
lines(loc01_fit[,1], loc01_fit[,2], lwd=2, col="red") ##keep all 45-55n America red; 55-65N america blue, 65-75 america orange; 45-55 eurasia purple, 55-65 eurasia green, 65-75 eurasia pink.

# Step 4 -- Find and add bootstrap CIs
yfit95 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
yfit05 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
lines(targage.df$x, yfit95, lwd=1, col="red")
lines(targage.df$x, yfit05, lwd=1, col="red")##keep all 45-55n America red; 55-65N america blue, 65-75 america orange; 45-55 eurasia purple, 55-65 eurasia green, 65-75 eurasia pink.




#Output
#The fitted curves are written out in the usual way.
ent <- unique(dat$ID_ENTITY)

curveout <- data.frame(cbind(targage.df$x, pred01$fit, yfit95, yfit05, ndec_tot, xspan, ninwin_tot))
colnames(curveout) <- c("age", "locfit", "cu95", "cl95", "nsites", "window", "nsamples")
curveout <- curveout[-c(28),] #remove additional row
curveout
curveout$Region <- 'MIDE'
write.csv(curveout, curvefile, row.names = F)
write.csv(ent, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Composite curves/Regional_analyses/MIDE_entities.csv")

#outputfile <- paste(curvecsvpath, curvefile, sep="")
#write.table(curveout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
proc.time() - ptm

curveout


warnings()


###################################################################################################
###############################   END ######################################################




