rm(list = ls())

library(plyr)
setwd("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/TENA_W/")
data_query2 <- ldply(list.files(), read.csv, header=TRUE)


data_query <- data_query2[,c("ID_ENTITY", 'ID_SAMPLE',     "depth"   ,      "est_age"   ,   "influx" ,  "xst_level" )] 
head(data_query)
names(data_query) <- c("ID_ENTITY"  ,'ID_SAMPLE', "sample_depth"  ,     "mean"  ,   "quantity"  ,"TYPE" ) #sample ID
head(data_query)

ls3<- unique(data_query$ID_ENTITY)
paste0(ls3, collapse=",")



library(RMySQL)
mydb = dbConnect(MySQL(), user='root', password='Vedde12171', dbname='RPD snapshot 29.7.21', host='localhost')
dbListTables(mydb)


site <-  dbGetQuery(mydb, "select e.* from entity e 
left join sample s on e.ID_ENTITY = s.ID_ENTITY
left join age_model am on am.ID_SAMPLE = s.ID_SAMPLE
left join chronology c on c.ID_SAMPLE = s.ID_SAMPLE
where e.ID_ENTITY in (1,4,5,6,7,9,13,19,21,23,25,37,38,41,42,44,48,54,55,57,58,67,68,129,137,138,139,141,142,145,163,192,193,194,195,196,253,254,319,321,347,361,425,429,544,791,792,822,823,824,839,840,841,842,845,846,850,851,852,856,937,1423,1427,1429,1441,1442,1443,1445,1465,2025,2026,2062,2063,2064,2065,2067,2098,2157,2232,2255,2256)
                    group by e.ID_ENTITY;") #longer site dataset is because of those from WP that have no entities.




########################################## CALC COMPOSITE CURVE AND PLOT ###################################################################


#########################
detach("package:plyr", unload=TRUE)

###############################################################################################################################
##########################################################################################################################################################

###################################################
#MINIMAX, BOX COX AND Z SCORE TRANSFORMATIONS
###################################################

##purpose: Transform the obtained influx data to minimax, box cox and z scores
##input data: dat_infl: one dataframe of all entities in long format with influx-transformed quantities (output of influx calculation section).

##############
# 1-parameter Box-Cox transformation of charcoal quantities for each entity
# (alpha (shift parameter) is specified, lambda (power transformation parameter) is estimated)

# input .csv files should contain at least the variables "est_age" and "quant",
#   identified by labels in a header row  
###############

rm(list =  ls()[!(ls() %in% c('data_query', 'dat', 'site', 'ls3'))])#should probably save other files you generated
# paths for input and output .csv files -- modify as appropriate
length(unique(data_query$ID_ENTITY))
datapath <- "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion"


#define baseperiod for z score calculations
basebeg <- -60
baseend <- 200 
basename <- "zt60_200"
date <- '20_09_2021'

##calculation - initial steps. create folders and files

# various path and filenames

sitecsvpath <- paste("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/TENA_W/",sep="") #where all the csvs that have the influx data are stored.
sitecsvpath
transcsvpath <- paste("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/ZTrans_TENA_W/",sep="") # 
transcsvpath

# if output folder does not exist, create it
queryname <- 'TENA_W'
dir.create(file.path(datapath, paste(queryname,"_trans_csv/",sep="")), showWarnings=FALSE)

statscsvpath <- file.path(datapath, paste(queryname,"_trans_csv/",sep="")) #to store the stats
# if output folder does not exist, create it
statsfile <- paste(statscsvpath,basename,"_stats_age_", date,".csv", sep="") #JUST 'stats' for the first run for mnx and box cox.
statsfile
statscsvpath
# read list of sites

ptm <- proc.time()
sites <- site[site$ID_ENTITY %in% data_query$ID_ENTITY,]
sites <- sites[sites$ID_ENTITY %in% ls3,] #this deletes all obs for some reason revise next round
ent_vec <-unique(sites$ID_ENTITY)
nsites <- length(ent_vec)
print (nsites)


#sites <- sites[!sites$ID_ENTITY == c(85,941,1001),]

##Main loop
#Loop over the individual sites, doing the following:

#  compose the site .csv file name
#read the input data
#discard samples with mssing ages
#discard samples with ages after 2020 CE
#initial minimax rescaling (for historical reasons) (minimax)
#set alpha the Box-Cox transformation shift parameter
#maximum likelihood estimation of lambda
#Box-Cox transformation of data
#minimax rescaling tall
#calculate mean and standard deviation of data over base period
#calculate z-scores ztrans
#write out transformed data for this site

# THE BELOW LOOP IS JUST TO CHECK FOR ANY MISSING DATA IN INFLUX FILES BUT IS NOT NECESSARY##leave out all the sites that get an error due to having 0 influx, depths etc. Its a good test for the file reading in before the main loop.----
print(ent_vec)
ls1 <- c()
#Transform loop ----
for (z in 1:nsites[]) {  #Do not change this counter to the name. very dangerous becuase the matrices you use are indexed by this counter as well.
  
  # 1. site .csv file name (input file)
  ent_n<- ent_vec[z]
  ent_nm <- sites[sites$ID_ENTITY == ent_n,"entity_name"]
  entnmsub <- substr(ent_nm, start = 1, stop = 3)
  
  entidchar <- as.character(ent_n)
  if (ent_n >= 1) entid <- paste("000", entidchar, entnmsub, sep="") 
  if (ent_n >= 10) entid <- paste("00", entidchar, entnmsub, sep="")
  if (ent_n >= 100) entid <- paste("0", entidchar, entnmsub, sep="")
  if (ent_n >= 1000) entid <- paste(    entidchar, entnmsub, sep="")
  
  
  inputfile <- paste(sitecsvpath, entid,  "_RPD_", ".csv", sep="")
  print(z);print(entid); print(inputfile) #printing names etc to see. - you just use the inputfile one later.
  sitedata <- read.csv(inputfile)
  if(sum(sitedata$influx) == 0){ls1 <- append(ls1,ent_vec[z])} 
}

ls1  
check <- data_query[data_query$ID_ENTITY %in% ls1,]
data_query <- data_query[!data_query$ID_ENTITY %in% ls1,]
sites <- sites[!sites$ID_ENTITY %in% check$ID_ENTITY,]

##remove what you defined here because yuo want to define them again in next loop

rm('ent_n', 'ent_nm', 'entnmsub','entidchar','entid','inputfile','sitedata')

##none with 0 for influx.-----

#####################################
#####################################

###main transforming loop - MAKE SURE THE FILES ARE WRITING TO THE CORRECT FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# storage for statistics (vector arrays to store lambda, likelihoods etc)
ent_vec[1000]
sn_save <- rep(0,nsites); lam_save <- rep(0,nsites); lik_save <- rep(0,nsites)
tmean_save <- rep(0,nsites); tstdev_save <- rep(0,nsites)

ent_vec

##for rough global loop (unfiltered) where I havent cleaned up the dataset: just removing all those with only 0 for quantity
lsq <- c()
for(p in 1:length(ent_vec)){
  checksub <- data_query[data_query$ID_ENTITY == ent_vec[p],]
  if(length(unique(checksub$quantity))<2){lsq <- append(lsq, ent_vec[p])}
}
lsq
# ent_vec <- ent_vec[!ent_vec %in% lsq]
# nsites <- length(ent_vec)

for (j in 1:nsites) {  #looping through each .csv file of each site.
  ent_n<- ent_vec[j]
  ent_nm <- sites[sites$ID_ENTITY == ent_n,"entity_name"]
  entnmsub <- substr(ent_nm, start = 1, stop = 3)
  
  entidchar <- as.character(ent_n)
  if (ent_n >= 1) entid <- paste("000", entidchar, entnmsub, sep="") 
  if (ent_n >= 10) entid <- paste("00", entidchar, entnmsub, sep="")
  if (ent_n >= 100) entid <- paste("0", entidchar, entnmsub, sep="")
  if (ent_n >= 1000) entid <- paste(    entidchar, entnmsub, sep="")
  
  inputfile <- paste(sitecsvpath, entid,  "_RPD_", ".csv", sep="")
  print(j);print(entid); print(inputfile) #printing names etc to see. - you just use the inputfile one later.
  sitedata <- read.csv(inputfile)
  nsamp <- length(sitedata[,"ID_SAMPLE"])
  nsampchar <- as.character(nsamp)
  
  
  # 3. discard samples with missing (-9999) ages
  sitedata <- sitedata[sitedata$est_age != -9999,]
  
  # 4. discard samples with ages > -70
  sitedata <- sitedata[sitedata$est_age > -70,]
  
  if(nrow(sitedata) >1){ #added this in to account for those entities that dont have any original est ages and are stopping the loop running.
    
    # 5. initial minimax rescaling of data
    minimax <- (sitedata$influx-min(sitedata$influx))/(max(sitedata$influx)-min(sitedata$influx))
    
    # 6. set `alpha` the Box-Cox transformation shift parameter
    alpha <- 0.01  # Box-Cox shift parameter
    # alternative alpha: 0.5 times the smallest nonzero value of influx
    # alpha <- 0.5*min(sitedata$influx[sitedata$influx != 0])  
    
    # 7. maximum likelihood estimation of lambda
    # derived from the boxcox.R function in the Venables and Ripley MASS library included in R 2.6.1
    
    npts <- 201 # number of estimates of lambda
    y <- minimax+alpha
    n <- length(y)
    logy <- log(y)
    ydot <- exp(mean(logy))
    lasave <- matrix(1:npts)
    liksave <- matrix(1:npts)
    for (i in 1:npts) {
      la <- -2.0+(i-1)*(4/(npts-1))
      if (la != 0.0) yt <- (y^la-1)/la else yt <- logy*(1+(la*logy)/2*(1+(la*logy)/3*(1+(la*logy)/4)))
      zt <- yt/ydot^(la-1)
      loglik <- -n/2*log(sum((zt - mean(zt))^2 ))
      lasave[i] <- la
      liksave[i] <- loglik
    }
    
    # save the maximum likelihood value and the associated lambda
    maxlh <- liksave[which.max(liksave)]
    lafit <- lasave[which.max(liksave)]
    print (c(entid, maxlh, lafit))
    
    # 8. Box-Cox transformation of data
    if (lafit == 0.0) tall <- log(y) else tall <- (y^lafit - 1)/lafit #because you loop over 700 something sites, need this if else statement even if theres only one value of lambda and ml for one iteration. #code is essentially saying that if lambda is 0 you dont use it. just use log y. avoids mathematical problems with 0 exponent.
    
    
    # 9. minimax rescaling
    tall <- (tall - min(tall))/(max(tall)-min(tall))
    
    # 10. calculate mean and standard deviation of data over base period
    tmean <- mean(tall[sitedata$est_age >= basebeg & sitedata$est_age <= baseend])
    tstdev <- sd(tall[sitedata$est_age >= basebeg & sitedata$est_age <= baseend])
    
    # 11. calculate z-scores
    ztrans <- (tall-tmean)/tstdev
    
    # 12. write out transformed data for this site
    siteout <- data.frame(cbind(sitedata$ID_ENTITY, sitedata$ID_SAMPLE, sitedata$est_age, sitedata$depth, sitedata$influx, minimax, tall, ztrans))
    colnames(siteout) <- c( "ID_ENTITY", "ID_SAMPLE", "EST_AGE", "DEPTH", "INFLUX", "influxmnx", "tall", "zt")
    
    outputfile <- paste(transcsvpath, entid, "_Ztrans_", basename, "_RPD_",".csv", sep="") 
    write.table(siteout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
    transcsvpath
    sn_save[j] <- entid
    lam_save[j] <- lafit
    lik_save[j] <- maxlh
    tmean_save[j] <- tmean
    tstdev_save[j] <- tstdev
    
  }
}

tail(ent_vec)
#that ran fine. I think, but would like bart to quality control this. 

head(siteout)
head(outputfile)

outdata[outdata$sed_rate <=0,] ##seems to have purged all problems with negative sed rates/ influxes.

# write out a file of statistics
stats <- data.frame(cbind(sn_save, lam_save, lik_save, tmean_save, tstdev_save))
colnames(stats) <- c("site", "lambda", "likelihood", "mean", "stdev")
write.table(stats, statsfile, col.names=TRUE, row.names=FALSE, sep=",")

head(statsfile)

proc.time() - ptm


write.csv(ent_vec, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/TENA_W/final_entities.csv")




library(plyr)
setwd("/Users/paullincoln/Dropbox/2021/Research/RPD time series/ZTrans_LGM/")
LGM_dataset <- ldply(list.files(), read.csv, header=TRUE)
dir.create('/Users/paullincoln/Dropbox/2021/Research/RPD time series/Regional_transformed_data')

write.csv(LGM_dataset, './LGM_dataset.csv')

########################################## CALC COMPOSITE CURVE AND PLOT ###################################################################


#########################
detach("package:plyr", unload=TRUE)
###BINNING Z SCORES:
rm(list = ls())

setwd("/Users/paullincoln/Dropbox/2021/Research/RPD time series/Regional_transformed_data/")

dat <- read.csv('LGM_dataset.csv')

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

Zbin_dat <- Zbinning(dt = dat, bins = seq(-60, 200,10), binhw = 10) #for now will go with  bin width of 250
colnames(Zbin_dat)
??Zbinning
#give ages of each bin as a column

bin_age <- seq(-60,200,10)


##assume it starts at 0.ive used the sequence from the function so it must be.

Zbin_dat <- cbind(bin_age, Zbin_dat)
Zbin_dat <- as.data.frame(Zbin_dat) 

###################################################################################################

###################################################################################################


###Now make a composite curve 



# names
basename <- "zt60_200"
binname <- "bw200" 

# paths for output .csv files -- modify as appropriate
dir.create( "/Users/paullincoln/Dropbox/2021/Research/RPD time series/Regional_transformed_data/Composite curves/")
datapath <- "/Users/paullincoln/Dropbox/2021/Research/RPD time series/Regional_transformed_data/Composite curves/"

#Load the locfit library and set parameters. The window width hw is, following convention, the half-window width (of the tricube weight function). The number of bootstrap resampling replications is given by nreps
library(locfit)


# locfit (half) window-width parameter
hw <- 10 # bandwidth (smoothing parameter)

# number of bootstrap samples/replications
nreps <- 1000

#Set the target ages for fitted values.

# target ages for fitted values
targbeg <- 20000
targend <- 22000
targstep <- 200 ## this is based on the previous binning procedure

# array sizes
maxrecs <- 2000
maxreps <- 1000

# plot output 
plotout <- "screen" # "screen" # "pdf"



##############
##CURVE GENERATION 

#########################

## Z-score run

curvecsvpath <- paste(datapath ,sep="")
curvecsvpath

# if output folder does not exist, create it
#dir.create(file.path(datapath, paste(queryname,"_curves/",sep="")), showWarnings=FALSE)
curvefile <- paste(curvecsvpath,"LGM_",basename,"_",binname,"_",as.character(hw),"_",as.character(nreps),'_',as.character(length(Zbin_dat)-1),'n', ".csv", sep="") #adjust to particular data you use now. .
curvefile

#Code block to implement writing .pdf to file: 

# .pdf plot of bootstrap iterations
if (plotout == "pdf") { ## when you define PDFfile over here, need to make sure you name it with the nsites and what baseperiod overlap level you insisted on.
  pdffile <- paste(curvecsvpath,"LGM_",basename,"_",binname,"_",as.character(hw),"_",as.character(nreps),'_',as.character(length(Zbin_dat)-1),'n','17_09_2021',".pdf", sep="")
  print(pdffile)
}

#Read the list of sites to be processed. Note that this is the site list may contain sites that LGMter transforming and presampling may have no useful data. Those sites are ignored when the data are read in


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

length(ns) #important to note that this vector is one less than the number of columns of the binned data datLGMrame. careful when looping that you dont miss any entities. I use the loop counter as the length of the binned datLGMrame instead of the length of ns.

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

if (plotout == "pdf") {pdf(file=pdffile)}
plot(x, y, xlab="Age (BP 1950)", ylab=fitname, xlim=c(-60, 200), xaxp  = c(0, 200, 10), ylim=c(-2,2), type="n") ##make sure xlim is right.
title(main = 'LGM_scores, bp = -60-200BP: min 50% overlap, bins = 27, hw = 10, nboot = 1000', cex.main = 0.6)

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


if (plotout == "pdf") {dev.off()}


#Output
#The fitted curves are written out in the usual way.

curveout <- data.frame(cbind(targage.df$x, pred01$fit, yfit95, yfit05, ndec_tot, xspan, ninwin_tot))
colnames(curveout) <- c("age", "locfit", "cu95", "cl95", "nsites", "window", "ninwin")
curveout
write.csv(curveout, curvefile, row.names = F)
#outputfile <- paste(curvecsvpath, curvefile, sep="")
#write.table(curveout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
proc.time() - ptm

curveout

warnings()


###################################################################################################
###############################   END ######################################################
