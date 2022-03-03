

###################################################
### MAIN TRANSFORMING AND COMPOSITING SCRIPT 
###################################################
#purpose: this script takes raw RPD data, filters it to desired criteria, transforms it to influx, including first transforming back to concentration for entities in influx calculated using an outdated age model, calculates z scores, then generates the desired composite curve.

######SECTIONS##########
#SELECT DATA FROM RPD
##REDUNDANT CHECKS
##INFLUX BACK TRANSFORMATION TO CONCENTRATION
#INFLUX CALCULATION 
#MINIMAX, BOX COX AND Z SCORE TRANSFORMATIONS
#FILTERING DATA AND BINNING 
#COMPOSITE CURVE GENERATION

##########################################################################################################################################################

###########################################
###SELECT DATA FROM RPD----
###########################################

##purpose: use RMySQL package to extract RPD data, desired subset e.g. >45N. Its intended that this data will have been cleaned of depths in cm, negative quantities, depths <-70 etc all these types of issues. but still good to redundantly check this (at the transforming stage? is it there.)
rm(list = ls())
install.packages('RMySQL')
library(RMySQL)
mydb = dbConnect(MySQL(), user='root', password='Vedde12171', dbname='RPD snapshot 29.7.21', host='localhost')
dbListTables(mydb)


#start with entities. we want all available entities that are in the RPD #SUBSETS DB

RPDent <-  dbGetQuery(mydb, "SELECT * FROM entity WHERE ID_ENTITY in(88, 241, 644, 987, 990, 994, 2287, 2444);") #longer site dataset is because of those from WP that have no entities.
str(RPDent)
length(unique(RPDent$ID_SITE))
length(unique(RPDent$entity_name))


RPDsite <- dbGetQuery(mydb, "SELECT s.* FROM site s
                      left join entity e on s.ID_SITE = e.ID_SITE 
                      WHERE e.ID_ENTITY in(88, 241, 644, 987, 990, 994, 2287, 2444 );") 

length(unique(RPDsite$ID_SITE))
length(unique(RPDsite$site_name))
str(RPDsite)

##narrow down site table to all sites that are in the entity table.

RPDsite <- RPDsite[RPDsite$ID_SITE %in% RPDent$ID_SITE,] ##404 sites.

#################

##now get the sample data for all these sites. should preserve sample IDs so you can eventually bind bacon chron.

RPDsample <-  dbGetQuery(mydb, "SELECT * FROM sample;")
RPDsample <- RPDsample[RPDsample$ID_ENTITY %in% RPDent$ID_ENTITY,]

length(unique(RPDsample$ID_ENTITY))

RPDchron <- dbGetQuery(mydb, "SELECT * FROM chronology;")
RPDchron <- RPDchron[RPDchron$ID_SAMPLE %in% RPDsample$ID_SAMPLE,] #wrote this myself PL may not be needed

#now merge into a single dataframe

RPDdat <- merge(RPDsample, RPDchron, by = "ID_SAMPLE", all.x = T)
RPDdat <- merge(RPDdat, RPDent, by = "ID_ENTITY", all.x = T)


##clean up unecessary columns that are not needed in influx and transforming scripts. we need the ID site, ID entity and TYPE variable.
head(RPDdat)
RPDdat <- RPDdat[,c('ID_SITE', 'ID_ENTITY', "ID_SAMPLE", 'entity_name',"avg_depth",'charcoal_measurement','original_est_age','TYPE')]

write.csv(RPDdat, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/MIDEdat.csv", row.names = F)
write.csv(RPDsite, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/MIDEsite.csv", row.names = F)
write.csv(RPDent, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/MIDEent.csv", row.names = F)
write.csv(RPDchron, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/MIDEchron.csv", row.names = F)


# check <- RPDdat[is.na(RPDdat$TYPE),]
# length(unique(check$ID_ENTITY))
# 
# check <- RPDsite[!(RPDsite$ID_SITE %in% check$ID_SITE),]
# unique(check$ID_SITE)
# 
# length(unique(RPDdat$ID_SITE))
# length(unique(RPDsite$ID_SITE))
# RPDsite[!(RPDsite$ID_SITE%in% RPDdat$ID_SITE),]



##############################################################################################----


ent <- read.csv('/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/MIDEent.csv', na.strings = 'NULL')

dat <- read.csv('/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/MIDEdat.csv' , na.strings = 'NULL')

chron <- read.csv('/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/MIDEchron.csv' , na.strings = 'NULL')

dat2 <- merge(ent, dat, by = 'ID_ENTITY')
dat3 <- merge(dat2, chron, by = 'ID_SAMPLE')
head(dat3)
length(unique(dat3$ID_ENTITY))
head(dat)
head(chron)
dat <- dat3


##############################################################################################

################################################-----
##INFLUX BACK TRANSFORMATION TO CONCENTRATION
################################################
#purpose: Transform all entities that are in influx and have an outdated age model, back to concentration using accumulation rate calculated from the model that was used to get it to influx in the first place (the outdated model). This is so that the section that transforms to influx for the purposes of z score calculation, uses the updated model to calculate influx from concentration.
#this section thus does the following:

#calculate accumulation rates for all entities in influx, using the outdated age model
#uses acc rate to back-transform these entities to concentration.
#changes TYPE variable to concentration.
#ensures that the updated chronology and the concentration data for these entities will feed into the influx transforming script.


site <- ent

data_query <- dat
head(data_query); tail(data_query)

# Parse the query
# The main part of the script loops over the individual entities, and does various checks and calculations, including
# calculation of sedimentation rates and deposition times
# checking for age or depth reversals
# setting of various indicator variables and flags
# calculation of concentration from influx
# further checking for anomalies


# number of samples in the file
nsamp <- length(data_query$ID_ENTITY); nsamp
length(unique(data_query$ID_ENTITY))
##lets check what 'TYPE' codes we have
unique(data_query$TYPE.x) #if any of your sites have a wierd one, see if you can figure out if its one you can use. in other words, look at all data you will use with TYPE != CONC or influx and make sure they defo are not either of these by checking sample size != UNKN and unit makes sense given the flag?) So you don’t lose data). 


##another problem is quantities that are less than 0. so remove these too - ideally this will be corrected in the database
data_query[data_query$charcoal_measurement <0,]
data_query <- data_query[!(data_query$charcoal_measurement < 0),]

neg.infl <- data.frame() #to save the data that are deleted because of negative sed rates. USE A DATAFRAME
# get names to check if renaming is necessary
names(data_query)
maxsites <- data_query$ID_ENTITY #changed from maxsites <- length(unique(data_query$ID_ENTITY)
class(data_query$ID_ENTITY)

##
length(unique(maxsites))
ls1 <- c() #check for 0 conc everywhere
#define empty df that youll bind new data to.
newdf <-data.frame(ID_ENTITY = as.numeric(), depth = as.numeric(), est_age = as.numeric(), sed_rate = as.numeric(), quant = as.numeric(), conc = as.numeric(),influx = as.numeric(), xst_level = as.character(), conc_source = as.character(), influx_source = as.character(), TYPE = as.character())

maxsites <- as.vector(unique(maxsites))
sites <- as.vector(unique(data_query$entity_name.x))
rm('j')
for (j in sites) { ##trying looping using entity names: changed from j in 1:maxsites
  #j <- 5
  nsamp <- 0
  sitedata <- data_query[data_query$entity_name.x ==  j, ]
  sitedata <- sitedata[order(sitedata$avg_depth),]
  nsamp <- length(sitedata$ID_ENTITY)

  # Process the data for the j-th site
  # Define some local variables, and replace NA's with a missing values code
  
if (nsamp > 1) { ##I changed this from 0 to 1, because we have entities that are only one row in v4.
  
  # local variables
  depth <- sitedata$avg_depth; age <- sitedata$original_est_age.x; quant <- sitedata$charcoal_measurement; xst_level <- sitedata$TYPE.x; sampid <- sitedata$ID_SAMPLE

  # recode NA's to missing
  miss <- -9999.0
  depth[is.na(depth)] <- miss
  age[is.na(age)] <- miss
  quant[is.na(quant)] <- miss
  
  # define some new variables
  thickness <- rep(miss, nsamp); dep_time <- rep(miss, nsamp); sed_rate <- rep(miss, nsamp)
  unit_dep_time <- rep(miss, nsamp) 
  

  # sed rate and deposition time
  # first (top) sample
  if (depth[1] != miss && depth[2] != miss) {
    thickness[1] <- (depth[2] - depth[1])*100.0 # meters to cm (depth in m, influx and conc in cm) - we corrected those where depths were in cm in gcd, so this will be rigorous hopefully
    dep_time[1] <- age[2] - age[1]
    if (dep_time[1] > 0.0) sed_rate[1] <- thickness[1]/dep_time[1]
    if (sed_rate[1] != miss) unit_dep_time[1] <- 1/sed_rate[1] 
  }
  # samples 2 to nsamp-1
  for (i in 2:(nsamp-1)) {
    if (depth[1] != miss && depth[2] != miss) {
      thickness[i] <- (depth[i+1] - depth[i])*100.0 
      dep_time[i] <- ((age[i+1] + age[i])/2.0) - ((age[i] + age[i-1])/2.0)
      if (dep_time[i] > 0.0) sed_rate[i] <- thickness[i]/dep_time[i]
      if (sed_rate[i] != miss) unit_dep_time[i] <- 1/sed_rate[i]  
    }
  }
  # last (bottom) sample
  if (depth[nsamp-1] != miss  && depth[nsamp] != miss) {
    thickness[nsamp] <- thickness[nsamp-1] # replicate thickness
    dep_time[nsamp] <- age[nsamp] - age[nsamp-1]
    sed_rate[nsamp] <- sed_rate[nsamp-1] # replicate sed_rate
    unit_dep_time[nsamp] <- unit_dep_time[nsamp-1] 
  }
  
  # counts of missing values
  depth_count <- 0; age_count <- 0; quant_count <- 0; sed_rate_count <- 0; sed_rate_flag <- 1
  depth_count <- sum(depth != miss)
  age_count <- sum(age != miss)
  quant_count <- sum(quant != miss)
  sed_rate_count <- sum(sed_rate != miss)
  if (sed_rate_count != nsamp) sed_rateflag = 0
  
  # check for age or depth reversals, and zero or negative sed rates (in nonmissing data)
  depth_reversal <- 0; age_reversal <- 0; sed_rate_zeroneg <- 0         
  for (i in 2:nsamp) {
    if (age[i] != miss && age[i-1] != miss && age[i] <= age[i-1]) age_reversal=1
    if (depth[i] != miss && depth[i-1] != miss) {
      if (depth[i] <= depth[i-1]) depth_reversal=1
    } 
  }
  for (i in 2:nsamp) {
    if (sed_rate[i] != miss && sed_rate[i] <= 0.0) sed_rate_zeroneg=1 #not sure why he put this in a second loop..
  }
  
 # alternative quantities
  
  conc <- rep(miss, nsamp); influx <- rep(miss, nsamp)
  influx_source <- rep("none", nsamp) ; conc_source <- rep("none", nsamp)
  
  # select case based on xst_level
  if (xst_level[1] == "influx")          # adopt influx values as they are, calculate concentration
  {  
    influx <- quant
    influx_source <- "data"
      if (influx != miss && unit_dep_time != miss && sed_rate != 0.0) {
      conc <- influx * unit_dep_time
      conc_source <- "calculated from influx "
    } else {
      conc <- quant
      conc_source <- "copied from quant "
    }
    
  }
  
  else if (xst_level[1] == "concentration")     #adopt conc values as they are
  {
    conc <- quant
    conc_source <- "data"
    if (conc != miss && sed_rate != miss && sed_rate != 0.0) {
      influx <- rep(NA, length(quant))
      influx_source <- "none, initial transformation of outdated influx back to conc, which doesnt apply if data is in conc/conc-like already"
    } else {
      influx <- rep(NA, length(quant))
      influx_source <- "none, initial transformation of outdated influx back to conc, which doesnt apply if data is in conc/conc-like already"
    }  
    
  }
  
  else if (xst_level[1] == "pollen concentration")     # assume quantity is concentration like
  {
    conc <- quant
    conc_source <- "C0P0"
    if (sed_rate != miss && sed_rate != 0.0) {
      influx <- rep(NA, length(quant))
      influx_source <- "none, initial transformation of outdated influx back to conc, which doesnt apply if data is in conc/conc-like already"
    } else {
      influx <- rep(NA, length(quant))
      influx_source <- "none, initial transformation of outdated influx back to conc, which doesnt apply if data is in conc/conc-like already"
    }    
    
  }
  
  else if (xst_level[1] == "per unit weight")     # for now we are assuming that this is like concentration, but meeting on 19/03 to discuss this. (I changed this from "SOIL") - yes its concentration like according to sandy.
  {
    conc <- quant
    conc_source <- "data"
    if (conc != miss && sed_rate != miss && sed_rate != 0.0) {
      influx <- rep(NA, length(quant))
      influx_source <- "none, initial transformation of outdated influx back to conc, which doesnt apply if data is in conc/conc-like already"
    } else {
      influx <- rep(NA, length(quant))
      influx_source <- "none, initial transformation of outdated influx back to conc, which doesnt apply if data is in conc/conc-like already"
    }  
    
  }
  
  else if (xst_level[1] == "count")     
  {
    conc <- NA
    conc_source <- "none - unit count"
    influx <- NA
    influx_source <- "none - unit count"
    
  }
  
  else if (xst_level[1] == "raw count")     
  {
    conc <- NA
    conc_source <- "none - unit count"
    influx <- NA
    influx_source <- "none - unit count"
    
  }
  
   else 
  {
    conc <- NA
    conc_source <- "none - couldn't verify unit"
    influx <- NA
    influx_source <- "none - couldn't verify unit"
  }
}
  
  # check for conc == 0.0 everywhere
  if(!all(is.na(conc))){
  nzero <- 0
  nzero <- sum(conc != 0.0)
  if (nzero == 0) {
    ls1 <- append(ls1, entname)
	}
  }
  #now bind all the data together.
  if (nsamp > 1 && nzero > 0) { ##Quality control bit after z score calculation: I changed nsamp from >0 to >1. not sure about nzero. dont think it applies to nzero. because its still writing out a csv file for those entities with only one row. and you dont want that!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    entidvec <- sitedata$ID_ENTITY
    entnamevec <- sitedata$entity_name.x
    #define type
    new_type <- rep(NA, nsamp)
    if(all(is.na(influx)) & any(!is.na(conc))){
    	new_type <- xst_level
    }else if (any(!is.na(influx)) & any(!is.na(conc))){new_type <- rep('CONC', nsamp)}
    #********need to account for 'COUN' that ends up NA after this loop***********
    
    outdata <- data.frame(entnamevec , entidvec, sampid, depth, age, sed_rate, quant, conc, influx, xst_level, conc_source, influx_source, new_type) 
    names(outdata) <- c('entity_name',"ID_ENTITY", 'ID_SAMPLE',"depth", "est_age", "sed_rate", "quant", "conc",
                        "influx", "xst_level", "conc_source", "influx_source", 'TYPE' ) #might change est age to original est age... for clarity. Also write out sample IDs, so you can bind new bacon chrons.
    
    
    #delete any rows in this data that have zero and negative sed rate, and save them so you can see how many are removed.
    neg.infl <- append(neg.infl, outdata[outdata$sed_rate <= 0,])
    outdata <- outdata[!(outdata$sed_rate <= 0),]
    newdf <- rbind(newdf, outdata)
    
      } 
  print(j)
}

ls1
outdata[outdata$sed_rate <= 0,]
names(newdf)
check <- newdf[newdf$ID_ENTITY %in% 5 , c('entity_name','depth', 'est_age','conc', 'influx', 'TYPE','sed_rate')]
head(data)
check2 <- dat[data_query$ID_ENTITY %in% 4, c('entity_name','sample_depth', 'original_est_age','quantity','TYPE')]
plot(check$est_age, check$conc, type = 'l')
plot(check$est_age, check$influx, type = 'l')
plot(check2$original_est_age, check2$quantity, type = 'l')

##values are slightly off, perhaps this is because, like in the palaeofire package, the sed rate calculation uses different depth differencing.

a <-unique(dat[dat$latitude > 45 & dat$TYPE %in% c('INFL', 'CONC'),'entity_name'])
a <- substr(a, start = 1, stop =  4)
dup <- a[duplicated(a)]
a <- a[!a %in% dup]




##########################################################################################################################################################

##################################             
#INFLUX CALCULATION                  
##################################
#purpose: calculate influx - write out a file for each entity with influx values, then read it in to bind into one dataframe. if there are new chronologies, use these.
#input data: RPD data cleaned and filtered to one entity per site.

#before doing this, will ensure that the data is all in concentration and has the new bacon age model, which will be used to generate influx. 

rm(list= ls()[!(ls() %in% c('dat','site','newdf'))])

# query label and path to query and query name
dir.create("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/MIDE")
csvpath <- "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/MIDE/"


#preserve only concentrations which you just got

data_query <- newdf[,c("entity_name"   ,"ID_ENTITY", 'ID_SAMPLE',     "depth"   ,      "est_age"   ,   "conc" ,  "TYPE" )] 
head(data_query)
names(data_query) <- c("entity_name", "ID_ENTITY"  ,'ID_SAMPLE', "sample_depth"  ,     "original_est_age"  ,   "quantity"  ,"TYPE" ) #sample ID
head(data_query)

# Parse the query for influx calcs----
# The main part of the script loops over the individual entities, and does various checks and calculations, including
# calculation of sedimentation rates and deposition times using the most up-to-date chronology.
# checking for age or depth reversals
# setting of various indicator variables and flags
# calculation of alternative quanties (e.g. influx, given concentrations)
# further checking for anomalies
# writing out a .csv file for each site

# number of samples in the file
nsamp <- length(data_query$ID_ENTITY); nsamp
length(unique(data_query$ID_ENTITY))

##lets check what 'TYPE' codes we have
unique(data_query$TYPE)
#data_query.m <- data_query
check <- data_query[is.na(data_query$TYPE),]
c2 <- dat[dat$ID_ENTITY %in% check$ID_ENTITY,]
unique(c2$TYPE)

#this is just an addition on 06/05/2021 to generate data z scores for data below 45N: need the type in code format: INFL, C0P0, COUN, CONC,
data_query[data_query$TYPE %in% "concentration",'TYPE'] <- 'CONC'
data_query[data_query$TYPE %in% "pollen concentration",'TYPE'] <- 'C0P0'
data_query[data_query$TYPE %in% "influx",'TYPE'] <- 'INFL'
data_query[data_query$TYPE %in% c("raw count",'count'),'TYPE'] <- 'COUN'
unique(data_query$TYPE)
head(data_query)

coun <- data_query[data_query$TYPE %in% c('COUN', 'soil','other') | is.na(data_query$TYPE),] ##for <45N:

#they are all counts, so should just fix it in the loop that gets conc (marked where). for now can leave it.

##in general in this script: at this point you want to exclude from this whole loop any entities that cant possibly be transformed to influx (e.g. type = COUNT.) per unit weight, conc, influx C0P0 you want to keep. so purge any that dont fit this HERE. - probably redundant given above getting conc and filtering for that.
coun <- data_query[data_query$TYPE %in% c('COUN', 'SOIL','OTHE') | is.na(data_query$TYPE),] ##soil is for global, the NA is temporary, should be resolved when resolving issue in conc transforming loop.
length(unique(coun$ID_ENTITY))
length(unique(data_query$ID_ENTITY))
unique(coun$TYPE)
site[site$ID_ENTITY %in% coun$ID_ENTITY,]
data_query <- data_query[!data_query$ID_ENTITY %in% coun$ID_ENTITY,]
unique(data_query$TYPE)
#shouldnt I subset the entity data to this now or do I do this below.
##so this should have all the entities I actually want. 

##another problem is quantities that are less than 0. so remove these too - ideally this will be corrected in the database

neg.dat <- data_query[data_query$quantity <0 | data_query$sample_depth <0,] #dont just delete this, some are sensical. they just took the core from a certain depth.
data_query <- data_query[!data_query$ID_ENTITY %in% neg.dat$ID_ENTITY,]
#length(unique(data_query$ID_ENTITY))
##########################


check <- data_query[is.na(data_query$ID_ENTITY),]
ls1 <- c() #to check which have 0 influx everywhere
ls2 <- c()
ls3 <- c()
neg.infl <- data.frame() #to save the data that are deleted because of negative sed rates. USE A DATAFRAME - need to add the right columns.
# get names to check if renaming is necessary
names(data_query)
maxsites <- unique(data_query$ID_ENTITY) #changed from maxsites <- length(unique(data_query$ID_ENTITY)
class(data_query$ID_ENTITY)




##ALWAYs CHANGE FILENAMES AND CHECK LOOP BEFORE RUNNING  - should actually change dates in loops with paste(date) so it automatically updates
sites <- as.vector(unique(data_query$entity_name))

for (j in maxsites) { ##trying looping using entity names: changed from j in 1:maxsites
  nsamp <- 0
  sitedata <- data_query[data_query$ID_ENTITY == j, ]
  sitedata <- sitedata[order(sitedata$sample_depth),]
  nsamp <- length(sitedata$ID_ENTITY)

  # Process the data for the j-th site
  # Define some local variables, and replace NA's with a missing values code
  
if (nsamp > 1) { ##I changed this from 0 to 1, because we have entities that are only one row in v4.
  jchar <- as.character(j)
  nsampchar <- as.character(nsamp)
  
  # local variables
  depth <- sitedata$sample_depth; age <- sitedata$original_est_age; quant <- sitedata$quantity; xst_level <- sitedata$TYPE 
  
  # recode NA's to missing
  miss <- -9999.0
  depth[is.na(depth)] <- miss
  age[is.na(age)] <- miss
  quant[is.na(quant)] <- miss
  
  # define some new variables
  thickness <- rep(miss, nsamp); dep_time <- rep(miss, nsamp); sed_rate <- rep(miss, nsamp)
  unit_dep_time <- rep(miss, nsamp) #removed x1st level here because added it as vector above
  

  # sed rate and deposition time
  # first (top) sample
  if (depth[1] != miss && depth[2] != miss) {
    thickness[1] <- (depth[2] - depth[1])*100.0 # meters to cm (depth in m, influx and conc in cm) - we corrected those where depths were in cm in gcd, so this will be rigorous hopefully
    dep_time[1] <- age[2] - age[1]
    if (dep_time[1] > 0.0) sed_rate[1] <- thickness[1]/dep_time[1]
    if (sed_rate[1] != miss) unit_dep_time[1] <- 1.0/sed_rate[1] #the 1/sed rate - 1 is a placeholder because it will be *quant
  }
  # samples 2 to nsamp-1
  for (i in 2:(nsamp-1)) {
    if (depth[1] != miss && depth[2] != miss) {
      thickness[i] <- (depth[i+1] - depth[i])*100.0 
      dep_time[i] <- ((age[i+1] + age[i])/2.0) - ((age[i] + age[i-1])/2.0)
      if (dep_time[i] > 0.0) sed_rate[i] <- thickness[i]/dep_time[i]
      if (sed_rate[i] != miss) unit_dep_time[i] <- 1.0/sed_rate[i] 
    }
  }
  # last (bottom) sample
  if (depth[nsamp-1] != miss  && depth[nsamp] != miss) {
    thickness[nsamp] <- thickness[nsamp-1] # replicate thickness
    dep_time[nsamp] <- age[nsamp] - age[nsamp-1]
    sed_rate[nsamp] <- sed_rate[nsamp-1] # replicate sed_rate
    unit_dep_time[nsamp] <- unit_dep_time[nsamp-1]
  }
  
  # counts of missing values
  depth_count <- 0; age_count <- 0; quant_count <- 0; sed_rate_count <- 0; sed_rate_flag <- 1
  depth_count <- sum(depth != miss)
  age_count <- sum(age != miss)
  quant_count <- sum(quant != miss)
  sed_rate_count <- sum(sed_rate != miss)
  if (sed_rate_count != nsamp) sed_rateflag = 0
  
  # check for age or depth reversals, and zero or negative sed rates (in nonmissing data)
  depth_reversal <- 0; age_reversal <- 0; sed_rate_zeroneg <- 0         
  for (i in 2:nsamp) {
    if (age[i] != miss && age[i-1] != miss && age[i] <= age[i-1]) age_reversal=1
    if (depth[i] != miss && depth[i-1] != miss) {
      if (depth[i] <= depth[i-1]) depth_reversal=1
    } 
  }
  for (i in 2:nsamp) {
    if (sed_rate[i] != miss && sed_rate[i] <= 0.0) sed_rate_zeroneg=1
  }
  
  # alternative quantities
  
  conc <- rep(miss, nsamp); influx <- rep(miss, nsamp)
  influx_source <- rep("none", nsamp) ; conc_source <- rep("none", nsamp)
  
  # select case based on xst_level
  if (xst_level[1] == "INFL")          # adopt influx values as they are, calculate concentration - should be no influxes as input, but leaving in anyway
  {  
    influx <- quant
    influx_source <- "data"
    if (influx != miss && unit_dep_time != miss && sed_rate != 0.0) {
      conc <- influx * unit_dep_time
      conc_source <- "calculated from influx "
    } else {
      conc <- quant
      conc_source <- "copied from quant "
    }
   
  }
  
  else if (xst_level[1] == "CONC")     # calculate influx, adopt conc values as they are
  {
    conc <- quant
    conc_source <- "data"
    if (conc != miss && sed_rate != miss && sed_rate != 0.0) {
      influx <- quant * sed_rate
      influx_source <- "calculated from conc "
    } else {
      influx <- quant
      influx_source <- "copied from quant "
    }  
   
  }
  
  else if (xst_level[1] == "C0P0")     # assume quantity is concentration like
  {
    conc <- quant
    conc_source <- "C0P0"
    if (sed_rate != miss && sed_rate != 0.0) {
      influx <- quant * sed_rate
      influx_source <- "calculated from C0P0 (conc) "
    } else {
      influx <- quant
      influx_source <- "copied from quant "
    }    
    
  }
  
  else if (xst_level[1] == "per unit weight")     # for now we are assuming that this is like concentration, but meeting on 19/03 to discuss this. (I changed this from "SOIL") - yes its concentration like according to sandy.
  {
    conc <- quant
    conc_source <- "data"
    if (conc != miss && sed_rate != miss && sed_rate != 0.0) {
      influx <- quant * sed_rate
      influx_source <- "calculated from per unit weight "
    } else {
      influx <- quant
      influx_source <- "copied from quant "
    }  
    
  }
  
  else if (xst_level[1] == "COUN")     
  {
    conc <- NA
    conc_source <- "none - unit count"
    influx <- NA
    influx_source <- "none - unit count"
    
  }
  
   else 
  {
    conc <- NA
    conc_source <- "none - couldn't verify unit"
    influx <- NA
    influx_source <- "none - couldn't verify unit"
    
  }
}
  
  # check for influx == 0.0 everywhere
    
  if (sum(influx) == 0) {
    ls1 <- append(ls1, j)
  }
  
  ##and NA everywhere - this is nonsensical because you define it as NA above for OTHE and COUN
  
  
  if (unique(xst_level) != 'COUN' & all(is.na(influx))) {
    ls2 <- append(ls2, j)
  }
  
  # .csv out
  if (nsamp > 1 && sum(influx) > 0) { ##Quality control bit after z score calculation: I changed nsamp from >0 to >1. 
    
    # get siteid string
    entname <- unique(sitedata$entity_name)
    entname <- entname[!is.na(entname)] ##for the global loop
    entname <- substr(entname, start = 1, stop = 3)
    entidchar <- as.character(j)
    if (j >= 1) entid <- paste("000", entidchar, entname, sep="") 
    if (j >= 10) entid <- paste("00", entidchar, entname, sep="")
    if (j >= 100) entid <- paste("0", entidchar, entname, sep="")
    if (j >= 1000) entid <- paste(    entidchar, entname, sep="")
    
    
    # assemble output data and write it out
    # also save a list of entity IDs that get successfully written out, for quality control.
     
     ls3 <- append(ls3,j)
   
    
    entidvec <- sitedata$ID_ENTITY
    sampidvec <- sitedata$ID_SAMPLE
    entname <-  unique(sitedata$entity_name)
    
    outdata <- data.frame( entidvec, sampidvec, depth, age, sed_rate, quant, conc, influx, xst_level, conc_source, influx_source) 
    names(outdata) <- c("ID_ENTITY", "ID_SAMPLE", "depth", "est_age", "sed_rate", "quant", "conc",
                        "influx", "xst_level", "conc_source", "influx_source" )
    
   #delete any rows in this data that have zero and negative sed rate, and save them so you can see how many are removed.
    neg.infl <- append(neg.infl, outdata[outdata$sed_rate <= 0,])
    outdata <- outdata[!(outdata$sed_rate <= 0),]
    csvfile <- paste("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/MIDE/",entid,"_RPD_", ".csv", sep="")
    write.csv(outdata, csvfile, row.names=FALSE)
    
      } 
  print(j)
}





ls1
check <- grep('Gar', data_query[,'entity_name'])
check <- data_query[check,]
ls2
length(unique(ls3))
length(unique(data_query$ID_ENTITY)) #4 didnt write out.
neg.infl
csvfile
nrow(site[site$ID_ENTITY %in% ls3,])

print(site$ID_ENTITY)
print(ls3,)

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
date <- '16_09_2021'

##calculation - initial steps. create folders and files

# various path and filenames

sitecsvpath <- paste("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/MIDE/",sep="") #where all the csvs that have the influx data are stored.
sitecsvpath
transcsvpath <- paste("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/ZTrans_MIDE/",sep="") # 
transcsvpath

# if output folder does not exist, create it
queryname <- 'MIDE'
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

ls1 <- c()
for (z in 1:nsites) {  #Do not change this counter to the name. very dangerous becuase the matrices you use are indexed by this counter as well.
 

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
ent_vec
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




##########################################  END OF MAIN LOOPING BIT ###################################################################


##now make a single file with all the transformed data (influx is also there).
#sites <- site
ent_vec
dat2 <- data.frame()
for (y in 1:length(ent_vec)) { 
	#y <- 3
ent_n<- ent_vec[y]
ent_nm <- sites[sites$ID_ENTITY == ent_n,"entity_name"]
entnmsub <- substr(ent_nm, start = 1, stop = 3)

entidchar <- as.character(ent_n)
if (ent_n >= 1) entid <- paste("000", entidchar, entnmsub, sep="") 
if (ent_n >= 10) entid <- paste("00", entidchar, entnmsub, sep="")
if (ent_n >= 100) entid <- paste("0", entidchar, entnmsub, sep="")
if (ent_n >= 1000) entid <- paste(    entidchar, entnmsub, sep="")


inputfile <- paste(transcsvpath, entid, "_Ztrans_RPD", basename, "_RPD_", ".csv", sep="")
print(y);print(entid); print(inputfile) #printing names etc to see. - you just use the inputfile one later.
if (file.exists(inputfile)){ sitedata <- read.csv(inputfile)}
dat2 <- rbind(dat2, sitedata)
}

length(unique(dat2$ID_ENTITY))
##less than 560. why.
check <- data_query[!(data_query$ID_ENTITY %in% dat2$ID_ENTITY),]
unique(check$ID_ENTITY)

##so some entities did not have ages.

##just check if any depths are in cm

cdepth <- dat2[dat2$DEPTH > 50,] #azzano decimo

datapath

#for a global dataset:
#  ne3045 <- read.csv('C:\\Users\\xr824576\\OneDrive - University of Reading\\PhD\\Data_analyses\\RPD\\RPD_trans_zt200_12000k_30_45N_06_05_2021.csv')
#  ne45 <- read.csv('C:\\Users\\xr824576\\OneDrive - University of Reading\\PhD\\Data_analyses\\RPD\\RPD_trans_zt200_12000k_45N__bacon2020_age_12_04_2021.csv')
# glob <- rbind(ne3045, ne45, dat2)
#  length(unique(glob$ID_ENTITY))

 #write.csv(glob, paste(datapath, '\\RPD_trans_',basename ,'_global_' ,date, '.csv', sep = ''), row.names = F)

#Remember to go back and repeat the loop filtering the entity IDs down to those that fit the criteria you choose. whether its 75% overlap or whatever.- or you could just select the sites after transforming.

dat <- dat2 #redefine your data to the subset you have got now.


check <- dat[dat$ID_ENTITY == 5,]
plot(check$EST_AGE, check$INFLUX, type = 'l')
plot(check$EST_AGE, check$influxmnx, type = 'l')
plot(check$EST_AGE, check$tall, type = 'l')
plot(check$EST_AGE, check$zt, type = 'l')
summary(dat2$EST_AGE)

ch <- site[site$ID_ENTITY %in% dat$ID_ENTITY,]
ch <- ch[ch$longitude > (-50) & ch$longitude < 50,]

##########################################################################################################################################################



#########################

###BINNING Z SCORES:

length(unique(dat$ID_ENTITY)) 
head(dat)
# Function that bins data:
#install.packages('dplyr')
library(dplyr)
#install.packages('tidyr')
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

Zbin_dat <- Zbinning(dt = dat, bins = seq(100,12000,250), binhw = 125) #for now will go with  bin width of 250
colnames(Zbin_dat)
#warnings()
#give ages of each bin as a column

bin_age <- seq(100,12000,250)


##assume it starts at 0.ive used the sequence from the function so it must be.

Zbin_dat <- cbind(bin_age, Zbin_dat)

## Limit this to the sites in lsover. ie. those with your decided cutoff of overlap with the baseperiod (in this case, 50% on 18/03/2020)



class(lsover)
lsoverch <- as.character(lsover)
class(lsoverch)
class(colnames(Zbin_dat))
Zbin_dat <- Zbin_dat[,c('bin_age',lsoverch)] #prob should define a different df so you can check the data is the same.
check <- colnames(Zbin_dat)

check[!(check %in% lsoverch)]
lsoverch

## need dframe for curve generation below

Zbin_dat <- as.data.frame(Zbin_dat)
binname <- "bw500" 


#############smooth the records, only if you have filtered to 50% overlap i.e. its the subset you used for biclustering WITHOUT taking grid cell means. ####################

rm('i')
entvec <- colnames(Zbin_dat[,2:ncol(Zbin_dat)])

##get into long format

library(reshape2)

 zbinlong <- melt(Zbin_dat, id.vars='bin_age', variable.name='ID_ENTITY',
        value.name='bxcx', measure.vars = colnames(Zbin_dat[,2:ncol(Zbin_dat)]))
check <- zbinlong[zbinlong$ID_ENTITY %in% 5,]
plot(check$bin_age, check$bxcx, type = 'l')
lines(Zbin_dat$bin_age, Zbin_dat[,2], col = 'red')

entvec <- unique(zbinlong$ID_ENTITY)
smthz <- zbinlong[0,]
names(smthz) <- c("bin_age" ,  "ID_ENTITY", "smth40.bxcx" )
for(i in 1:length(entvec)){
#  i <- 6
sub <- zbinlong[zbinlong$ID_ENTITY %in% entvec[i],]
loessMod40 <- loess(bxcx ~ bin_age, data=sub, span=0.40) # 40% smoothing span

#join the smoothed output to the data
sub.tmp <- sub[!is.na(sub$bxcx),]
sub.tmp$smth40.bxcx <- predict(loessMod40)
sub.tmp <- sub.tmp[,c('bin_age','smth40.bxcx')]
sub <- merge(sub, sub.tmp, by = 'bin_age', all.x = T)
smthz <- rbind(smthz, sub)
}
warnings()

check <- df2[df2$ID_ENTITY %in% entvec[4],]
plot(check$year, check$mean, type = 'l', xlim = c(9000,0), lwd = 2)
lines(check$smth40.bxcx, x=check$year, col="red", lwd = 2)

check <- df[df$ID_ENTITY %in% entvec[2],]
plot(check$year, check$mean, type = 'l', xlim = c(9000,0), lwd = 2)
lines(check$smth40.bxcx, x=check$year, col="red", lwd = 2)

##

#bindat2.orig <- bindat2
df2$bxcx <- df2$smth40.bxcx #for ease of biclustering. dont have to rewrite stuff

check <- df2[df2$ID_ENTITY %in% entvec[2],]
plot(check$year, check$bxcx, type = 'l', xlim = c(9000,0), lwd = 2)
lines(check$smth40.bxcx, x=check$year, col="red", lwd = 2)

#anything else? 

#################################
#centre again before biclustering ******* this only makes sense to do if they are all the same length i.e., if you have one site that is mostly 9-4.5k and another that is 4.5-0k, and they would be from the same cluster and you centre on 0, that variation is not going to come out. It only makes sense if from the perspective of a fixed time point, all records start at 0.
#################################

##need df2 in wide format

recent <- as.data.frame(Zbin_dat$bin_age)
colnames(recent) <- 'bin_age'

entvec <- colnames(Zbin_dat[,!colnames(Zbin_dat) %in% 'bin_age'])
for(e in 1:length(entvec)){
	#e <- 1
  sub <- Zbin_dat[,colnames(Zbin_dat) %in% entvec[e]] #subset entity
  sub <- sub-mean(sub, na.rm = T) #subtract mean
  sub <- as.data.frame(sub)
  colnames(sub) <- entvec[e]
  recent <- cbind(recent, sub)
}

bindat <- recent

#should check the before and after here. do a plot. actually plot unbinned, binned, and recentered binned to see they make sense.
c1 <- as.vector(Zbin_dat[,7])
c2 <- as.vector(bindat[,7])
c3 <- as.vector(bindat$bin_age)

plot(c3,c1, type = 'l')
plot(c3,c2, type = 'l')

##now need to get it in long format again.


Zbin_dat <- as.data.frame(Zbin_dat[,1:20])
Zbin_dat.full <- Zbin_dat




#############


#####################################################



##########################################################################################################################################################

##################################
#COMPOSITE CURVE GENERATION
###################################

#purpose: select desired entities, then generate a composite curve for them.
#input data: binned z score data. 



#Load the locfit library and set parameters. The window width hw is, following convention, the half-window width (of the tricube weight function). The number of bootstrap resampling replications is given by nreps

#install.packages('locfit')
library(locfit)


# locfit (half) window-width parameter
hw <- 750 # bandwidth (smoothing parameter)

# number of bootstrap samples/replications
nreps <- 200

#Set the target ages for fitted values.

# target ages for fitted values
targbeg <- 100
targend <- 11850
targstep <- 250## this is based on the previous binning procedure

# array sizes
maxrecs <- 2000
maxreps <- 1000

# plot output 
plotout <- "pdf" # "screen" # "pdf"


##############

##filter to desired location:

#you have transformed data for the globe from dec 2020, so if you want a composite > 30N, you can just append 30-45 from there, before binning.
date <- '_06_05_2021'
site.filt <- mult[mult$ID_ENTITY %in% depo2,]
#must remove 605 because its 
head(meth2)
#site.filt <- site.filt[!site.filt$ID_ENTITY %in% c(127,126,434,435),]
site.filt <- site.filt[site.filt$depositional_context %in% 'lake sediment',] #'bog sediment'
s2 <- site.filt[site.filt$measurement_method %in% 'pollen slide',]
clustno <- 11
site.filt <- ent[ent$cluster %in% clustno,]
head(site.filt)
library(maps)
library(sf)
# world outline
world <- st_as_sf(map("world", plot = FALSE, fill = TRUE))
length(unique(site.filt$sitenum))
length(unique(site.filt$ID_ENTITY))
library(ggplot2)
names(site.filt)
range(site.filt$latitude.x);range(site.filt$longitude.x)
ggplot() +
  geom_sf(data = world, fill="grey80", color="black") +
  geom_point(data = site.filt, aes(x = longitude.x, y = latitude.x), shape = 19,col='red',  size = 4)+
	 scale_x_continuous(limits = c(15,40), breaks = seq(-180, 180, by=30)) + #-179, -50;-15,150
  theme(legend.text=element_text(size=18),legend.position = 'bottom',legend.title=element_blank() , panel.background = element_rect(fill = "aliceblue"), axis.text = element_text(size=14, face='bold', colour = 'black'), axis.title = element_text(size=14, face='bold', colour = 'black'), title = element_text(size = 10, colour = 'black'))+
  scale_y_continuous(limits = c(45,65), breaks = seq(30, 90, by=10)) +
  labs(x = "", y =  "") + #, title = "Site Clusters"
  coord_sf(expand = FALSE)+
	ggtitle('Sites with paired central/marginal entities (no other metadata differences)')
#sites <- ent[ent$ID_ENTITY %in% dat$ID_ENTITY,]

#site.filt <- sites[sites$longitude < (-50) & sites$latitude > 30,] 
site.filt <- as.character(site.filt$ID_ENTITY)
##to avoid overwriting original binned data - rename original - should actually go back and do this in the binning bit

#say you want to filter to a sub-region: read in data you saved that indicates which region each entity belongs to, and limit to the given region.
# datapath
# pwt <- read.csv('C:\\Users\\xr824576\\OneDrive - University of Reading\\PhD\\Data_analyses\\RPD\\Ambiclust_dat_bestRC_08_01_2021.csv')
# 
# site.filt <- unique(pwt[pwt$Site_cluster %in% 5, 'sitenum'])
# site.filt <- as.character(site.filt)
	#for clusters, need the values of the time clusters.
	
#Eurasia:
#America: 
Zbin_dat <- out[,colnames(out) %in%  c('bin_age',site.filt)]
z2 <- out[,colnames(out) %in%  c('bin_age',site.filt)]

#original one 	
# Zbin_dat <- Zbin_dat.full[,colnames(Zbin_dat.full) %in%  c('bin_age',site.filt)]
# Zbin_dat <- Zbin_dat[Zbin_dat$bin_age < 8000,]


##########################

## Z-score run

#define where to save the curve
curvecsvpath <- paste(datapath,"Regional_analyses\\",sep="")
curvecsvpath <- 'C:\\Users\\xr824576\\OneDrive - University of Reading\\Documents\\'

#define the file name of the statistics of the curve.
clustno <- 'miaro'
curvefile <- paste(curvecsvpath,"k", clustno,"_onlysizclass_var",basename, '0prcnt',"_",binname,"_",as.character(hw),"_",as.character(nreps),'_',as.character(length(Zbin_dat)-1),'n',date, ".csv", sep="") #adjust to particular data you use now. see README in curve folder in RPD folder for how to name if confused.
print(curvefile)


#Code block to implement writing .pdf to file: 

# .pdf plot of bootstrap iterations
if (plotout == "pdf") { ## when you define PDFfile over here, need to make sure you name it with the nsites and what baseperiod overlap level you insisted on.
  pdffile <- paste(curvecsvpath,"k", clustno,"_onlysizclass_var",basename, '0prcnt',"_",binname,"_",as.character(hw),"_",as.character(nreps),'_',as.character(length(Zbin_dat)-1),'n',date, ".pdf", sep="")
  print(pdffile)
}


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

length(ns) #important to note that this vector is one less than the number of columns of the binned data dataframe. careful when looping that you dont miss any entities. I use the loop counter as the length of the binned dataframe instead of the length of ns.
ls1 <- c()
for(i in 2:length(Zbin_dat)){ ##important that its length. not the actual character. NB that it skips the AGE COLUMN!!!!!!!!!!!!! which is why the ages might be better as last column of dframe.
 indata <- Zbin_dat[,c(i,1)] #1 is the age column. might be better last, and the loop starts at 1?
 indata <- na.omit(indata)  ##remove NaNs.
    nsamp <-  length(indata$bin_age) 
    if (nsamp > 0) {
      ii <- ii+1
      age[1:nsamp,ii] <- indata$bin_age # bins in >45n subset, age in original globe curve
      influx[1:nsamp,ii] <- indata[,1] # 
      nsamples[ii] <- nsamp
    } else{ls1 <- append(ls1, colnames(Zbin_dat[,i]))}
    print(i)
}

nsites <- ii
ls1 #which didnt have data. find out why.

# number of sites with data
nsites #why lose sites.. theres a site in Eurasia that gets lost.

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


head(cbind(targage,ndec_tot,xspan,ninwin_tot)) ##why is it cumulative?

#so is it calculating a point of the smoothed curve at 20 year intervals, using a 500 year window (is this 250 either side of the point) using all samples of all sites that fall within that window to calculate that point?


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
x <- as.vector(age) #ages for each entity
y <- as.vector(influx) #influx values of each bin for each entity
lfdata <- data.frame(x,y)
lfdata <- na.omit(lfdata) #exclude missing bins
x <- lfdata$x; y <- lfdata$y

# 2. locfit
# initial fit, unresampled (i.e. all) data
loc01 <- locfit(y ~ lp(x, deg=1, h=hw), maxk=800, family="qrgauss") #fit line through all pool of records.
summary(loc01)
plot(loc01)
#?locfit

# 3. get  fitted values
pred01 <- predict(loc01, newdata=targage.df, se.fit=TRUE)
loc01_fit <- data.frame(targage.df$x, pred01$fit)
fitname <- paste("locfit_",as.character(hw), sep="")
colnames(loc01_fit) <- c("age", fitname)
head(loc01_fit)
plot(loc01_fit$targage.df.x, loc01_fit$pred01.fit) #predicted vals are identical to 'initial fit' because same age values. it gives you 'se.fit = T' which I suppose gives you ses around the fit....

proc.time() - ptm


#Bootstrap-by-site confidence intervals
#The bootstrap confidence intervals are obtained by sampling with replacement sites (as opposed to individual samples), fitting a composite curve using locfit(), and saving these. 

#The 95-percent confidence intervals are then determined for each target age using the quantile() function.

# bootstrap samples
ptm <- proc.time()
# Step 1 -- Set up to plot individual replications - grey lines I assume
curvecsvpath
pdffile

##check y limits by scaling max and min of fitted vals

max(loc01_fit$locfit_750)
min(loc01_fit$locfit_750)
nsites

if (plotout == "pdf") {pdf(file=pdffile)} #png(filename=pdffile,res=600, width = 4000, height = 4000)
plot(x, y, xlab="Age (BP 1950)", ylab=fitname, xlim=c(5850,-100), xaxp  = c(5850, 0, 6), ylim=c(-6.5,1.8), type="n") 
#abline(v = c(4200 ,6600), lty = 2)
##make sure xlim is right.
#title(main = 'RPD 45N Eurasia filtered, bacon2020, Z-scores, n = 153, bp = 200-4200BP: min 80% overlap, bins = 250, hw = 750, nboot = 200', cex.main = 0.6)

#?par
#?set.seed
# Step 2 -- Do the bootstrap iterations, and plot each composite curve
set.seed(5) # do this to get the same sequence of random samples for each run

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
fitname <- paste("Zlocfit_",as.character(hw), sep="")
colnames(loc01_fit) <- c("age", fitname)
lines(loc01_fit[,1], loc01_fit[,2], lwd=2, col="red") ##keep all 45-55n America red; 55-65N america blue, 65-75 america orange; 45-55 eurasia purple, 55-65 eurasia green, 65-75 eurasia pink. #coral2; darkolivegreen3, turquoise3, lightgoldenrod4 royalblue2; deepskyblue3

# Step 4 -- Find and add bootstrap CIs
yfit95 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
yfit05 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
lines(targage.df$x, yfit95, lwd=1, col="red")
lines(targage.df$x, yfit05, lwd=1, col="red")##keep all 45-55n America red; 55-65N america blue, 65-75 america orange; 45-55 eurasia purple, 55-65 eurasia green, 65-75 eurasia pink.


if (plotout == "pdf") {dev.off()}

	curveout
curvefile
#Output
#The fitted curves are written out in the usual way.

curveout <- data.frame(cbind(targage.df$x, pred01$fit, yfit95, yfit05, ndec_tot, xspan, ninwin_tot))
colnames(curveout) <- c("age", "locfit", "cu95", "cl95", "nsites", "window", "ninwin")
curveout
write.csv(curveout, curvefile, row.names = F)
#outputfile <- paste(curvecsvpath, curvefile, sep="")
#write.table(curveout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
proc.time() - ptm


###################################################################################################################################
###############################   END ######################################################