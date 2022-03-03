

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

RPDent <-  dbGetQuery(mydb, "SELECT e.* FROM entity e 
                                 WHERE ID_ENTITY in(295, 527, 528, 645, 857, 939, 941, 967, 1007, 2154, 2222, 2224)
                                     and TYPE != 'other' and 
                                     TYPE !='raw count';") #longer site dataset is because of those from WP that have no entities.
str(RPDent)
length(unique(RPDent$ID_SITE))
length(unique(RPDent$entity_name))


RPDsite <- dbGetQuery(mydb, "SELECT s.* FROM site s
                      left join entity e on s.ID_SITE = e.ID_SITE 
                        WHERE ID_ENTITY in(295, 527, 528, 645, 857, 939, 941, 967, 1007, 2154, 2222, 2224)
                                     and TYPE != 'other' and 
                                     TYPE !='raw count';")
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

write.csv(RPDdat, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/SSHAFat.csv", row.names = F)
write.csv(RPDsite, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/SHAFsite.csv", row.names = F)
write.csv(RPDent, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/SHAFent.csv", row.names = F)
write.csv(RPDchron, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/SHAFchron.csv", row.names = F)


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


ent <- read.csv('/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/SHAFent.csv', na.strings = 'NULL')

dat <- read.csv('/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/SSHAFat.csv' , na.strings = 'NULL')

chron <- read.csv('/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/SHAFchron.csv' , na.strings = 'NULL')

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

str(data_query)
data_query$original_est_age.x <- as.numeric(data_query$original_est_age.x)
data_query$avg_depth <- as.numeric(data_query$avg_depth)

#loop 1 ----
for (j in maxsites[]) { ##trying looping using entity names: changed from j in 1:maxsites
   nsamp <- 0
  sitedata <- data_query[data_query$ID_ENTITY ==  j, ]
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
    influx <- as.numeric(quant)
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
    ls1 <- append(ls1, entnamevec)
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
    
    outdata <- data.frame(entnamevec , entidvec, sampid, depth, age, sed_rate, quant, conc, influx, xst_level, conc_source, influx_source, xst_level)#'* note PL has changed the output here from new_type to repeat xst_level to prevent overwriting original RPD types * 
    names(outdata) <- c('entity_name',"ID_ENTITY", 'ID_SAMPLE',"depth", "est_age", "sed_rate", "quant", "conc",
                        "influx", "xst_level", "conc_source", "influx_source", 'TYPE' ) #might change est age to original est age... for clarity. Also write out sample IDs, so you can bind new bacon chrons.
    
    
    #delete any rows in this data that have zero and negative sed rate, and save them so you can see how many are removed.
    neg.infl <- append(neg.infl, outdata[outdata$sed_rate <= 0,])
    outdata <- outdata[!(outdata$sed_rate <= 0),]
    newdf <- rbind(newdf, outdata)
    
      } 
  print(j)
}

unent<-unique(newdf$ID_ENTITY)
length(unent)


##########################################################################################################################################################

##################################             
#INFLUX CALCULATION                  
##################################
#purpose: calculate influx - write out a file for each entity with influx values, then read it in to bind into one dataframe. if there are new chronologies, use these.
#input data: RPD data cleaned and filtered to one entity per site.

#before doing this, will ensure that the data is all in concentration and has the new bacon age model, which will be used to generate influx. 

rm(list= ls()[!(ls() %in% c('dat','site','newdf'))])

# query label and path to query and query name
dir.create("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/SHAF")
csvpath <- "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/SHAF/"


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
data_query[data_query$TYPE %in% c("concentration", 'Concentration'), 'TYPE'] <- 'CONC'
data_query[data_query$TYPE %in% "pollen concentration",'TYPE'] <- 'C0P0'
data_query[data_query$TYPE %in% "INFLUX",'TYPE'] <- 'INFL'
data_query[data_query$TYPE %in% "influx",'TYPE'] <- 'INFL'
data_query[data_query$TYPE %in% c("raw count",'count', 'Count'),'TYPE'] <- 'COUN'
unique(data_query$TYPE)
head(data_query)

coun <- data_query[data_query$TYPE %in% c('COUN', 'soil','other') | is.na(data_query$TYPE),] ##for <45N:

#they are all counts, so should just fix it in the loop that gets conc (marked where). for now can leave it.

##in general in this script: at this point you want to exclude from this whole loop any entities that cant possibly be transformed to influx (e.g. type = COUNT.) per unit weight, conc, influx C0P0 you want to keep. so purge any that dont fit this HERE. - probably redundant given above getting conc and filtering for that.
coun <- data_query[data_query$TYPE %in% c('raw count','other') | is.na(data_query$TYPE),] ##soil is for global, the NA is temporary, should be resolved when resolving issue in conc transforming loop.
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
print(maxsites[])
match(1786, maxsites)
data_query$quantity <- as.numeric(data_query$quantity) #new piece of code to convert quantities into numeric rather than character values
str(data_query)
#INFLUX LOOP----

for (j in maxsites[]) { ##trying looping using entity names: changed from j in 1:maxsites
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
    
    outdata <- data.frame( entidvec, sampidvec, depth, age, sed_rate, quant, conc, influx, xst_level, conc_source, influx_source) 
    names(outdata) <- c("ID_ENTITY", "ID_SAMPLE", "depth", "est_age", "sed_rate", "quant", "conc",
                        "influx", "xst_level", "conc_source", "influx_source" )
    
   #delete any rows in this data that have zero and negative sed rate, and save them so you can see how many are removed.
    neg.infl <- append(neg.infl, outdata[outdata$sed_rate <= 0,])
    outdata <- outdata[!(outdata$sed_rate <= 0),]
    csvfile <- paste("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/SHAF/",entid,"_RPD_", ".csv", sep="")
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
length(ls3)
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

sitecsvpath <- paste("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/SHAF/",sep="") #where all the csvs that have the influx data are stored.
sitecsvpath
transcsvpath <- paste("/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/ZTrans_SHAF/",sep="") # 
transcsvpath

# if output folder does not exist, create it
queryname <- 'SHAF'
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
#Transform loop ----
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


write.csv(ent_vec, "/Users/paullincoln/Dropbox/2021/September 2021/Composite code v2/Influx conversion/SHAF/final_entities.csv")

##########################################  END OF MAIN LOOPING BIT ###################################################################
