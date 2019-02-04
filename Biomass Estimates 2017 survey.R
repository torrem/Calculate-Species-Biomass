###### bootstrap with boot package ##########

# load the package
library(boot)
library(data.table)
library(geosphere)

setwd("C:/Users/Mike/Documents/Scallop/Data/NGOM Biomass Estimates")

options(scipen=999)


#set scallop size cutoff ad dredge efficiency
size = 88.9



#create stratum area table
ngom.strat = matrix(c('MSI', 'MDI', 'PB', 'FL', 'NCAH', 'NCAM', 'NCAL', 'SJH', 'SJM','SJL','NSBH','NSBM','NSBL',
                      426, 1163, 148, 94, 142, 118, 145, 69, 33, 590, 196, 155, 526), ncol = 2)
colnames(ngom.strat) = c('STRATUM', 'AREA')


#-----------Establish Shell-Height Meat Weight Relationships--------------------

scallop  <- read.csv("C:/Users/Mike/Documents/Scallop/Data/NGOM Survey/LF.csv")



scallop = subset(scallop, !(TOW_STATION %in% c('MDI31','MDI32','MDI33','MDI34','MDI35',
                                               'MSI36','MSI37','MSI38','MSI39','MSI4','MSI40',  
                                               'MSI41','MSI42','MSI43','MSI44','MSI45','MSI46','MSI47',
                                               'NSBBC1','MSIBC1')))



shmw  <- read.csv('C:/Users/Mike/Documents/Scallop/Data/NGOM Survey/NGOM 2016 shmw.csv', header = TRUE)

shmw = subset(shmw, !(TOW_STATION %in% c('MDI31','MDI32','MDI33','MDI34','MDI35',
                                         'MSI36','MSI37','MSI38','MSI39','MSI4','MSI40',  
                                         'MSI41','MSI42','MSI43','MSI44','MSI45','MSI46','MSI47',
                                         'NSBBC1','MSIBC1')))
#--------------MSI
stlf = subset(scallop, STRATUM=="MSI")
stlf$TOW_STATION = factor(stlf$TOW_STATION) 
stmt = subset(shmw, STRATUM == "MSI")

lm2 = lm(log(SAMPLE_WEIGHT) ~ log(SAMPLE_LENGTH), data=stmt)

all <- split(stlf, stlf$TOW_STATION)

for (i in 1:length(all)) {
  x = all[[i]]
  
  lengths = as.numeric(x$SAMPLE_LENGTH)
  hist = hist(lengths,breaks = seq(0,200, by=5))
  breaks = hist$breaks
  breaks = breaks[2:length(breaks)]-2.5
  counts = hist$counts
  LF = data.table(cbind(breaks, counts))
  abundance = unique(x$REPORTED_QUANTITY)
  LF$proportion = LF$counts/sum(LF$counts)
  new = data.frame(SAMPLE_LENGTH = LF$breaks)
  
  n1L = predict(lm2, newdata = new)
  MTpred = exp(n1L)
  LF = cbind(LF,MTpred)
  LF$binTotal = LF$proportion*abundance
  LF$mass = LF$binTotal * LF$MTpred 
  LF = subset(LF, breaks > size)#scallop cutoff size
  biomass = sum(LF$mass)
  biomass = rep(biomass, nrow(x))
  
  x = cbind(x,biomass)
  all[[i]]=x
}

#recombine list into single data.frame
stlf <- do.call("rbind", all)
row.names(stlf)<-NULL

MSItows <- setDT(stlf)[,list(unique(EFFORT_START_DATE), unique(STRATUM), unique(DENSITY),unique(STRAT),
                             unique(TOW_START_LATDD), unique(TOW_START_LONDD), unique(HAULBACK_LATDD), unique(HAULBACK_LONDD),
                             unique(DEPTH), unique(BOTTOM_TYPE), unique(SPECIES_ITIS_CODE),
                             unique(REPORTED_QUANTITY), unique(biomass)),  by=list(TOW_STATION)]

setnames(MSItows,2,names(stlf)[1])
setnames(MSItows,1, names(stlf)[2])
setnames(MSItows,c(3:13),names(stlf)[3:13])
setnames(MSItows,14, names(stlf)[15])


#--------------MDI
stlf = subset(scallop, STRATUM=="MDI")
stlf$TOW_STATION = factor(stlf$TOW_STATION) 
stmt = subset(shmw, STRATUM == "MDI")

MDItows <- setDT(stlf)[,list(unique(EFFORT_START_DATE), unique(STRATUM), unique(DENSITY),unique(STRAT),
                             unique(TOW_START_LATDD), unique(TOW_START_LONDD), unique(HAULBACK_LATDD), unique(HAULBACK_LONDD),
                             unique(DEPTH), unique(BOTTOM_TYPE), unique(SPECIES_ITIS_CODE),
                             unique(REPORTED_QUANTITY)),  by=list(TOW_STATION)]

setnames(MDItows,2,names(stlf)[1])
setnames(MDItows,1, names(stlf)[2])
setnames(MDItows,c(3:13),names(stlf)[3:13])

MDItows$biomass = c(rep(0,16), 27, rep(0,26))

#--------------FL
stlf = subset(scallop, STRATUM=="FL")
stlf$TOW_STATION = factor(stlf$TOW_STATION) 
stmt = subset(shmw, STRATUM == "FL")

lm2 = lm(log(SAMPLE_WEIGHT) ~ log(SAMPLE_LENGTH), data=stmt)

all <- split(stlf, stlf$TOW_STATION)

for (i in 1:length(all)) {
  x = all[[i]]
  
  lengths = as.numeric(x$SAMPLE_LENGTH)
  hist = hist(lengths,breaks = seq(0,200, by=5))
  breaks = hist$breaks
  breaks = breaks[2:length(breaks)]-2.5
  counts = hist$counts
  LF = data.table(cbind(breaks, counts))
  abundance = unique(x$REPORTED_QUANTITY)
  LF$proportion = LF$counts/sum(LF$counts)
  new = data.frame(SAMPLE_LENGTH = LF$breaks)
  
  n1L = predict(lm2, newdata = new)
  MTpred = exp(n1L)
  LF = cbind(LF,MTpred)
  LF$binTotal = LF$proportion*abundance
  LF$mass = LF$binTotal * LF$MTpred 
  LF = subset(LF, breaks > size)#scallop cutoff size
  biomass = sum(LF$mass)
  biomass = rep(biomass, nrow(x))
  x = cbind(x,biomass)
  all[[i]]=x
}

#recombine list into single data.frame
stlf <- do.call("rbind", all)
row.names(stlf)<-NULL

FLtows <- setDT(stlf)[,list(unique(EFFORT_START_DATE), unique(STRATUM), unique(DENSITY),unique(STRAT),
                            unique(TOW_START_LATDD), unique(TOW_START_LONDD), unique(HAULBACK_LATDD), unique(HAULBACK_LONDD),
                            unique(DEPTH), unique(BOTTOM_TYPE), unique(SPECIES_ITIS_CODE),
                            unique(REPORTED_QUANTITY), unique(biomass)),  by=list(TOW_STATION)]

setnames(FLtows,2,names(stlf)[1])
setnames(FLtows,1, names(stlf)[2])
setnames(FLtows,c(3:13),names(stlf)[3:13])
setnames(FLtows,14, names(stlf)[15])

#--------------PB
stlf = subset(scallop, STRATUM=="PB")
stlf$TOW_STATION = factor(stlf$TOW_STATION) 
stmt = subset(shmw, STRATUM == "PB")

lm2 = lm(log(SAMPLE_WEIGHT) ~ log(SAMPLE_LENGTH), data=stmt)

all <- split(stlf, stlf$TOW_STATION)

for (i in 1:length(all)) {
  x = all[[i]]
  
  lengths = as.numeric(x$SAMPLE_LENGTH)
  hist = hist(lengths,breaks = seq(0,200, by=5))
  breaks = hist$breaks
  breaks = breaks[2:length(breaks)]-2.5
  counts = hist$counts
  LF = data.table(cbind(breaks, counts))
  abundance = unique(x$REPORTED_QUANTITY)
  LF$proportion = LF$counts/sum(LF$counts)
  new = data.frame(SAMPLE_LENGTH = LF$breaks)
  
  n1L = predict(lm2, newdata = new)
  MTpred = exp(n1L)
  LF = cbind(LF,MTpred)
  LF$binTotal = LF$proportion*abundance
  LF$mass = LF$binTotal * LF$MTpred 
  LF = subset(LF, breaks > size)#scallop cutoff size
  biomass = sum(LF$mass)
  biomass = rep(biomass, nrow(x))
  x = cbind(x,biomass)
  all[[i]]=x
}

#recombine list into single data.frame
stlf <- do.call("rbind", all)
row.names(stlf)<-NULL

PBtows <- setDT(stlf)[,list(unique(EFFORT_START_DATE), unique(STRATUM), unique(DENSITY),unique(STRAT),
                            unique(TOW_START_LATDD), unique(TOW_START_LONDD), unique(HAULBACK_LATDD), unique(HAULBACK_LONDD),
                            unique(DEPTH), unique(BOTTOM_TYPE), unique(SPECIES_ITIS_CODE),
                            unique(REPORTED_QUANTITY), unique(biomass)),  by=list(TOW_STATION)]

setnames(PBtows,2,names(stlf)[1])
setnames(PBtows,1, names(stlf)[2])
setnames(PBtows,c(3:13),names(stlf)[3:13])
setnames(PBtows,14, names(stlf)[15])

#------------NCA
stlf = subset(scallop, STRATUM=="NCA")
stlf$TOW_STATION = factor(stlf$TOW_STATION) 
stmt = subset(shmw, STRATUM == "NCA")

lm2 = lm(log(SAMPLE_WEIGHT) ~ log(SAMPLE_LENGTH), data=stmt)

all <- split(stlf, stlf$TOW_STATION)

for (i in 1:length(all)) {
  x = all[[i]]
  
  lengths = as.numeric(x$SAMPLE_LENGTH)
  hist = hist(lengths,breaks = seq(0,200, by=5))
  breaks = hist$breaks
  breaks = breaks[2:length(breaks)]-2.5
  counts = hist$counts
  LF = data.table(cbind(breaks, counts))
  abundance = unique(x$REPORTED_QUANTITY)
  LF$proportion = LF$counts/sum(LF$counts)
  new = data.frame(SAMPLE_LENGTH = LF$breaks)
  
  n1L = predict(lm2, newdata = new)
  MTpred = exp(n1L)
  LF = cbind(LF,MTpred)
  LF$binTotal = LF$proportion*abundance
  LF$mass = LF$binTotal * LF$MTpred 
  LF = subset(LF, breaks > size)#scallop cutoff size
  biomass = sum(LF$mass)
  biomass = rep(biomass, nrow(x))
  x = cbind(x,biomass)
  all[[i]]=x
}

#recombine list into single data.frame
stlf <- do.call("rbind", all)
row.names(stlf)<-NULL

NCAtows <- setDT(stlf)[,list(unique(EFFORT_START_DATE), unique(STRATUM), unique(DENSITY),unique(STRAT),
                             unique(TOW_START_LATDD), unique(TOW_START_LONDD), unique(HAULBACK_LATDD), unique(HAULBACK_LONDD),
                             unique(DEPTH), unique(BOTTOM_TYPE), unique(SPECIES_ITIS_CODE),
                             unique(REPORTED_QUANTITY), unique(biomass)),  by=list(TOW_STATION)]

setnames(NCAtows,2,names(stlf)[1])
setnames(NCAtows,1, names(stlf)[2])
setnames(NCAtows,c(3:13),names(stlf)[3:13])
setnames(NCAtows,14, names(stlf)[15])


#------SJ
stlf = subset(scallop, STRATUM=="SJ")
stlf$TOW_STATION = factor(stlf$TOW_STATION) 
stmt = subset(shmw, STRATUM == "SJ")

lm2 = lm(log(SAMPLE_WEIGHT) ~ log(SAMPLE_LENGTH), data=stmt)

all <- split(stlf, stlf$TOW_STATION)

for (i in 1:length(all)) {
  x = all[[i]]
  
  lengths = as.numeric(x$SAMPLE_LENGTH)
  hist = hist(lengths,breaks = seq(0,200, by=5))
  breaks = hist$breaks
  breaks = breaks[2:length(breaks)]-2.5
  counts = hist$counts
  LF = data.table(cbind(breaks, counts))
  abundance = unique(x$REPORTED_QUANTITY)
  LF$proportion = LF$counts/sum(LF$counts)
  new = data.frame(SAMPLE_LENGTH = LF$breaks)
  
  n1L = predict(lm2, newdata = new)
  MTpred = exp(n1L)
  LF = cbind(LF,MTpred)
  LF$binTotal = LF$proportion*abundance
  LF$mass = LF$binTotal * LF$MTpred 
  LF = subset(LF, breaks > size)#scallop cutoff size
  biomass = sum(LF$mass)
  biomass = rep(biomass, nrow(x))
  x = cbind(x,biomass)
  all[[i]]=x
}

#recombine list into single data.frame
stlf <- do.call("rbind", all)
row.names(stlf)<-NULL

SJtows <- setDT(stlf)[,list(unique(EFFORT_START_DATE), unique(STRATUM), unique(DENSITY),unique(STRAT),
                            unique(TOW_START_LATDD), unique(TOW_START_LONDD), unique(HAULBACK_LATDD), unique(HAULBACK_LONDD),
                            unique(DEPTH), unique(BOTTOM_TYPE), unique(SPECIES_ITIS_CODE),
                            unique(REPORTED_QUANTITY), unique(biomass)),  by=list(TOW_STATION)]

setnames(SJtows,2,names(stlf)[1])
setnames(SJtows,1, names(stlf)[2])
setnames(SJtows,c(3:13),names(stlf)[3:13])
setnames(SJtows,14, names(stlf)[15])

#----------NSB
stlf = subset(scallop, STRATUM=="NSB")
stlf$TOW_STATION = factor(stlf$TOW_STATION) 
stmt = subset(shmw, STRATUM == "NSB")

lm2 = lm(log(SAMPLE_WEIGHT) ~ log(SAMPLE_LENGTH), data=stmt)

all <- split(stlf, stlf$TOW_STATION)

for (i in 1:length(all)) {
  x = all[[i]]
  
  lengths = as.numeric(x$SAMPLE_LENGTH)
  hist = hist(lengths,breaks = seq(0,200, by=5))
  breaks = hist$breaks
  breaks = breaks[2:length(breaks)]-2.5
  counts = hist$counts
  LF = data.table(cbind(breaks, counts))
  abundance = unique(x$REPORTED_QUANTITY)
  LF$proportion = LF$counts/sum(LF$counts)
  new = data.frame(SAMPLE_LENGTH = LF$breaks)
  
  n1L = predict(lm2, newdata = new)
  MTpred = exp(n1L)
  LF = cbind(LF,MTpred)
  LF$binTotal = LF$proportion*abundance
  LF$mass = LF$binTotal * LF$MTpred 
  LF = subset(LF, breaks > size)#scallop cutoff size
  biomass = sum(LF$mass)
  biomass = rep(biomass, nrow(x))
  x = cbind(x,biomass)
  all[[i]]=x
}

#recombine list into single data.frame
stlf <- do.call("rbind", all)
row.names(stlf)<-NULL

NSBtows <- setDT(stlf)[,list(unique(EFFORT_START_DATE), unique(STRATUM), unique(DENSITY),unique(STRAT),
                             unique(TOW_START_LATDD), unique(TOW_START_LONDD), unique(HAULBACK_LATDD), unique(HAULBACK_LONDD),
                             unique(DEPTH), unique(BOTTOM_TYPE), unique(SPECIES_ITIS_CODE),
                             unique(REPORTED_QUANTITY), unique(biomass)),  by=list(TOW_STATION)]

setnames(NSBtows,2,names(stlf)[1])
setnames(NSBtows,1, names(stlf)[2])
setnames(NSBtows,c(3:13),names(stlf)[3:13])
setnames(NSBtows,14, names(stlf)[15])



#----------Combine all stratum--------------------

dat16 = rbind(MSItows, MDItows, FLtows, PBtows, NCAtows,SJtows,NSBtows)
dat16$biomass[is.na(dat16$biomass)] <- 0


#-----------Calculate tow area and convert biomass per tow into biomass per meter----

tow_length = matrix(NA, nrow = nrow(dat16), ncol = 1)

for (i in 1:nrow(dat16)){
  n = distm(c(dat16$TOW_START_LONDD[i], dat16$TOW_START_LATDD[i]), 
            c(dat16$HAULBACK_LONDD[i], dat16$HAULBACK_LATDD[i]), fun=distHaversine)
  tow_length[[i]] = n
}

dat16 = cbind(dat16,tow_length)
colnames(dat16)[15] = "tow_length"



##Alter tow values
dat16[22,15] = 500.0495
dat16[65,15] = 470.7662
dat16[1,15] = 457.015
dat16[34,15] = 240.2215
dat16[188,15] = 523.3174
dat16[189,15] = 497.1139
dat16[159,15] = 487.9227
dat16[203,15] = 615.9446
dat16[213,15] = 464.9962
dat16[230,15] = 526.4861
dat16[185,15] = 437.0247

dat16$tow_area = dat16$tow_length * 2.13
dat16$BM_M2 = dat16$biomass/dat16$tow_area



#----------------------- Separate out NSB & SJ tows inside 2017 survey area--------------

stat17sep  <- read.csv("C:/Users/Mike/Documents/Scallop/Data/NGOM Biomass Estimates/2016_2017 IN_OUT Stations.csv")

dat16$INOUT17 = NA
dat16$INOUT17 = as.factor(dat16$INOUT17)
f = match(stat17sep$TOW_STATION,dat16$TOW_STATION,nomatch=-999)

for (i in 1:length(f))
{
y=f[i]
h = as.factor(stat17sep[i,1])
dat16[y,18] = h
}




#------------------Bootstrap--------------


SJ_IN = subset(dat16, STRATUM=="SJ" & INOUT17 == "SJ_IN")
SJ_OUT = subset(dat16, STRATUM=="SJ" & INOUT17 == "SJ_OUT")

NSB_IN = subset(dat16, STRATUM=="NSB" & INOUT17 == "NSB_IN")
NSB_OUT = subset(dat16, STRATUM=="NSB" & INOUT17 == "NSB_OUT")

# make the number B
B=10000
# make blank lists 10,000 long

SJ_IN.v=numeric(B)
SJ_OUT.v=numeric(B)

NSB_IN.v=numeric(B)
NSB_OUT.v=numeric(B)


for (i in 1:B){
 
  SJ_IN.v[i]=mean(sample(SJ_IN$BM_M2,size=nrow(SJ_IN),replace=TRUE))
  SJ_OUT.v[i]=mean(sample(SJ_OUT$BM_M2,size=nrow(SJ_OUT),replace=TRUE))

  NSB_IN.v[i]=mean(sample(NSB_IN$BM_M2,size=nrow(NSB_IN),replace=TRUE))
  NSB_OUT.v[i]=mean(sample(NSB_OUT$BM_M2,size=nrow(NSB_OUT),replace=TRUE))
 
}



#--------------------Tables-------------------------------------


areas = read.csv('C:/Users/Mike/Documents/Scallop/Data/NGOM Biomass Estimates/areas2017inout.csv')
areas = areas[c(3,4,1,2),] #reorder areas to match order in x
areas[c(3:4,1:2),2]=c(189,688,81,611)

#summary list
x = list(SJ_IN.v, SJ_OUT.v, NSB_IN.v, NSB_OUT.v)

summary_stats = matrix(NA, nrow = length(x), ncol = 8) 
colnames(summary_stats) = c('q0.025', 'q0.10', 'q0.15', 'q0.20', 'q0.25', 'Mean', 'q0.975')
rownames(summary_stats) = c('SJ_IN', 'SJ_OUT', 'NSB_IN', 'NSB_OUT')


dredge = 0.4

for (i in 1:length(x)){
  y=x[[i]]
  a = as.numeric(areas[i,2])
  y = y*1000*1000*a/dredge/1000000
  summary_stats[i,1] = quantile(y, 0.025)
  summary_stats[i,2] = quantile(y, 0.10)
  summary_stats[i,3] = quantile(y, 0.15)
  summary_stats[i,4] = quantile(y, 0.20)
  summary_stats[i,5] = quantile(y, 0.25)
  summary_stats[i,6] = mean(y)
  summary_stats[i,7] = quantile(y, 0.975)
   summary_stats[i,8] = sd(y)/sqrt(length(y))
}

# s = barplot(height =  summary_stats[,4], ylim = c(0,2000), xlab = '', ylab = '')
# segments(x0 = s, y0 = summary_stats[,7], x1 = s, y1 = summary_stats[,8])

#computation of the standard error of the mean
sem<-sd(x)/sqrt(length(x))
#95% confidence intervals of the mean
c(mean(x)-2*sem,mean(x)+2*sem)

