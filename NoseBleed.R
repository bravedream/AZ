#######################################################################
#Technical Scenario modeling for
#Evaluation of Treatment effect of Nasal Mucus Viscosity on Nose Bleeds
#Date: April 2016
#Author: Xiang Ji
#######################################################################


##read in csv
set.seed(1985)
filepath<-"C:/Users/xji/Desktop/Nasal"  
eff<- read.csv(file.path(filepath, "efficacy.csv"),header=T, na.strings=c("NULL",""),stringsAsFactors=F)
sub<- read.csv(file.path(filepath, "subject.csv"),header=T, na.strings=c("NULL",""),stringsAsFactors=F)
random<- read.csv(file.path(filepath, "randomization.csv"),header=T, na.strings=c("NULL",""),stringsAsFactors=F)
str(eff)
str(sub)
str(random)

#merge files
dat<-merge(random,merge(eff,sub,by="subject",all=T),by="subject",all=T)
str(dat)
#derive variables to be used
dat$mucus.viscosity<-as.numeric(dat$mucus.viscosity)
dat$country<-as.factor(dat$country)
dat<-na.omit(dat) #delete one  observation with missing mucus viscosity
names(dat)[names(dat)=="previous.year"]="hosp.bleed"
summary(dat)
apply(dat, 2, function(x) sum(is.na(x)))
piefun<-function(x) {
	lbls<-paste(names(table(x)),round(prop.table(table(x))*100), "%")
	pie(table(x), label=lbls, main=deparse(substitute(x)))
}
piefun(dat$arm)
ggplot(dat,aes(x=nosebleeds))+geom_density(alpha=0.5)
ggplot(dat,aes(x=duration))+geom_density(alpha=0.5)
ggplot(dat,aes(x=nosebleeds))+geom_density(alpha=0.5)
piefun(dat$country)
piefun(dat$eye.colour)
piefun(dat$tissue.use)
ggplot(dat,aes(x=hosp.bleed))+geom_density(alpha=0.5)
ggplot(dat,aes(x=mucus.viscosity))+geom_density(alpha=0.5)
#normality test for continuous
apply(dat[,c("nosebleeds","duration","hosp.bleed","mucus.viscosity")],2,function(x) shapiro.test(x)) #all violated normality
#test independence of baseline continuous variables between arms (wilcoxon for non-normal variables)
apply(dat[,c("duration","hosp.bleed","mucus.viscosity")],2,function(x) wilcox.test(x ~arm,data=dat))
#test difference between arms for baseline categorical variables
apply(dat[,c("country","eye.colour","tissue.use")],2,function(x) summary(table(dat$arm,x)))

##Exploratory Plot
#install.packages("ggplot2")
library(ggplot2)
d<-ggplot(dat,aes(x=nosebleeds))
d+geom_density(aes(group=arm, colour=arm, fill=arm), alpha=0.3)
#install.packages("lattice")
library(lattice)
#bleed to viscosity by arm
xyplot(nosebleeds ~ mucus.viscosity, group=arm, data=dat, 
       auto.key=list(space="right"), 
       jitter.x=TRUE, jitter.y=TRUE,type = c("p","r"))
summary(lm(nosebleeds~arm+duration+mucus.viscosity+arm:mucus.viscosity,data=dat))
t.test(nosebleeds~arm,data=dat[dat$mucus.viscosity<=summary(dat$mucus.viscosity)[[2]],]) #viscosity<=25 quartile
t.test(nosebleeds~arm,data=dat[dat$mucus.viscosity>summary(dat$mucus.viscosity)[[2]] & dat$mucus.viscosity<=summary(dat$mucus.viscosity)[[5]],]) #viscosity between 25 and 75 quartile
t.test(nosebleeds~arm,data=dat[dat$mucus.viscosity>summary(dat$mucus.viscosity)[[5]],]) #viscosity above 75 quartile

#bleed to tissue use by arm
d<-ggplot(dat,aes(x=nosebleeds))
d+geom_density(aes(group=arm, colour=arm, fill=arm), alpha=0.3)+facet_wrap(~tissue.use,nrow=2)
summary(lm(nosebleeds~arm+tissue.use+arm:tissue.use,data=dat))
t.test(nosebleeds~arm,data=dat[dat$tissue.use=="MEDIUM",])
t.test(nosebleeds~arm,data=dat[dat$tissue.use=="HIGH",])

#hospitalization due to nose bleed by country
d<-ggplot(dat,aes(x=hosp.bleed))
d+geom_density(aes(group=country, colour=country, fill=country), alpha=0.3)+
facet_wrap(~country,ncol= nlevels(dat$country))

#bleed by eye color by arm
d<-ggplot(dat,aes(x=nosebleeds))
d+geom_density(aes(group=arm, colour=arm, fill=arm), alpha=0.3)+facet_wrap(~eye.colour,nrow=2)

