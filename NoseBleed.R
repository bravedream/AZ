#######################################################################
#Technical Scenario modeling for the
#Evaluation of Treatment effect of Nasal Mucus Viscosity on Nose Bleeds
#Date: April 2016
#Author: Xiang Ji
#######################################################################

#load/install libraries
library(gridExtra)
library(car)
library(lattice)
library(ggplot2)
library(MASS)
library(caret)
library(pscl)


##read in csv files
set.seed(1985)
filepath<-"C:/Users/xji/Desktop/Nasal"    ##please change the file folder to where files are saved
eff<- read.csv(file.path(filepath, "efficacy.csv"),header=T, na.strings=c("NA",""),stringsAsFactors=F)
sub<- read.csv(file.path(filepath, "subject.csv"),header=T, na.strings=c("NA",""),stringsAsFactors=F)
random<- read.csv(file.path(filepath, "randomization.csv"),header=T, na.strings=c("NA",""),stringsAsFactors=F)
str(eff)
str(sub)
str(random)

#merge files
dat<-merge(random,merge(eff,sub,by="subject",all=T),by="subject",all=T)
str(dat)  #examine data structure
#derive variables to be used
dat$mucus.viscosity<-as.numeric(dat$mucus.viscosity)  #convert to numeric 
dat$country<-as.factor(dat$country)  #convert to factor
dat$arm<-as.factor(dat$arm)  #convert to factor
#check possible values of each variable to see if there is unexpected values
apply(dat, 2, function(x) table(x)) 
#check number of missing for each variable
apply(dat, 2, function(x) sum(is.na(x))) ##found 60 NAs in eye color, which is fine because it can be treated as a special category. 1 missing in continuous variable viscosity
#impute/re-categorize missing
dat[is.na(dat$eye.colour),"eye.colour"]<-"UNKNOWN"   #keep the NAs and replace with unknown to avoid losing 60 sample
dat[is.na(dat$mucus.viscosity),"mucus.viscosity"]<-median(dat[,"mucus.viscosity"], na.rm=T)  #simple impute the median for the missing
#check summary statistics
summary(dat)  
apply(dat, 2, function(x) sum(is.na(x)))  #recheck, no missing
names(dat)[names(dat)=="previous.year"]="hosp.bleed"  #change to more meaningful name for hospitalization due to bleed

piefun<-function(x,title) {
	lbls<-paste(names(table(x)),round(prop.table(table(x))*100), "%")
	pie(table(x), label=lbls, main=title)
}
par(mfrow=c(2,2))
piefun(dat$arm,"Treatment Arms")
piefun(dat$country,"Country")
piefun(dat$eye.colour, "Eye Colour")
piefun(dat$tissue.use, "Baseline Tissue Usage")
boxplot(dat$nosebleeds,dat$hosp.bleed,dat$mucus.viscosity)
p1<-ggplot(dat,aes(x=nosebleeds))+geom_density(alpha=0.5)+ggtitle("Number of Nose Bleeds")
p2<-ggplot(dat,aes(x=duration))+geom_density(alpha=0.5)+ggtitle("Trial Duration of Follow Up")
p3<-ggplot(dat,aes(x=hosp.bleed))+geom_density(alpha=0.5)+ggtitle("Number of Hosp from Bleeds Previous Year")
p4<-ggplot(dat,aes(x=mucus.viscosity))+geom_density(alpha=0.5)+ggtitle("Mucus Viscosity at Baseline")
grid.arrange(p1, p2, p3, p4, ncol=2, nrow =2)

#normality test for continuous
apply(dat[,c("nosebleeds","duration","hosp.bleed","mucus.viscosity")],2,function(x) shapiro.test(x)) #all violated normality
#test independence of baseline continuous variables between arms (wilcoxon for non-normal variables)-- no difference found between arms
aggregate(. ~ arm,data = dat[,c("arm","duration","hosp.bleed","mucus.viscosity")], FUN=function(x) c(mean =round(mean(x),2), SD=round(sd(x),2)))  #extract mean and SD
apply(dat[,c("duration","hosp.bleed","mucus.viscosity")],2,function(x) wilcox.test(x ~arm,data=dat)) #wilcoxon non-parametric test 
#test difference between arms for baseline categorical variables  -- no difference found between arms
apply(dat[,c("country","eye.colour","tissue.use")],2,function(x) round(prop.table(table(dat$arm,x),1)*100,2)) #frequency/p -value
apply(dat[,c("country","eye.colour","tissue.use")],2,function(x) summary(table(dat$arm,x))) #chisq p value

##Exploratory Plot
#correlation matrix
scatterplotMatrix(~nosebleeds+hosp.bleed+duration+mucus.viscosity|arm, by.group=T,lower.panel=NULL,legend.pos="center", data=dat,main="")
#by treatment difference
summary(lm(nosebleeds~arm,data=dat))
round(wilcox.test(nosebleeds ~arm,data=dat)$p.value,4) #viscosity<=25 quartile
d<-ggplot(dat,aes(x=nosebleeds))
d+geom_density(aes(group=arm, colour=arm, fill=arm),size=2, alpha=0.3)+ggtitle("Number of Nose Bleeds by Arm")+ annotate("text", x = 4, y = 1, label=paste0("p value: ",pval))

#bleed to viscosity by arm
wilcox.test(nosebleeds ~arm,data=dat[dat$mucus.viscosity<=summary(dat$mucus.viscosity)[[2]],]) #viscosity<=25 quartile
wilcox.test(nosebleeds ~arm,data=dat[dat$mucus.viscosity>summary(dat$mucus.viscosity)[[2]] & dat$mucus.viscosity<=summary(dat$mucus.viscosity)[[5]],])  #viscosity between 25 and 75 quartile
wilcox.test(nosebleeds ~arm,data=dat[dat$mucus.viscosity>summary(dat$mucus.viscosity)[[5]],])  #viscosity above 75 quartile
xyplot(nosebleeds ~ mucus.viscosity, group=arm, data=dat, 
       auto.key=list(corner = c(0.95,0.9 ), pch=10,cex=1.0), 
       jitter.x=TRUE, jitter.y=TRUE,type = c("p","r"), main="Bleeds and Mucus Viscosity" )
summary(lm(nosebleeds~arm+mucus.viscosity+arm:mucus.viscosity,data=dat))

#bleed to tissue use by arm
pvalM<-round(wilcox.test(nosebleeds~arm,data=dat[dat$tissue.use=="MEDIUM",])$p.value,4)
pvalH<-round(wilcox.test(nosebleeds~arm,data=dat[dat$tissue.use=="HIGH",])$p.value,4)
ggplot(dat,aes(x=nosebleeds))+geom_density(aes(group=arm, colour=arm, fill=arm), size=2,alpha=0.3)+facet_wrap(~tissue.use,ncol=2)+ggtitle("Nose Bleeds by Tissue Usage")+ ylab("Count")

summary(lm(nosebleeds~arm+tissue.use+arm:tissue.use,data=dat))

#bleed by eye color by arm
wilcox.test(nosebleeds~arm,data=dat[dat$eye.colour=="BLUE",])
wilcox.test(nosebleeds~arm,data=dat[dat$eye.colour=="BLACK",])
wilcox.test(nosebleeds~arm,data=dat[dat$eye.colour=="BROWN",])
wilcox.test(nosebleeds~arm,data=dat[dat$eye.colour=="UNKNOWN",])
d<-ggplot(dat,aes(x=nosebleeds))
d+geom_density(aes(group=arm, colour=arm, fill=arm),size=2, alpha=0.3)+facet_wrap(~eye.colour,nrow=2)
#bleed by country by arm

wilcox.test(nosebleeds~arm,data=dat[dat$country=="A",])
wilcox.test(nosebleeds~arm,data=dat[dat$eye.colour=="BLACK",])
wilcox.test(nosebleeds~arm,data=dat[dat$eye.colour=="BROWN",])
wilcox.test(nosebleeds~arm,data=dat[dat$eye.colour=="UNKNOWN",])

d<-ggplot(dat,aes(x=nosebleeds))
d+geom_density(aes(group=arm, colour=arm, fill=arm),size=2, alpha=0.3)+facet_wrap(~country,nrow=2)+ggtitle("Nose Bleeds by country")

#bleed by previous hospitalizations by arm
wilcox.test(nosebleeds ~arm,data=dat[dat$hosp.bleed==2,]) #hospitalization lack of variability (majority=2)

#hospitalization due to nose bleed by country
kruskal.test(hosp.bleed ~ country, data = dat)   #non parametric test (kruskal) of between country hospitalization
d<-ggplot(dat,aes(x=hosp.bleed))
d+geom_density(aes(group=country, colour=country, fill=country), alpha=0.3)+facet_wrap(~country,ncol= nlevels(dat$country))+ggtitle("Hospitalizations by Country")

round(prop.table(table(dat$hosp.bleed,dat$country)*100),2)
aggregate(. ~ country,data = dat[,c("country","hosp.bleed")], FUN=function(x) c(mean =round(mean(x),2), SD=round(sd(x),2)))  #extract mean and SD

##########################
#Prediction Models
##########################
##bleed prediction\
mean(dat$nosebleeds)
var(dat$nosebleeds)   #overdispersed (mean<variance)
#negative binomial   
summary(nb <- glm.nb(nosebleeds ~  arm +offset(log(duration))+ hosp.bleed+tissue.use+mucus.viscosity+hosp.bleed+country
	+arm:mucus.viscosity+hosp.bleed:mucus.viscosity+arm:hosp.bleed
	,data = dat,link=log))
prednb<-exp(predict(nb,dat)) ##predict number of events
plot(dat$nosebleeds,prednb)
sqrt(sum((prednb-dat$nosebleeds)^2)/nrow(dat))
summary(prednb-dat$nosebleeds)

#zero-inflated binomial model
summary(zinb <- zeroinfl(nosebleeds ~ arm + offset(log(duration))+hosp.bleed+tissue.use+mucus.viscosity+country
	+arm:mucus.viscosity+hosp.bleed:mucus.viscosity
	|tissue.use
	,data = dat,dist="negbin"))
predzinb<-predict(zinb,dat,type="response")
plot(dat$nosebleeds,predzinb)
sqrt(sum((predzinb-dat$nosebleeds)^2)/nrow(dat))
vuong(update(zinb,.~1),zinb)  #test if intercept-only zero part model is sig. different than the zero model with tissue--answer is yes
vuong(zinb,nb)  #test if there is significant difference between standard negative binomial and zero-inflated model--answer is very marginally

#try standard poisson model 
summary(zip <- zeroinfl(nosebleeds ~ arm + offset(log(duration))+hosp.bleed+tissue.use+mucus.viscosity+country
	+arm:mucus.viscosity+hosp.bleed:mucus.viscosity
	|tissue.use, data = dat,dist="poisson"))
predzip<-predict(zip,dat,type="response")
sqrt(sum((predzip-dat$nosebleeds)^2)/nrow(dat))
#try zero inflated poisson model
summary(poi <- glm(nosebleeds ~ arm + offset(log(duration))+hosp.bleed+tissue.use+mucus.viscosity+country
	+arm:mucus.viscosity+hosp.bleed:mucus.viscosity
	, family = poisson, data = dat))
predpoi<-predict(poi,dat,type="response")
sqrt(sum((predpoi-dat$nosebleeds)^2)/nrow(dat))
vuong(zinb,zip)  #test if zero inflated NB is different than zero inflated Poisson--answer is no
vuong(nb,poi)  #test if regular negative binomial is better than poisson (due to overdispersion), answer is yes. Should choose regular negative binomial
vuong(zinb,poi)
vuong(nb,zip)
vuong(poi,zip)

####################
#Standard negative binomial model is selected as final GLM-based model for prediction
####################
#add nb predicted nose bleeds to result dataset
dat$nbPred<-prednb
head(dat,10)

##machine learning methods for prediction
attach(dat)
inTraining <- createDataPartition(nosebleeds, p = .75, list = FALSE)
training <- dat[ inTraining,-c(which(colnames(dat) %in% c("subject","nbPred")))]
testing  <- dat[-inTraining,-c(which(colnames(dat) %in% c("subject","nbPred")))]
str(training)
str(testing)
#random forest
rfFit <- train(nosebleeds~., 
		data=training,method = "rf",importance = TRUE,
		trControl=trainControl(method = "oob"))
rfFit$finalModel
plot(rfFit$finalModel)
rfPred<-predict(rfFit,testing)
plot(testing$nosebleeds,rfPred,xlim=c(0,max(testing$nosebleeds)),ylim=c(0,max(rfPred)),main="Random Forest")
RMSErf<-sqrt(sum((rfPred-testing$nosebleeds)^2)/nrow(testing))
RMSErf
plot(varImp(rfFit,scale=TRUE),main="Random Forest")

#gradient boost machine
gbmControl <- trainControl(## k-fold CV
                           method = "repeatedcv",
                           number = 5,
                           ## repeated x times
                           repeats = 5)
gbmGrid <-  expand.grid(interaction.depth = (1:10)*1,
                        n.trees = (1:30)*3,
                        shrinkage =0.1,
                        n.minobsinnode = 10)
gbmFit1 <- train(nosebleeds~., 
		data=training,
            method = "gbm",
           	trControl = gbmControl,
		tuneGrid=gbmGrid,
            verbose = FALSE
#		,preProc = c("center", "scale")
)
gbmFit1
gbmPred=predict(gbmFit1,testing)
RMSEgbm<-sqrt(sum((gbmPred-testing$nosebleeds)^2)/nrow(testing))
RMSEgbm
plot(testing$nosebleeds,gbmPred,xlim=c(0,max(testing$nosebleeds)),ylim=c(0,max(gbmPred)))
plot(varImp(gbmFit1,scale=TRUE),main="Gradient Bossting Machine")
ggplot(gbmFit1)


#svm
svmControl <- trainControl(## k-fold CV
                           method = "repeatedcv",
                           number = 5,
                           ## repeated x times
                           repeats = 5)
#
svmFit<- train(nosebleeds~.,
		data=training,
            method = "svmRadial",
            trControl = svmControl,
		tuneLength=10
		,preProc = c("center", "scale")
)
svmFit$finalModel
svmPred=predict(svmFit,testing)
RMSEsvm<-sqrt(sum((svmPred-testing$nosebleeds)^2)/nrow(testing))
RMSEsvm
plot(testing$nosebleeds,svmPred,xlim=c(0,max(testing$nosebleeds)),ylim=c(0,max(svmPred)))
plot(varImp(svmFit,scale=TRUE))
ggplot(svmFit)

##ensemble  of GBM, Random Forest 
combine<-as.data.frame(cbind(gbmPred,rfPred))
combineN<-cbind(gbmPred,rfPred)
cor(combineN)
combine$avg<-apply(combine,1,mean) 
sqrt(sum((combine$avg-testing$nosebleeds)^2)/nrow(testing))
plot(testing$nosebleeds,combine$avg,xlim=c(0,max(testing$nosebleeds)),ylim=c(0,max(combine$avg)))


#extract distribution (percentiles) of mucus viscosity
viscosity_q<-as.data.frame(quantile(dat$mucus.viscosity,seq(0,1,0.1)))

#From the Negative Binomial Model, calculate mean predicted treatment effect (reduction of bleeds) by percentile/cumulatively of viscosity
TxbyVisc<-data.frame()
for (i in 1:10) {
	TxbyVisc<-rbind(TxbyVisc,
		as.data.frame(cbind(paste0("Top ", i*10,"%"),as.data.frame(TukeyHSD(aov(nbPred~arm,data=dat[dat$mucus.viscosity<=viscosity_q[[1]][12-i],]))[[1]])))
	)	
}
names(TxbyVisc)[1]<-"Top x%"
plot(TxbyVisc[,1],TxbyVisc$diff,
    ylim=range(c(TxbyVisc$lwr, TxbyVisc$upr)),
    pch=19, xlab="Segment by Mucus Viscosity", ylab="Mean Predicted Reduction (Active-Placebo) ",
    main="Predicted Bleeds Reduction by Mucus Viscosity"
)
arrows(c(1:10), TxbyVisc$lwr, c(1:10), TxbyVisc$upr, length=0.05, angle=90, code=3)
abline(h=0,lwd=2,col="red")

#Repeat the plot for the Machine Learning-Predicted (ensemble) 
TxbyViscML<-data.frame()
for (i in 1:10) {
	TxbyViscML<-rbind(TxbyViscML,
		as.data.frame(cbind(paste0("Top ", i*10,"%"),as.data.frame(TukeyHSD(aov(nbPred~arm,data=dat[dat$mucus.viscosity<=viscosity_q[[1]][12-i],]))[[1]])))
	)	
}
names(TxbyViscML)[1]<-"Top x%"
plot(TxbyViscML[,1],TxbyViscML$diff,
    ylim=range(c(TxbyViscML$lwr, TxbyViscML$upr)),
    pch=19, xlab="Segment by Mucus Viscosity", ylab="Mean Predicted Reduction (Active-Placebo) ",
    main="Predicted Bleeds Reduction by Mucus Viscosity"
)
arrows(c(1:10), TxbyViscML$lwr, c(1:10), TxbyViscML$upr, length=0.05, angle=90, code=3)
abline(h=0,lwd=2,col="red")











