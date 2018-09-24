#Import tidyverse

library("tidyverse")


#Dataset 1: Pharmacokinetics of Indomethacin
data(Indometh)

Inm<-as.tibble(Indometh)

#Create additional variables to calculate mean and standard deviation of drug concentrations across all six subjects at each timepoint
Inm%>%group_by(time)%>%mutate(mean_conc=mean(conc),sd_conc=sd(conc))->Inm_change

#Remove subject IDs and conc to only retain the mean and standard deviation of drug concentrations for all six subjects at each timepoint
Inm_change%>%select(-Subject,-conc)%>%distinct()->Inm_trunc

#Create plot of all subjects' data, with mean and standard deviations superimposed 

ggplot()+
  geom_point(data=Inm,mapping=aes(x=time,y=conc,colour=Subject),size=1.5)+scale_colour_hue(breaks=c("1","2","3","4","5","6"))+
  geom_point(data=Inm_trunc,mapping=aes(x=time,y=mean_conc),shape=21,fill="white")+
  geom_line(data=Inm_trunc,mapping=aes(x=time,y=mean_conc),size=2,alpha=0.1)+
  geom_errorbar(data=Inm_trunc,mapping=aes(x=time,ymin=mean_conc-sd_conc,ymax=mean_conc+sd_conc),width=0.25,size=1.5,alpha=0.1)+
  scale_y_continuous(name="Plasma Concentration (ug/mL)",breaks=seq(0,3,by=0.5))+
  scale_x_continuous(name="Time of Blood Draw (hr)",breaks=seq(0,8,by=1))+
  labs(title="Pharmacokinetics of Indomethacin in 6 Subjects over 8 Hours",
       caption = "Drug plasma concentrations for each subject are plotted as individual points.
       \n The line and white points are used to represent the mean plasma concentration at each time point. 
       \n The error bars represent a range of two standard deviations of the plasma concentrations.")+theme(plot.title=element_text(hjust=0.5))




#Dataset 2: Puromycin and Reaction Velocity of an Enzymatic Reaction
data(Puromycin)

Puro<-as.tibble(Puromycin)

#Create variables for 1/concentration and 1/rate and create separate tibbles (for cells treated and untreated with puromycin)

Puro%>%mutate(Inv_rate=1/rate,Inv_conc=1/conc)->Puro_mod

Puro_mod%>%filter(state=="treated")->Puro_LB_treated

Puro_mod%>%filter(state=="untreated")->Puro_LB_untreated

#Calculating the Michaelis-Menten constant and maximum reaction velocity.
#Y-intercept is 1/(maximum reaction velocity)
#X-intercept is -1/(Michaelis-Menten constant)

#For untreated

lm(Inv_rate~Inv_conc,Puro_LB_untreated)->lm_untreated

#summary(lm_untreated) #Gives summary statistics

#names(summary(lm_untreated)) #Find names of the elements I want

as.data.frame(summary(lm_untreated)$coefficients)$Estimate -> lm_untreated_int_slope

Vmax_untreated=round(1/(lm_untreated_int_slope[1]),digits=2)


Km_untreated=round(-1/(-lm_untreated_int_slope[1]/lm_untreated_int_slope[2]),digits=4)


untreated_eqn<-paste("y = ",formatC(lm_untreated_int_slope[2],format='e',digits=2)," x + ",formatC(lm_untreated_int_slope[1],format='e',digits=2),sep="")


#For treated 

lm(Inv_rate~Inv_conc,Puro_LB_treated)->lm_treated

#summary(lm_treated) #Gives summary statistics

#names(summary(lm_treated)) #Find names of the elements I want

as.data.frame(summary(lm_treated)$coefficients)$Estimate -> lm_treated_int_slope

Vmax_treated=round(1/(lm_treated_int_slope[1]),digits=2)

Km_treated=round(-1/(-lm_treated_int_slope[1]/lm_treated_int_slope[2]),digits=4)

treated_eqn<-paste("y = ",formatC(lm_treated_int_slope[2],format='e',digits=2)," x + ",formatC(lm_treated_int_slope[1],format='e',digits=2),sep="")



#Plotting
ggplot()+
  geom_point(data=Puro_mod,mapping=aes(x=Inv_conc,y=Inv_rate,colour=state))+
  scale_x_continuous(name="Reciprocal Substrate Concentration (1/ppm)",limits=c(-40,50),breaks=seq(-40,50,by=10))+
  scale_y_continuous(name="Reciprocal Reaction Rate (1/counts/min/min)",limits=c(-0.01,0.025),breaks=seq(-0.01,0.03,by=0.005))+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=c("red","blue"))+
  geom_smooth(data=Puro_LB_treated,mapping=aes(x=Inv_conc,y=Inv_rate),colour="red",method="lm",se=FALSE,fullrange=TRUE)+
  annotate('text',label=treated_eqn,x=30,y=0.005,colour="red")+
  geom_smooth(data=Puro_LB_untreated,mapping=aes(x=Inv_conc,y=Inv_rate),colour="blue",method="lm",se=FALSE,fullrange=TRUE)+
  annotate('text',label=untreated_eqn,x=20,y=0.017,colour="blue")+
  labs(title="Lineweaver-Burk Plot of an Enzymatic Reaction in Cells \n with or without Puromycin Treatment",
       caption = "Datapoints are grouped by treatment status.
        \n The reciprocals of the substrate concentrations and reaction rates are then plotted.
       \n The blue and red lines represent the Lineweaver-Burk equation for each treatment group. 
       \n The lines are extrapolated to graphically show important terms.")+theme(plot.title=element_text(hjust=0.5))





#Dataset 3: ELISA Assay of DNase

data(DNase)

dna<-as.tibble(DNase)

#Convert concentrations to log10 scale and create summary statistics for both all runs together and individual runs.

dna%>%mutate(log_conc=log10(conc))%>%group_by(Run,log_conc)%>%summarise(mean_density=mean(density),sd_density=sd(density))->summed_log_dna

dna%>%mutate(log_conc=log10(conc))%>%group_by(log_conc)%>%summarise(mean_density=mean(density),sd_density=sd(density))->summed_mean_log_dna

summed_log_dna%>%filter(log_conc >=-0.25, log_conc<=1.25)->linear_log_dna #Create another tibble with just the data that falls within a predicted linear dynamic range

lm(mean_density~log_conc,linear_log_dna)->lm_linear

as.data.frame(summary(lm_linear)$coefficients)$Estimate->lm_linear_int_slope

lin_range_eqn<-paste("y = ",formatC(lm_linear_int_slope[2],format='e',digits=2)," x + ",formatC(lm_linear_int_slope[1],format='e',digits=2),sep="")


ggplot()+
  geom_point(data=summed_mean_log_dna,mapping=aes(x=log_conc,y=mean_density),fill="brown",size=4,shape=23,alpha=0.3)+
  geom_point(data=summed_log_dna,mapping=aes(x=log_conc,y=mean_density,colour=Run))+
  scale_colour_hue(breaks=c("1","2","3","4","5","6","7","8","9","10","11"))+
  geom_line(data=summed_mean_log_dna,mapping=aes(x=log_conc,y=mean_density),colour="brown",size=2,alpha=0.3)+
  geom_errorbar(data=summed_mean_log_dna,mapping=aes(x=log_conc,ymin=mean_density-sd_density,ymax=mean_density+sd_density),colour="brown",width=0.1,size=1.5,alpha=0.2)+
  geom_smooth(data=linear_log_dna,mapping=aes(x=log_conc,y=mean_density),method="lm",se=FALSE,size=1,alpha=0.2)+
  geom_vline(xintercept=-0.25)+
  geom_vline(xintercept=1.25)+
  annotate('text',label="Reasonable dynamic \n linear range",x=0.5,y=0.25)+
  annotate('text',label=lin_range_eqn,x=0.25,y=1.5,colour="blue")+
  xlab('Protein Concentration (log10 Scale)')+
  ylab('Measured Optical Density')+
  labs(title="Standard Curve of ELISA Assay for Recombinant DNase in Rat Serum",
       caption = "Log10 Concentrations were first calculated for each protein sample.
       \n Averaged optical densities measured for each protein concentration were plotted for run. 
       \n Measured optical densities for the same concentration across all runs were then averaged and plotted.
       \n The error bars represent a range of two standard deviations of the optical densities.
       \n A dynamic linear range was assessed visually, and a linear regression was performed.")+
  theme(plot.title=element_text(hjust=0.5))

#Dataset 4:  The Effect of Vitamin C on Tooth Growth in Guinea Pigs

data(ToothGrowth)

tooth<-as.tibble(ToothGrowth)


#Create summary statistics for each combination of delibery method and dose to compare means to median graphically.

tooth%>%group_by(supp,dose)%>%summarise(mean_length=mean(len),median_length=median(len), IQR_half=IQR(len)/2,sd_length=sd(len))->summed_tooth



ggplot()+
  geom_boxplot(data=tooth,mapping=aes(x=as.factor(dose),y=len,colour=supp))+
  geom_point(data=summed_tooth,mapping=aes(x=as.factor(dose),y=mean_length,colour=supp),shape=17,size=2.5)+
  xlab('Vitamin C Dose (mg/day)')+
  ylab("Length of Odontoblasts")+
  labs(title="Effect of Vitamin C Dose and Delivery Method on Tooth Growth in Guinea Pigs",
       caption="Means of the corresponding datasets are plotted (as triangles) to compare to the median.")+
  theme(plot.title=element_text(hjust=0.5))


  
  
  
