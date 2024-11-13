setwd("C:\\Users\\Lenovo\\Desktop\\corn wheat mixture analyses")
library(ggplot2)
library(gridExtra)
library(basicTrendline)
library(lmodel2) 
library(scales)
library(sqldf)
library(dplyr)
library(patchwork)
df=read.csv("cornwheatmix.csv")
df$total=df$rootDW+df$stemDW+df$leafDW
df$leafraction=(df$leafDW)/(df$total)

fig1F=ggplot(data = df,aes(x=log10(total),y=log10(leafDW),group=density,color=density))+geom_point(size=4)+geom_smooth(method="lm",span=1)+
  labs(title="Corn wheat mixture",colour="",x="Log total biomass (g)",y=expression(paste("Log leaf biomass (g)")))+
  annotate("text",x=-1.95, y=-0.1, label="F",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.8,0.3), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(-2,0,by=0.5),limits=c(-2,0))+
  scale_y_continuous(breaks=seq(-4,0,by=1),limits=c(-4,0))

result10=lmodel2(log10(df$leafDW)~log10(df$total))
#################
#leaf fraction vs ontogeny
dfLowdense=subset(df,df$density=="Low")
dfMiddledense=subset(df,df$density=="Middle")
dfHighdense=subset(df,df$density=="High")

MLowdense=sqldf("select repeated, IDtime,avg(leafraction) as leafraction_avg,avg(NOtime) as time_avg from dfLowdense group by repeated, IDtime")
MMiddledense=sqldf("select repeated, IDtime,avg(leafraction) as leafraction_avg,avg(NOtime) as time_avg from dfMiddledense group by repeated, IDtime")
MHighdense=sqldf("select repeated, IDtime,avg(leafraction) as leafraction_avg,avg(NOtime) as time_avg from dfHighdense group by repeated, IDtime")

MLowdense=MLowdense %>%  
  group_by(time_avg) %>%  summarise(sd = sd(leafraction_avg, na.rm = TRUE),MeanLT=mean(leafraction_avg, na.rm = TRUE))
MMiddledense=MMiddledense %>%  
  group_by(time_avg) %>%  summarise(sd = sd(leafraction_avg, na.rm = TRUE),MeanLT=mean(leafraction_avg, na.rm = TRUE)) 
MHighdense=MHighdense %>%  
  group_by(time_avg) %>%  summarise(sd = sd(leafraction_avg, na.rm = TRUE),MeanLT=mean(leafraction_avg, na.rm = TRUE)) 

NN=length(MLowdense$time_avg);class=c(rep("Low",NN),rep("Medium",NN),rep("High",NN))
dffig2A=data.frame(leafraction=c(MLowdense$MeanLT,MMiddledense$MeanLT,MHighdense$MeanLT),
                   SDLT=c(MLowdense$sd,MMiddledense$sd,MHighdense$sd),
                   time=rep(MLowdense$time_avg,3),class=class)

#leaf fraction vs ontogeny
fig2A=ggplot(data = dffig2A,aes(x=time,y=leafraction,group=class,color=class))+geom_point(size=2)+geom_line()+
  geom_errorbar(aes(x = time, ymin = leafraction - SDLT, ymax = leafraction + SDLT, color = class),    
                width = 0.5, size = 1)+
  labs(title="Corn wheat mixture",colour="",x="Time (days)",y=expression(paste("Mean leaf fraction")))+
  annotate("text",x=5, y=0.5, label="A",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.8,0.3), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(5,45,by=10),limits=c(5,45))+
  scale_y_continuous(breaks=seq(0,0.5,by=0.1),limits=c(0,0.5))
#################

df$Dense=(df$number)/(0.32*0.22)# calculate density per unit area (g/m2)
Meandf1=df %>% group_by(species,repeated,IDtime,density) %>%  
  summarise(Dense=mean(Dense, na.rm = TRUE),.groups="drop")  
Meandf2=df %>% group_by(repeated,IDtime,density) %>%  
  summarise(cv = sd(leafraction, na.rm = TRUE)/mean(leafraction, na.rm = TRUE),
            MeanTotal=mean(total, na.rm = TRUE),Dense=mean(Dense, na.rm = TRUE),.groups="drop") 
Meandf2$wholemass=(Meandf2$MeanTotal)*(Meandf2$Dense)
# CV in response to joint competition (both intra- and interspecific competition) denoted by whole biomass per unit area
fig5=ggplot(data=Meandf2,aes(x=wholemass,y=cv,group=density,color=density))+geom_point(size=3)+geom_smooth(method = "lm")+
  labs(title="",y=expression("CVLT"),x=expression("Whole biomass (g/" ~ m^{2}~")"))+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.85,0.2), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(1.2,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0.0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.7, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(0,900,300),limits=c(0,900))+
  scale_y_continuous(breaks=seq(0,1,0.25),limits=c(0,1))

#cv in response to q(interspecific competition)
Corn_df<- Meandf1 %>%  filter(species == "Corn") 
Wheat_df<- Meandf1 %>%  filter(species == "Wheat") 
meancvLow=q1Low=q2Low=meancvM=q1M=q2M=meancvH=q1H=q2H=NA;Kcorn=10000;Kwheat=15000;iL=iM=iH=1
for (i in 1:54){
  if (Wheat_df$density[i]=="Low"){
    meancvLow[iL]=Meandf2$cv[i]#(Corn_df$cv[i]+Wheat_df$cv[i])/2
    q1Low[iL]=(Kcorn-Corn_df$Dense[i])*Kwheat/(Wheat_df$Dense[i]*Kcorn) 
    q2Low[iL]=(Kwheat-Wheat_df$Dense[i])*Kcorn/(Corn_df$Dense[i]*Kwheat)
    iL=iL+1
  }
  else if (Wheat_df$density[i]=="Middle"){
    meancvM[iM]=Meandf2$cv[i]#(Corn_df$cv[i]+Wheat_df$cv[i])/2
    q1M[iM]=(Kcorn-Corn_df$Dense[i])*Kwheat/(Wheat_df$Dense[i]*Kcorn) 
    q2M[iM]=(Kwheat-Wheat_df$Dense[i])*Kcorn/(Corn_df$Dense[i]*Kwheat)
    iM=iM+1
  }
  else if (Wheat_df$density[i]=="High"){
    meancvH[iH]=Meandf2$cv[i]#(Corn_df$cv[i]+Wheat_df$cv[i])/2
    q1H[iH]=(Kcorn-Corn_df$Dense[i])*Kwheat/(Wheat_df$Dense[i]*Kcorn) 
    q2H[iH]=(Kwheat-Wheat_df$Dense[i])*Kcorn/(Corn_df$Dense[i]*Kwheat)
    iH=iH+1
  }
  
}
dfcvq=data.frame(meancvLow=meancvLow,meancvM=meancvM,meancvH=meancvH,
                 q1Low=q1Low,q1M=q1M,q1H=q1H,  q2Low=q2Low,q2M=q2M,q2H=q2H)
fig3A=ggplot(data=dfcvq,aes(x=q1M,y=meancvM))+geom_point(size=3)+geom_smooth(method = "lm")+
  annotate("text",x=4, y=0.6, label="A",colour="black",angle = 0,size=10 ,family="serif")+ 
  annotate("text",x=5.55, y=0.6, label="Medium density",colour="black",angle = 0,size=10 ,family="serif")+ 
  labs(title="Empirical observations",x="Interspecific competition (q1)",y="CVLT")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        #text = element_text(family="Times New Roman",color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.1,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(1.2,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0.0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.7, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(4,6,0.5),limits=c(4,6))+
  scale_y_continuous(breaks=seq(0.2,0.6,0.1),limits=c(0.2,0.6))

figS1A=ggplot(data=dfcvq,aes(x=q1Low,y=meancvLow))+geom_point(size=3)+geom_smooth(method = "lm")+
  labs(title="Low density",x="Interspecific competition (q1)",y="CVLT")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        #text = element_text(family="Times New Roman",color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.1,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(1.2,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0.0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.7, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(10,18,2),limits=c(10,18))+
  scale_y_continuous(breaks=seq(0,1,0.25),limits=c(0,1))

figS1B=ggplot(data=dfcvq,aes(x=q1H,y=meancvH))+geom_point(size=3)+geom_smooth(method = "lm")+
  labs(title="High density",x="Interspecific competition (q1)",y="CVLT")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        #text = element_text(family="Times New Roman",color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.1,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(1.2,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0.0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.7, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(1,1.4,0.1),limits=c(1,1.4))+
  scale_y_continuous(breaks=seq(0,1,0.25),limits=c(0,1))

figS1C=ggplot(data=dfcvq,aes(x=q2Low,y=meancvLow))+geom_point(size=3)+geom_smooth(method = "lm")+
  labs(title="Low density",x="Interspecific competition (q2)",y="CVLT")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        #text = element_text(family="Times New Roman",color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.1,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(1.2,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0.0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.7, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(6,12,2),limits=c(6,12))+
  scale_y_continuous(breaks=seq(0,1,0.25),limits=c(0,1))

figS1D=ggplot(data=dfcvq,aes(x=q2M,y=meancvM))+geom_point(size=3)+geom_smooth(method = "lm")+
  labs(title="Medium density",x="Interspecific competition (q2)",y="CVLT")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        #text = element_text(family="Times New Roman",color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.1,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(1.2,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0.0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.7, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(3,5,0.5),limits=c(3,5))+
  scale_y_continuous(breaks=seq(0.2,0.6,0.1),limits=c(0.2,0.6))

figS1E=ggplot(data=dfcvq,aes(x=q2H,y=meancvH))+geom_point(size=3)+geom_smooth(method = "lm")+
  labs(title="High density",x="Interspecific competition (q2)",y="CVLT")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        #text = element_text(family="Times New Roman",color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.1,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(1.2,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0.0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.7, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(1,1.3,0.1),limits=c(1,1.3))+
  scale_y_continuous(breaks=seq(0,1,0.25),limits=c(0,1))

result1=lmodel2(meancvLow~q1Low)
result2=lmodel2(meancvLow~q2Low)
result3=lmodel2(meancvM~q1M)
result4=lmodel2(meancvM~q2M)
result5=lmodel2(meancvH~q1H)
result6=lmodel2(meancvH~q2H)

L_df<- Meandf2 %>%  filter(density == "Low") 
M_df<- Meandf2 %>%  filter(density == "Middle") 
H_df<- Meandf2 %>%  filter(density == "High") 
result7=lmodel2(L_df$cv~L_df$wholemass)
result8=lmodel2(M_df$cv~M_df$wholemass)
result9=lmodel2(H_df$cv~H_df$wholemass)

fitresults=list(result1,result2,result3,result4,result5,result6,result7,result8,result9,result10)
OLSreg=SMAreg=matrix(nrow=10,ncol=9)
for (ij in 1:10){
  result=fitresults[[ij]]
  OLSreg[ij,]=c(result$n,result$rsquare,result$regression.results[1,3],result$confidence.intervals[1,4],
                result$confidence.intervals[1,5],result$regression.results[1,2],result$confidence.intervals[1,2],
                result$confidence.intervals[1,3],result$P.param)
  SMAreg[ij,]=c(result$n,result$rsquare,result$regression.results[3,3],result$confidence.intervals[3,4],
                result$confidence.intervals[3,5],result$regression.results[3,2],result$confidence.intervals[3,2],
                result$confidence.intervals[3,3],result$P.param)
}
#write.csv(OLSreg,"OLSreg CV vs q.csv")
#write.csv(SMAreg,"SMAreg CV vs q.csv")
#ggsave("fig1F.pdf",fig1F,width = 20, height = 20, units = "cm", dpi = 300) 
#ggsave("fig2A.pdf",fig2A,width = 20, height = 20, units = "cm", dpi = 300) 
#ggsave("fig3.pdf",(fig3A|transient1)/(transient2|Resil1),width = 40, height = 40, units = "cm", dpi = 300) 
#ggsave("fig5.pdf",fig5,width = 20, height = 20, units = "cm", dpi = 300) 
#ggsave("figS1.pdf",figS1A+figS1B+figS1C+figS1D+figS1E,width = 70, height = 40, units = "cm", dpi = 300) 

###################################
library(ggplot2)
library(gridExtra)
library(patchwork)
g=1;q2=1;p2=2;p1=2;q1=seq(0.1,0.9,0.05)
resilience1=-g*((p1-q2)*q1+(q2-p1)*p2)/(p1*p2-q1*q2)
g=2
resilience2=-g*((p1-q2)*q1+(q2-p1)*p2)/(p1*p2-q1*q2)
g=3
resilience3=-g*((p1-q2)*q1+(q2-p1)*p2)/(p1*p2-q1*q2)
resilience=c(resilience1,resilience2,resilience3)

qq1=c(q1,q1,q1)
class=c(rep("g=1",length(q1)),rep("g=2",length(q1)),rep("g=3",length(q1)))
df=data.frame(resilience=resilience,qq1=qq1,class=class)

Resil1=ggplot(data = df,aes(x=qq1,y=resilience,group=class,color=class))+geom_point(size=3)+scale_color_manual(values = c("red","blue","orange"))+
  labs(title="Theoretical prediction",colour="",x="Interspecific competition",y="Resilience")+
  annotate("text",x=0.0, y=1.8, label="D",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        #text = element_text(family="Times New Roman",color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.85,0.85), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 30,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(3,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(0.0,0.9,0.3),limits=c(0.0,0.9))+
  scale_y_continuous(breaks=seq(0,1.8,0.6),limits=c(0,1.8))


g=1;q2=2;p2=2;q1=0.1;p1=seq(0.25,1,0.01)
resilience1=-g*((p1-q2)*q1+(q2-p1)*p2)/(p1*p2-q1*q2)
g=2
resilience2=-g*((p1-q2)*q1+(q2-p1)*p2)/(p1*p2-q1*q2)
g=3
resilience3=-g*((p1-q2)*q1+(q2-p1)*p2)/(p1*p2-q1*q2)
resilience=c(resilience1,resilience2,resilience3)

pp1=c(p1,p1,p1)
class=c(rep("g=1",length(p1)),rep("g=2",length(p1)),rep("g=3",length(p1)))
df=data.frame(resilience=resilience,pp1=pp1,class=class)

Resil2=ggplot(data = df,aes(x=pp1,y=resilience,group=class,color=class))+geom_point(size=3)+scale_color_manual(values = c("red","blue","orange"))+
  labs(title="Theoretical prediction",colour="",x="Plant intraspecific competition",y="")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        #text = element_text(family="Times New Roman",color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.1,0.85), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 30,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(3,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(0.25,1,0.25),limits=c(0.25,1))+
  scale_y_continuous(breaks=seq(-9,0,3),limits=c(-9,0))

Resil1|Resil2

#ggsave("resilience.pdf",plot=grid.arrange(plot1,plot2,ncol=2),width = 40, height = 20, units = "cm",dpi=300)
####################################
g=1;q2=1;p2=2;p1=2;q1=seq(0.1,0.9,0.05)
r1=g*(1-q1/p2)
g=2
r2=g*(1-q1/p2)
g=3
r3=g*(1-q1/p2)

r=c(r1,r2,r3)
qq1=c(q1,q1,q1)
class=c(rep("g=1",length(q1)),rep("g=2",length(q1)),rep("g=3",length(q1)))
df=data.frame(r=r,qq1=qq1,class=class)

transient1=ggplot(data = df,aes(x=qq1,y=r,group=class,color=class))+geom_point(size=3)+scale_color_manual(values = c("red","blue","orange"))+
  labs(title="Theoretical prediction",colour="",x="Interspecific competition",y=expression("Transient metric"~r))+
  annotate("text",x=0, y=3, label="B",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        #text = element_text(family="Times New Roman",color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.85,0.85), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 30,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(3,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=c(0.0,0.3,0.6,0.9),limits=c(0.0,0.9))+
  scale_y_continuous(breaks=seq(0,3,1),limits=c(0,3))


t=1;
g=1
A=g-q1*g/p2;B=0;C=-q2*g/p2;D=-g
X=t*sqrt((A-D)^2+4*B*C)/2
K1=exp((A+D)*t/2)*(cosh(X)+(A-D)*sinh(X)/sqrt((A-D)^2+4*B*C))

g=2
A=g-q1*g/p2;B=0;C=-q2*g/p2;D=-g
X=t*sqrt((A-D)^2+4*B*C)/2
K2=exp((A+D)*t/2)*(cosh(X)+(A-D)*sinh(X)/sqrt((A-D)^2+4*B*C))

g=3
A=g-q1*g/p2;B=0;C=-q2*g/p2;D=-g
X=t*sqrt((A-D)^2+4*B*C)/2
K3=exp((A+D)*t/2)*(cosh(X)+(A-D)*sinh(X)/sqrt((A-D)^2+4*B*C))



rho=c(K1,K2,K3)
qq1=c(q1,q1,q1)
class=c(rep("g=1",length(q1)),rep("g=2",length(q1)),rep("g=3",length(q1)))

df=data.frame(rho=rho,qq1=qq1,class=class)
transient2=ggplot(data = df,aes(x=qq1,y=rho,group=class,color=class))+geom_point(size=3)+scale_color_manual(values = c("red","blue","orange"))+
  labs(title="Theoretical prediction",colour="",x="Interspecific competition",y=expression("Transient metric"~rho))+
  annotate("text",x=0.0, y=18, label="C",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        #text = element_text(family="Times New Roman",color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.position = c(0.85,0.85), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 30,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(3,'line'), # space between legend text
        legend.key = element_blank(),
        plot.title = element_text(size=30, hjust=0.4, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=c(0.0,0.3,0.6,0.9),limits=c(0.0,0.9))+
  scale_y_continuous(breaks=seq(0,18.0,6),limits=c(0,18))

transient1|transient2
#ggsave("transientmetric.pdf",plot=grid.arrange(plot1,plot2,plottfig4,ncol=3),width = 60, height = 20, units = "cm",dpi=300)
