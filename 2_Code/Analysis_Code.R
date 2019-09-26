###############################################
### Code for analysis of FOCE dataset ###
### Authors:  D.K. Okamoto & D. Kline ###
### Updated 10.6.18                   ###
#########################################

library(lme4)
library(nlme)
library(dplyr)
library(ggplot2)

### set nonlinear parameters
k1 = 1e-4
k2 = 20
k3= 3

### read and format data ###
FOCE <- read.csv("1_Data/Experimental_Data.csv")%>%
  mutate(S1= ifelse(Status!="Live",1,0),
         S2= ifelse(Status=="Live",1,0),
         EID = factor(FOCE),
         k1=k1,
         k2=k2,
         k3=k3,
         TRT = factor(ifelse(Omega_AR<3,"TRT","CONTROL")),
         grp= factor(1),
         FOCE= ifelse(FOCE=="F1 ","F1",FOCE),
         color2= ifelse(TRT=="TRT",Color,0),
         color3= ifelse(Status!="Dead",Color,0))

FOCE <- FOCE%>%mutate(TRT_FAC = factor(TRT:Status))
levels(FOCE$TRT_FAC) <- c("Control\nDead","Control\nLive","Treatment\nDead","Treatment\nLive")
FOCE$TRT_FAC = factor(FOCE$TRT_FAC,levels(FOCE$TRT_FAC)[c(2,4,1,3)])

ggplot(aes(Change*100,x=TRT_FAC),data= FOCE)+
  geom_boxplot(size=0.75,outlier.shape= NA,aes(colour= TRT))+
  geom_point(alpha= 0.5,hjust= 0.1,aes(group= TRT_FAC,colour= TRT),position= position_nudge(x = -0.1))+
  stat_summary(geom= "point",aes(shape= TRT,fill= TRT),size= 3)+
  scale_colour_manual(values= c("#00AEEF","#EE2F24"))+
  scale_fill_manual(values= c("#00AEEF","#EE2F24"))+
  scale_shape_manual(values= c(21,22))+
  scale_y_continuous(breaks= c(seq(-20,15,by= 5)))+
  ylab("% Weight Change (200 day experiment)")+
  geom_hline(yintercept= 0)+
  theme(legend.position= "none",
        axis.title.x= element_blank())

### fit linear mixed model 
lme_mod <- lmer(I(Change*100)~Status*TRT+(1+Status+Status:TRT|Color)+(1|EID),data= FOCE,REML=FALSE)

cont <- lsmeans(lme_mod,pairwise ~ TRT*Status)


FOCE_summary <- FOCE%>%
  group_by(Status,TRT)%>%
  summarize(Omega_AR = mean(Omega_AR),
            Change = mean(Change))%>%
  data.frame()%>%
  mutate(est= summary(cont$lsmeans)$lsmean,
         SE = summary(cont$lsmeans)$SE)


### test using type II LR tests (because no interaction)
Anova(lme_mod,type= "II")
contrast(cont, list("Dead- Cntrl v Trt" = c(1,-1,0,0),
                 "Live - Cntrl v Trt" = c(0,0,1,-1),
                 "Cntrl - Live vs Dead" = c(1,0,-1,0),
                 "Trt - Live vs Dead"= c(0,1,0,-1)),
         options = list(adjust= "bonferroni"))

### check assumption of residual normality 
shapiro.test(resid(lme_mod))

### check for bias in residual against fitted values
ggplot(aes(predict(lme_mod),resid(lme_mod,type="pearson",scaled= TRUE)),data= FOCE)+
  geom_point()+
  stat_smooth()+
  xlab("fitted")+
  ylab("residuals (pearson, standardized)")

### fit nonlinear model alternatives ### 
### linear
fit1 <- nlme(Change~-S1*exp(d1+b1*log(Omega_AR))+S2*exp(b2)*(Omega_AR-d2)/(1/k1/exp(b2)+(Omega_AR-d2)), data= FOCE,
             fixed= list(b1+b2+d1+d2~1), random= d1+b2~1|Color, start= c(b1=1,d1=1,b2=1,d2=1))

### strong nonlinearity 
fit2 <- nlme(Change~-S1*exp(d1+b1*log(Omega_AR))+S2*exp(b2)*(Omega_AR-d2)/(1/k2/exp(b2)+(Omega_AR-d2)), data= FOCE,
             fixed= list(b1+b2+d1+d2~1), random= d1+b2~1|Color, start= c(b1=1,d1=1,b2=1,d2=1))

### weak nonlinearity
fit3 <- nlme(Change~-S1*exp(d1+b1*log(Omega_AR))+S2*exp(b2)*(Omega_AR-d2)/(1/k3/exp(b2)+(Omega_AR-d2)), data= FOCE,
             fixed= list(b1+b2+d1+d2~1), random= d1+b2~1|Color, start= c(b1=1,d1=1,b2=1,d2=1))

## show models are equivalent fits
AIC(fit1,fit2,fit3)

### generate preduction dataframe
FOCE2 <- expand.grid(Omega_AR= seq(2.2,3.5,by= 0.01),
                     Status= c("Live","Dead"))%>%
  mutate(S1= ifelse(Status!="Live",1,0),
         S2= ifelse(Status=="Live",1,0),
         k1=k1,
         k2=k2,
         k3=k3)

### plot the expectations ###
p1a <- ggplot(aes(Omega_AR,Change*100),data= FOCE)+
  geom_point(aes(colour= Status,shape= Color))+
  geom_pointrange(aes(y=est,ymin= (est-SE),ymax= (est+SE),colour= Status),data= FOCE_summary)+
  geom_line(aes(y=predict(fit1,FOCE2,level=0)*100,colour= Status,linetype= "linear (k -> 0)"),size= 0.5,data= FOCE2)+
  geom_line(aes(y=predict(fit3,FOCE2,level=0)*100,colour= Status,linetype= "shallow saturation (k= 3)"),size= 0.5,data= FOCE2)+
  geom_line(aes(y=predict(fit2,FOCE2,level=0)*100,colour= Status,linetype= "steep saturation (k= 20)"),size= 0.5,data= FOCE2)+
  ylab("net change (%)")+
  scale_x_continuous(breaks= seq(2.25,3.5,by= 0.25),
                     labels= c(2.25,"2.5",2.75,"3",3.25,"3.5"))+
  scale_colour_manual(values= c("#00AEEF","#EE2F24"),name="")+
  guides(shape= FALSE,linetype= FALSE)+
  scale_linetype_manual(values= c(3,2,1),name= "")+
  xlab(expression(paste(Omega[AR])))+
  theme(legend.position= c(0.2,0.9),
        legend.spacing = unit(-0.6, "cm"))

### generate preduction dataframe
FOCE3 <- expand.grid(Omega_AR= seq(2.2,3.5,length.out= 30),
                     S1= seq(0,1,length.out=30),
                     k= c(k1,k2,k3),
                     k1=k1,
                     k2=k2,
                     k3=k3)%>%
  mutate(S2= 1-S1,
         pred= NA,
         k4= factor(k,labels= c("steep (k=20)","shallow (k=3)","linear (k= 0.0001)")))

FOCE3$pred[FOCE3$k==k1] <- predict(fit1,subset(FOCE3,k==k1),level=0)
FOCE3$pred[FOCE3$k==k2] <- predict(fit2,subset(FOCE3,k==k2),level=0)
FOCE3$pred[FOCE3$k==k3] <- predict(fit3,subset(FOCE3,k==k3),level=0)


### plot the expectations ###
p1 <- ggplot(aes(Omega_AR,S2),data= FOCE3)+
  geom_tile(aes(fill= pred*100))+
  geom_contour(aes(z= pred*100,linetype= factor(k)),breaks= 0)+
  facet_wrap(~k,ncol=1)+
  theme(legend.position= "top")+
  guides(fill= guide_colorbar(title.position= "top"))+
  ylab("proportion of coral living")+
  scale_fill_gradient2(low=brewer.pal(5,"RdYlBu")[[1]],mid=brewer.pal(5,"RdYlBu")[[3]],high=brewer.pal(5,"RdYlBu")[[5]],limits= c(-12,10),midpoint= 0,name= "net change (%)")+
  coord_cartesian(expand=c(0,0))+
  xlab(expression(paste(Omega[AR])))+
  scale_x_continuous(breaks= seq(2.25,3.5,by= 0.25),
                     labels= c(2.25,"",2.75,"",3.25,""))+
  guides(linetype= FALSE)+
  scale_linetype_manual(values= c(1,2,3),name= "dissolution isocline")+
  theme(plot.title=element_text(hjust=0),
        strip.text=element_blank(),
        legend.key.width= unit(0.8,"cm"))


pdf(width= 12,height =3.5,file= "contours.pdf")
p4
dev.off()

p5 <- ggplot(aes(Omega_AR,S2),data= FOCE3)+
  geom_contour(aes(z= pred,linetype= k4),breaks= 0)+
  ylab("proportion of coral living")+
  theme(legend.position= c(0.4,0.8))+
  scale_x_continuous(breaks= seq(2.25,3.5,by= 0.25),
                     labels= c(2.25,"2.5",2.75,"3",3.25,"3.5"))+
  xlab(expression(paste(Omega[AR])))+
  scale_linetype_manual(values= c(1,2,3),name= "")


pdf(width= 6,height =6,file= "Figure_3.pdf")
plot_grid(plot_grid(p1a,p5,rel_widths= c(2,1),ncol=1,align= "hv"),
          p1,ncol=2,rel_widths= c(1,0.75))
dev.off()

FOCE4 <- expand.grid(k1=k1,
                     k2=k2,
                     k3=k3)%>%
  mutate(pred= NA,
         S1=NA,
         S2=NA,
         Omega_AR= NA)


f1 <- function(Live,k,b1,b2,d1,d2){
  f_sub <- function(x) {
    -(1-Live)*exp(d1+b1*log(x))+Live*exp(b2)*(x-d2)/(1/k/exp(b2)+(x-d2))
   }
  uniroot(f_sub,interval = c(2.25,5))$root
}

boot_df <- expand.grid(k= c(k1,k2,k3),Live=seq(0.5,1,by= 0.05))

boot_df <- cbind(boot_df,rbind(fixef(fit1),
  fixef(fit2),
  fixef(fit3)))

for(i in 1:nrow(boot_df)){
  boot_df$N0[i] <- f1(boot_df$Live[i],boot_df$k[i],boot_df$b1[i],boot_df$b2[i],boot_df$d1[i],boot_df$d2[i])
}

pars <- fixef(fit1)

boot_df <- expand.grid(k= c(k1,k2,k3),Live=seq(1,by= 0.05))

boot_df <- cbind(boot_df,rbind(fixef(fit1),
                               fixef(fit2),
                               fixef(fit3)))

write.csv(boot_df,"net_dissolution.csv")

