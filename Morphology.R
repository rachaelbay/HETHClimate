setwd("~/Documents/IndividualProjects/Nicole/HETH/")

library(tidyverse)
#library(lme4)
library(viridis)
library(shades)
library(patchwork)

############
## All samples
############
dat <- read.csv("Data/Data_for_Rachael.csv")
spr <- dat %>% filter(Season=="Spring" & NewAge=="AHY" & !is.na(Sex))


pdf("Figures/YearHist.pdf",width=5,height=4)
ggplot(spr,aes(x=Year)) + geom_histogram(fill="grey80",color="grey20") + theme_classic()
dev.off()

##Color palette
cols <- viridis(10)
light <- lightness(cols, scalefac(1.3))
dark <- cols
m.cols <- light[3]
f.cols <- light[5]
a.cols <- alpha(c(f.cols,m.cols),alpha=0.7)



##Bill
B1 <- lm(Bill~Year*Sex,data=spr)
B2 <- lm(Bill~Year+Sex,data=spr)
anova(B1,B2)
B3 <- lm(Bill~Year,data=spr[!is.na(spr$Sex),])
anova(B2,B3) #No need for Sex, use B3
summary(B3)
bl <- ggplot(spr,aes(x=Year,y=Bill)) + 
  geom_jitter(size=2,aes(color=Sex)) + 
  theme_classic() + geom_smooth(color=viridis(1),method="lm") +
  scale_color_manual(labels=c("Female","Male"),values=a.cols) + ylab("Bill Length")

spr %>% filter(Year<1985) %>% summarize(mean(Bill))  
spr %>% filter(Year>2010) %>% summarize(mean(Bill,na.rm=T))  


##Tarsus
T1 <- lm(Tarsus~Year*Sex,data=spr)
T2 <- lm(Tarsus~Year+Sex,data=spr)
anova(T1,T2)
T3 <- lm(Tarsus~Year, data=spr[!is.na(spr$Sex),])
anova(T2,T3) # We do need sex use T2
summary(T2)
tr <- ggplot(spr,aes(x=Year,y=Tarsus,color=Sex)) + 
  geom_jitter(size=2) + 
  theme_classic() + geom_smooth(method="lm",show.legend=F) + 
  scale_color_manual(labels=c("Female","Male"),values=a.cols) +ylab("Tarsus Length")

##Wing
W1 <- lm(Wing~Year*Sex,data=spr)
W2 <- lm(Wing~Year+Sex,data=spr)
anova(W1,W2)
W3 <- lm(Wing~Year, data=spr[!is.na(spr$Sex),])
W4 <- lm(Wing~Sex, data=spr)
anova(W2,W4) # 
summary(W4)
wi <- ggplot(spr,aes(x=Year,y=Wing,color=Sex)) + 
  geom_jitter(size=2) + 
  theme_classic() + ylim(82,106) +
  scale_color_manual(labels=c("Female","Male"),values=a.cols) + ylab("Wing Length")

##Relative Bill
spr$relBill <- spr$Bill/spr$Tarsus
RB1 <- lm(relBill~Year*Sex,data=spr)
RB2 <- lm(relBill~Year+Sex,data=spr)
anova(RB1,RB2)
RB3 <- lm(relBill~Year,data=spr)
anova(RB2,RB3)
summary(RB2)
RB4 <- lm(relBill~Sex,data=spr)
rb <- ggplot(spr,aes(x=Year,y=relBill,color=Sex)) + 
  geom_jitter(size=2) + 
  theme_classic() + geom_smooth(method="lm",show.legend=F) + 
  scale_color_manual(labels=c("Female","Male"),values=a.cols) +ylab("Relative Bill Length")

##Relative Wing
spr$relWing <- spr$Wing/spr$Tarsus
RW1 <- lm(relWing~Year*Sex,data=spr)
RW2 <- lm(relWing~Year+Sex,data=spr)
anova(RW1,RW2)
RW3 <- lm(relWing~Year,data=spr)
anova(RW2,RW3)
summary(RW2)
rw <- ggplot(spr,aes(x=Year,y=relWing,color=Sex)) + 
  geom_jitter(size=2) + 
  theme_classic() + geom_smooth(method="lm",show.legend=F) + 
  scale_color_manual(labels=c("Female","Male"),values=a.cols) +ylab("Relative Wing Length")

###For supplementary figure
pdf("Figures/New_07.24/AllMorph.pdf",width=14,height=7)
tr + bl + wi +plot_spacer() + rb + rw + plot_annotation(tag_levels='A')  +plot_layout(guides="collect",nrow=2) 
dev.off()

##Just Bill for main figure
bl2 <- ggplot(spr,aes(x=Year,y=Bill)) + 
  geom_jitter(size=2,color=alpha("grey50",alpha=0.5)) + 
  theme_classic() + geom_smooth(color=viridis(1),method="lm") +
   ylab("Bill Length")

pdf("Figures/Bill.pdf",width=5,height=4)
bl2
dev.off()

