getwd()
setwd("H:/ej index data/version1.0") 


################################################################
#Load data
cardio <- read.csv("Cardio.csv")
diesel <- read.csv("Diesel emission.csv")
education <- read.csv("Education.csv")
enviro_e <- read.csv("Environmental Effects.csv")
expo <- read.csv("Environmental Exposures.csv")
final_score <- read.csv("Environmental Health Disparities.csv")
housing <- read.csv("Housing burden.csv")
lead <- read.csv("Lead_risk.csv")
linguistic <- read.csv("Linguistic isolation.csv")
lbw <- read.csv("Low birth weight.csv")
npl <- read.csv("NPL sites.csv")
ozone <- read.csv("Ozone.csv")
pm25 <- read.csv("PM2.5.csv")
POC <- read.csv("POC.csv")
poverty <- read.csv("Poverty.csv")
rmp <- read.csv("RMP facility.csv")
sensitive_pop <- read.csv("Sensitive Populations.csv")
socioeco <- read.csv("Socioeconomic Factors.csv")
toxic_release <- read.csv("Toxic release into air.csv")
traffic <- read.csv("Traffic density.csv")
transportation <- read.csv("Transportation burden.csv")
haz_waste <- read.csv("TSDF.csv")
unemployment <- read.csv("Unemployment.csv")
wastewater <- read.csv("Wastewater discharge.csv")

################################################################
#merge files
#Merge Sensitive population
sensitive_pop <- merge(sensitive_pop,lbw,by="Census_Tract")
sensitive_pop <- merge(sensitive_pop,cardio,by="Census_Tract")
#write.csv(sensitive_pop,file="./merged/sensitive_pop.csv")

#Merge Socioeco
socioeco <- merge(socioeco,education,by="Census_Tract")
socioeco <- merge(socioeco,housing,by="Census_Tract")
socioeco <- merge(socioeco,linguistic,by="Census_Tract")
socioeco <- merge(socioeco,poverty,by="Census_Tract")
socioeco <- merge(socioeco,POC,by="Census_Tract")
socioeco <- merge(socioeco,unemployment,by="Census_Tract")
socioeco <- merge(socioeco,transportation,by="Census_Tract")
#write.csv(socioeco,file="./merged/socioeco.csv")

#Merge Exposure
expo <- merge(expo, diesel, by="Census_Tract")
expo <- merge(expo,ozone,by="Census_Tract")
expo <- merge(expo, pm25, by="Census_Tract")
expo <- merge(expo, toxic_release, by="Census_Tract")
expo <- merge(expo, traffic, by="Census_Tract")
#write.csv(expo,file="./merged/expo.csv")

#Merge Environmental effect
enviro_e <- merge(enviro_e,haz_waste,by="Census_Tract")
enviro_e <- merge(enviro_e,lead,by="Census_Tract")
enviro_e <- merge(enviro_e,npl,by="Census_Tract")
enviro_e <- merge(enviro_e,rmp,by="Census_Tract")
enviro_e <- merge(enviro_e,wastewater,by="Census_Tract")
#write.csv(enviro_e,file="./merged/enviro_e.csv")

#Merge all
pop_char <- merge(sensitive_pop,socioeco,by="Census_Tract")
poll_bur <- merge(expo,enviro_e,by="Census_Tract")
final_score <- merge(final_score,pop_char,by="Census_Tract")
final_score <- merge(final_score,poll_bur, by="Census_Tract")

#write.csv(final_score,file="./merged/final_score.csv")

final_score <- read.csv("./merged/final_score.csv")
################################################################
#Spearman's correlation
#Raw score
cor <- subset(final_score,select=c("Avg_PM25","Avg_Ozone","Toxic_release_air","Diesel_emission","Traffic_density_percent",
                                     "Lead_Percent_Risk","Avg_PWDIS","Avg_PNPL","Avg_PRMP","Avg_PTSDF",
                                     "Cardio_Mortality_Rate","LBW_Percent_of_live_births_modified",
                                     "Poverty_Percent","POC_Percent","Education_Percent","Linguistic_percent","Housing_Percent",
                                     "Unemployed_Percent","Transportation_percent"))

#Ranking
cor <- subset(df1,select=c("expo_IBL_rank","Ozone_IBL_Ranking","PM25_IBL_Ranking","Diesel_IBL_Rank","Toxic_Release_IBL_Ranking","Traffic_Density_IBL_Ranking",
                            "enviro_e_IBL_rank","Lead_IBL_Ranking", "NPL_IBL_Ranking","TSDF_IBL_Ranking","RMP_IBL_Ranking","WDIS_IBL_Ranking",
                            "sensitive_pop_IBL_rank", "Cardio_IBL_Ranking", "LBW_IBL_Ranking",
                            "socioeco_IBL_rank","Education_IBL_Ranking","Linguistic_IBL_Ranking","Poverty_IBL_Ranking",
                            "Unemployed_IBL_Ranking","Housing_IBL_Ranking","POC_IBL_Ranking","Transportation_IBL_Ranking","final_IBL_rank"))


correlation <- cor(cor,use="na.or.complete",method=c("spearman"))
#write.csv(correlation, file="./analysis/spear_ranking.csv")

################################################################


################################################################
#Detailed analysis
#df1 <- read.csv("./analysis/final_ruca.csv")

#Racial categories: White/Non-Hispanic, Black, American Indian/Alaskan Native, Asian, Native Hawaiian-Other Pacific Islander and Two or more races
#income categories: (??? 185% below Federal Poverty Level (FPL), between 185% below FPL and 185% above FPL; ??? above 185% FPL)



























############ Analysis
#merge RUCA code
library(dplyr)
library(readxl)
ruca <- read_excel("./analysis/RUCA2010.xls")
df1 <- merge(ruca,final_score,by="Census_Tract")
check <- anti_join(ruca,df1,by="Census_Tract")
View(check)
#confirmed 12 census tracts with 0 population 
rm(check)
#write.csv(df1,file="./analysis/final_ruca.csv")

####Mapping
library(rgdal)
library(maptools)
library(maps)
library(shapefiles)
library(RColorBrewer)

tracts <- readOGR(dsn="H:/ej index data/ej indicator data/sensitivity analysis/cb_2016_53_tract_500k/cb_2016_53_tract_500k.shp")
names(tracts)
tracts$Census_Tract <- tracts$ctcode

tracts <- merge(tracts,df1,by="Census_Tract",all.x=T)

plot(tracts) #map of Washington Census Tracts

#Separate out each census tract excluded by missing indicator
which_ozone <- which(tracts$final_IBL_rank=="10" & !(tracts$Ozone_IBL_Ranking=="10"))
which_pm <- which(tracts$final_IBL_rank=="10" & !(tracts$PM25_IBL_Ranking=="10"))
which_diesel <- which(tracts$final_IBL_rank=="10" & !(tracts$Diesel_IBL_Rank=="10"))
which_tox <- which(tracts$final_IBL_rank=="10" & !(tracts$Toxic_Release_IBL_Ranking=="10"))
which_traffic <- which(tracts$final_IBL_rank=="10" & !(tracts$Traffic_Density_IBL_Ranking=="10"))
which_npl <- which(tracts$final_IBL_rank=="10" & !(tracts$NPL_IBL_Ranking=="10"))
which_tsdf <- which(tracts$final_IBL_rank=="10" & !(tracts$TSDF_IBL_Ranking=="10"))
which_rmp <- which(tracts$final_IBL_rank=="10" & !(tracts$RMP_IBL_Ranking=="10"))
which_lead <- which(tracts$final_IBL_rank=="10" & !(tracts$Lead_IBL_Ranking=="10"))
which_wdis <- which(tracts$final_IBL_rank=="10" & !(tracts$WDIS_IBL_Ranking=="10"))
which_cardio <- which(tracts$final_IBL_rank=="10" & !(tracts$Cardio_IBL_Ranking=="10"))
which_lbw <- which(tracts$final_IBL_rank=="10" & !(tracts$LBW_IBL_Ranking=="10"))
which_edu <- which(tracts$final_IBL_rank=="10" & !(tracts$Education_IBL_Ranking=="10"))
which_ling <- which(tracts$final_IBL_rank=="10" & !(tracts$Linguistic_IBL_Ranking=="10"))
which_pov <- which(tracts$final_IBL_rank=="10" & !(tracts$Poverty_IBL_Ranking=="10"))
which_unemp <- which(tracts$final_IBL_rank=="10" & !(tracts$Unemployed_IBL_Ranking=="10"))
which_housing <- which(tracts$final_IBL_rank=="10" & !(tracts$Housing_IBL_Ranking=="10"))
which_transportation <- which(tracts$final_IBL_rank=="10" & !(tracts$Transportation_IBL_Ranking=="10"))
which_poc <- which(tracts$final_IBL_rank=="10" & !(tracts$POC_IBL_Ranking=="10"))

# highlight CT of interest  
#expo
plot(tracts, border = "gray50",main="Census tracts excluded by no ozone indicator")
plot(tracts[which_ozone, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no PM2.5 indicator")
plot(tracts[which_pm, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no diesel indicator")
plot(tracts[which_diesel, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no toxic release indicator")
plot(tracts[which_tox, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no traffic indicator")
plot(tracts[which_traffic, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no NPL indicator")
plot(tracts[which_npl, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no hazardous waste indicator")
plot(tracts[which_tsdf, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no Risk Management Plan indicator")
plot(tracts[which_rmp, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no lead indicator")
plot(tracts[which_lead, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no wastewater discharge indicator")
plot(tracts[which_wdis, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no cardiovascular disease indicator")
plot(tracts[which_cardio, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no low birth weight indicator")
plot(tracts[which_lbw, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no education indicator")
plot(tracts[which_edu, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no linguistic isolation indicator")
plot(tracts[which_ling, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no poverty indicator")
plot(tracts[which_pov, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no housing indicator")
plot(tracts[which_housing, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no transportation indicator")
plot(tracts[which_transportation, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no unemployment indicator")
plot(tracts[which_unemp, ], col = "green", add = T)

plot(tracts, border = "gray50",main="Census tracts excluded by no POC indicator")
plot(tracts[which_poc, ], col = "green", add = T)
################################################################
#PCA
df1 <- read.csv("./merged/final_score.csv")

library(psych)
library(ggplot2)
library(ggfortify)
library(data.table)
library(GPArotation)

#Decile data
data <- subset(df1,select=c("Census_Tract","Cardio_IBL_Ranking", "LBW_IBL_Ranking",
                            "Education_IBL_Ranking","Linguistic_IBL_Ranking","Poverty_IBL_Ranking",
                            "Unemployed_IBL_Ranking","Housing_IBL_Ranking","POC_IBL_Ranking","Transportation_IBL_Ranking",
                            "Ozone_IBL_Ranking","PM25_IBL_Ranking","Diesel_IBL_Rank","Toxic_Release_IBL_Ranking","Traffic_Density_IBL_Ranking",
                            "Lead_IBL_Ranking", "NPL_IBL_Ranking","TSDF_IBL_Ranking","RMP_IBL_Ranking","WDIS_IBL_Ranking",
                            "socioeco_IBL_rank", "sensitive_pop_IBL_rank", "expo_IBL_rank","enviro_e_IBL_rank",
                            "final_IBL_rank"))

colnames(data)[colnames(data)=="Cardio_IBL_Ranking"] <- "cardio"
colnames(data)[colnames(data)=="LBW_IBL_Ranking"] <- "LBW"
colnames(data)[colnames(data)=="Education_IBL_Ranking"] <- "low_education_attainment"
colnames(data)[colnames(data)=="Linguistic_IBL_Ranking"] <- "linguistic_isolation"
colnames(data)[colnames(data)=="Poverty_IBL_Ranking"] <- "poverty"
colnames(data)[colnames(data)=="Unemployed_IBL_Ranking"] <- "unemployment"
colnames(data)[colnames(data)=="Housing_IBL_Ranking"] <- "housing_burden"
colnames(data)[colnames(data)=="POC_IBL_Ranking"] <- "race_ethnicity"
colnames(data)[colnames(data)=="Transportation_IBL_Ranking"] <- "transportation_expense"
colnames(data)[colnames(data)=="Ozone_IBL_Ranking"] <- "ozone"
colnames(data)[colnames(data)=="PM25_IBL_Ranking"] <- "pm2.5"
colnames(data)[colnames(data)=="Diesel_IBL_Rank"] <- "diesel"
colnames(data)[colnames(data)=="Toxic_Release_IBL_Ranking"] <- "toxic_release"
colnames(data)[colnames(data)=="Traffic_Density_IBL_Ranking"] <- "traffic_density"
colnames(data)[colnames(data)=="Lead_IBL_Ranking"] <- "lead_risk"
colnames(data)[colnames(data)=="NPL_IBL_Ranking"] <- "Superfund_sites"
colnames(data)[colnames(data)=="TSDF_IBL_Ranking"] <- "hazardous_waste"
colnames(data)[colnames(data)=="RMP_IBL_Ranking"] <- "facilities_with_toxic_substances"
colnames(data)[colnames(data)=="WDIS_IBL_Ranking"] <- "wastewater_discharge"

#all 19 indicators
subset <- subset(data,select=c("cardio","LBW","low_education_attainment","linguistic_isolation","poverty","transportation_expense",
                               "unemployment","housing_burden","race_ethnicity","ozone","pm2.5",
                               "diesel","toxic_release","traffic_density","lead_risk","Superfund_sites",
                               "hazardous_waste","facilities_with_toxic_substances","wastewater_discharge"))

#18 indicators with LBW dropped
subset <- subset(data,select=c("cardio","low_education_attainment","linguistic_isolation","poverty","transportation_expense",
                               "unemployment","housing_burden","race_ethnicity","ozone","pm2.5",
                               "diesel","toxic_release","traffic_density","lead_risk","Superfund_sites",
                               "hazardous_waste","facilities_with_toxic_substances","wastewater_discharge"))

#17 indicators with cardio and LBW 
subset <- subset(data,select=c("Census_Tract","low_education_attainment","linguistic_isolation","poverty","transportation_expense",
                               "unemployment","housing_burden","race_ethnicity","ozone","pm2.5",
                               "diesel","toxic_release","traffic_density","lead_risk","Superfund_sites",
                               "hazardous_waste","facilities_with_toxic_substances","wastewater_discharge"))


#subset for pollution burden
subset <- subset(data,select=c("ozone","pm2.5","diesel","toxic_release","traffic_density","lead_risk",
                                    "Superfund_sites","hazardous_waste","facilities_with_toxic_substances","wastewater_discharge"))

#subset for population characteristics
subset <- subset(data,select=c("low_education_attainment","linguistic_isolation","poverty","transportation_expense",
                               "unemployment","housing_burden","race_ethnicity"))


#prcomp code PCA#
pca <- prcomp(na.omit(subset),scale=TRUE,rank.=5)
pca
summary(pca)
str(pca)
plot(pca,type="lines")
autoplot(pcaresults5,loadings=TRUE,loadings.label=TRUE,colour="gray")
#autoplot only shows pc1 (29% of variance) and pc2 (15% of var)



#principal code PCA
#Factor 3 analysis
pcaresults3 = principal(na.omit(subset),nfactors = 3)

zero.vals<-predict(pcaresults3,rep(0, length(pca.vars)),subset[,pca.vars])

adj.scores <- pcaresults3$scores
adj.scores[,1]<-pcaresults3$scores[,1]-zero.vals[1]
adj.scores[,2]<-pcaresults3$scores[,2]-zero.vals[2]
adj.scores[,3]<-pcaresults3$scores[,3]-zero.vals[3]

weightspca = data.table(pcaresults3$weights)
weightspca$vars = rownames(pcaresults3$weights)
#weightspca[vars=="pnc_background", vars := "PNCback"]
weightspca =  melt(weightspca, id = "vars")

facet_names <- c(
  `RC1` = paste0("RC1 ", "(Variance Accounted ",
                 round(pcaresults3$Vaccounted[2,1],4)*100,
                 "%)"),
  
  `RC2` = paste0("RC2", "(Variance Accounted ",
                 round(pcaresults3$Vaccounted[2,2],4)*100,
                 "%)"),
  `RC3` = paste0("RC3 ", "(Variance Accounted ",
                 round(pcaresults3$Vaccounted[2,3], 4)*100,
                 "%)")
)


p3 = ggplot(weightspca) + geom_bar(aes(vars, value), stat="identity") +
  theme_light(18) +
  ylab("Weights") +
  xlab("indicators") +
  facet_wrap(~variable,ncol=1, labeller = as_labeller(facet_names)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p3
fa.plot(pcaresults3,labels=col,show.points=TRUE,pos=1)
cluster.plot(pcaresults3)
fa.diagram(pcaresults3)

#Factors 4
pcaresults4 = principal(na.omit(subset),rotate="varimax",nfactors = 4)

zero.vals<-predict(pcaresults4,rep(0, length(pca.vars)),subset[,pca.vars])

adj.scores <- pcaresults4$scores
adj.scores[,1]<-pcaresults4$scores[,1]-zero.vals[1]
adj.scores[,2]<-pcaresults4$scores[,2]-zero.vals[2]
adj.scores[,3]<-pcaresults4$scores[,3]-zero.vals[3]
adj.scores[,4]<-pcaresults4$scores[,4]-zero.vals[4]

weightspca = data.table(pcaresults4$weights)
weightspca$vars = rownames(pcaresults4$weights)
weightspca =  melt(weightspca, id = "vars")

facet_names <- c(
  `RC1` = paste0("RC1 ", "(Variance Accounted ",
                 round(pcaresults4$Vaccounted[2,1],4)*100,
                 "%)"),
  
  `RC2` = paste0("RC2", "(Variance Accounted ",
                 round(pcaresults4$Vaccounted[2,2],4)*100,
                 "%)"),
  `RC3` = paste0("RC3 ", "(Variance Accounted ",
                 round(pcaresults4$Vaccounted[2,3], 4)*100,
                 "%)"),
  `RC4` = paste0("RC4 ", "(Variance Accounted ",
                 round(pcaresults4$Vaccounted[2,4], 4)*100,
                 "%)")
)

p4 = ggplot(weightspca) + geom_bar(aes(vars, value), stat="identity") +
  theme_light(18) +
  ylab("Weights") +
  xlab("indicators") +
  facet_wrap(~variable,ncol=1, labeller = as_labeller(facet_names)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p4
fa.plot(pcaresults4,labels=col,show.points=TRUE,pos=1)
cluster.plot(pcaresults4)
fa.diagram(pcaresults4)

#Factors 5
pca.vars <- colnames(subset)
fa.parallel(subset[,pca.vars]) #parallel analysis for PCA and FA

#pcaresults5 = principal(na.omit(subset),rotate="quartimax",nfactors = 5)
##################ES method of PCA###############################
#pcaresult5 <- prcomp(~ diesel + ozone + pm2.5 + traffic_density +
#                      lead_risk + hazardous_waste + Superfund_sites + facilities_with_toxic_substances + toxic_release + wastewater_discharge + 
#                      low_education_attainment + race_ethnicity + linguistic_isolation + poverty + unemployment + transportation_expense + housing_burden +
#                      cardio + LBW, data=subset)
#exclude cardio and lbw

pcaresults5 <- prcomp(~ diesel + ozone + pm2.5 + traffic_density +
                        lead_risk + hazardous_waste + Superfund_sites + facilities_with_toxic_substances + toxic_release + wastewater_discharge + 
                        low_education_attainment + race_ethnicity + linguistic_isolation + poverty + unemployment + transportation_expense + housing_burden,
                        data=subset,rank=5)


pcaresults5
screeplot(pcaresults5, npcs=length(pcaresults5$sdev))
summary(pcaresults5)

#merging data for pca predictions
pcadata <- subset
pcapredictions <- predict(pcaresults5, newdata=pcadata)
mergeall <- cbind(pcadata, pcapredictions)
#write.csv(mergeall, file = "pcaresults.csv")
#to plot PCA results on map
pcamap <- merge(tracts,mergeall,by="Census_Tract",all.x=T)

#PC1
spplot(pcamap,"PC1",colorkey=list(space="bottom"),scales=list(draw=TRUE),
       main="Urbanized areas",cuts=8,par.settings = list(axis.line = list(col = "transparent")),
       col.regions=brewer.pal(9,"BuGn"),col=NA)

#PC2
spplot(pcamap,"PC2",colorkey=list(space="bottom"),scales=list(draw=TRUE),
       main="SES",cuts=8,par.settings = list(axis.line = list(col = "transparent")),
       col.regions=brewer.pal(9,"BuGn"),col=NA)

#PC3
spplot(pcamap,"PC3",colorkey=list(space="bottom"),scales=list(draw=TRUE),
       main="Traffic-related pollution",cuts=8,par.settings = list(axis.line = list(col = "transparent")),
       col.regions=brewer.pal(9,"BuGn"),col=NA)

#PC4
spplot(pcamap,"PC4",colorkey=list(space="bottom"),scales=list(draw=TRUE),
       main="PC4",cuts=8,par.settings = list(axis.line = list(col = "transparent")),
       col.regions=brewer.pal(9,"BuGn"),col=NA)

#PC5
spplot(pcamap,"PC5",colorkey=list(space="bottom"),scales=list(draw=TRUE),
       main="PC5",cuts=8,par.settings = list(axis.line = list(col = "transparent")),
       col.regions=brewer.pal(9,"BuGn"),col=NA)

####################Principal method factor 5###############################
zero.vals<-predict(pcaresults5,rep(0, length(pca.vars)),subset[,pca.vars])

adj.scores <- pcaresults5$scores
adj.scores[,1]<-pcaresults5$scores[,1]-zero.vals[1]
adj.scores[,2]<-pcaresults5$scores[,2]-zero.vals[2]
adj.scores[,3]<-pcaresults5$scores[,3]-zero.vals[3]
adj.scores[,4]<-pcaresults5$scores[,4]-zero.vals[4]
adj.scores[,5]<-pcaresults5$scores[,5]-zero.vals[5]

weightspca = data.table(pcaresults5$weights)
weightspca$vars = rownames(pcaresults5$weights)
weightspca =  melt(weightspca, id = "vars")
weightspca

facet_names <- c(
  `RC1` = paste0("RC1 ", "(Variance Accounted ",
                 round(pcaresults5$Vaccounted[2,1],4)*100,
                 "%)"),
  
  `RC2` = paste0("RC2", "(Variance Accounted ",
                 round(pcaresults5$Vaccounted[2,2],4)*100,
                 "%)"),
  `RC3` = paste0("RC3 ", "(Variance Accounted ",
                 round(pcaresults5$Vaccounted[2,3], 4)*100,
                 "%)"),
  `RC4` = paste0("RC4 ", "(Variance Accounted ",
                 round(pcaresults5$Vaccounted[2,4], 4)*100,
                 "%)"),
  `RC5` = paste0("RC5 ", "(Variance Accounted ",
                 round(pcaresults5$Vaccounted[2,5], 4)*100,
                 "%)")
  
)

p5 = ggplot(weightspca) + geom_bar(aes(vars, value), stat="identity") +
  theme_light(18) +
  ylab("Weights") +
  xlab("indicators") +
  facet_wrap(~variable,ncol=1, labeller = as_labeller(facet_names)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p5

facet_names_assign <- c(
  `RC1` = paste0("Urban ", "(Variance Accounted ",
                 round(pcaresults5$Vaccounted[2,1],4)*100,
                 "%)"),
  
  `RC2` = paste0("Socioeconomics status ", "(Variance Accounted ",
                 round(pcaresults5$Vaccounted[2,2],4)*100,
                 "%)"),
  `RC3` = paste0("Traffic  ", "(Variance Accounted ",
                 round(pcaresults5$Vaccounted[2,3], 4)*100,
                 "%)"),
  `RC4` = paste0("Hazardous waste ", "(Variance Accounted ",
                 round(pcaresults5$Vaccounted[2,4], 4)*100,
                 "%)"),
  `RC5` = paste0("Peri-urban Superfund ", "(Variance Accounted ",
                 round(pcaresults5$Vaccounted[2,5], 4)*100,
                 "%)")
  
)

p5_assign = ggplot(weightspca) + geom_bar(aes(vars, value), stat="identity") +
  theme_light(18) +
  ylab("Weights") +
  xlab("indicators") +
  facet_wrap(~variable,ncol=1, labeller = as_labeller(facet_names_assign)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p5_assign

fa.plot(pcaresults5,labels=col,show.points=TRUE,pos=1)
cluster.plot(pcaresults5)
fa.diagram(pcaresults5)

################################################################
#Plot
df1$final_EHD <- as.numeric(df1$final_IBL_rank)
df1$socio_rank <- as.factor(df1$socioeco_IBL_rank)
socio_final <- ggplot(df1, aes(x=socio_rank,y=final_EHD)) + geom_boxplot()
socio_final

df1$race_rank <- as.factor(df1$POC_IBL_Ranking)
race_final <- ggplot(df1,aes(race_rank,final_EHD)) + geom_boxplot()
race_final

boxplot(final_EHD~race_rank, data=df1, xlab="Race (POC) ranking", ylab="EHD Rank")

df1$income_rank <- as.factor(df1$Poverty_IBL_Ranking)
income_final <- ggplot(df1,aes(income_rank, final_EHD)) + geom_boxplot()
income_final
boxplot(final_EHD~income_rank, data=df1, xlab="Income ranking", ylab="EHD Rank")


df1$housing_rank <- as.factor(df1$Housing_IBL_Ranking)
ggplot(df1,aes(housing_rank, final_EHD)) + geom_boxplot()

df1$cardio_rank <- as.factor(df1$Cardio_IBL_Ranking)
cardio_rank <- ggplot(df1,aes(cardio_rank, final_EHD)) + geom_boxplot()
cardio_rank

df1$lbw_rank <- as.factor(df1$LBW_IBL_Ranking)
lbw_rank <- ggplot(df1,aes(lbw_rank, final_EHD)) + geom_boxplot()
lbw_rank

ruca_final <- ggplot(df1,aes(Scheme1, final_EHD)) + geom_boxplot()
ruca_final



