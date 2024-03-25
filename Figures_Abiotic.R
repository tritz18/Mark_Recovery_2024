library(tidyverse)
library(readxl)
library(writexl)
library(patchwork)

#### Set working directory and files 
setwd("/Users/thorn/OneDrive/Desktop/Ritz_YOY_NP_2023/Abiotic/")
f = list.files(pattern="*.xlsx")

#### Clean logger data from folder. Collects all logger files
#### and defines different variables in correct r format. 
#### Filters data to select stocking to catch ####
Abiotic_Raw <- purrr::map_df(f, function(x) {
  mydata <- read_excel(x)
  mydata$Date_Time <- as.POSIXct(mydata$Date_Time,  format="%Y/%m/%d %H:%M")
  mydata$Location<- as.factor(mydata$Location)
  mydata$Wetland<- as.factor(mydata$Wetland)
  mydata$Type <-as.factor(mydata$Type)
  mydata |>
    filter(Date_Time < "2022-07-01" & Date_Time > "2022-04-29") |>
    mutate(Date = as.Date(Date_Time)) |> 
    select(-Date_Time)
})

Abiotic<- Abiotic_Raw |>
  filter(Location %in% c("CHIP_REF3", "CHIP_REF4", "CHIP_SP2",
         "CM_REF4", "CM_REF5", "CM_SP1", "CM_SP2",
         "FC_REF7", "FC_SP5", "FC_SP6", 
         "PV_REF1", "PV_REF2", "PV_SP4", "PV_SP5"))


TEMP<- ggplot(Abiotic, aes(Wetland, TEMP))+
  stat_boxplot(geom="errorbar")+
         geom_boxplot(fill="#999999")+
  theme_bw()+
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=14),
  axis.title.x= element_blank())+
  scale_y_continuous(breaks=seq(4,30,6), limits=c(3.9, 31))+
  ylab("Water Temperature (C°)") 

TEMP

DO<- ggplot(Abiotic, aes(Wetland, DO))+
  stat_boxplot(geom="errorbar")+
  geom_boxplot(fill="#999999")+
  theme_bw()+
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=14),
        axis.title.x= element_blank())+
  scale_y_continuous(breaks=seq(0,20,4), limits=c(0,20))+
  ylab("Dissolved Oxygen (mg/l)") 


Figure_4<- TEMP/DO

ggsave("Figure_4.png", dpi = 300)






