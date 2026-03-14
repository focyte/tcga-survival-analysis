# Load required packages
library("survival")
library("survminer")
library("stringr")
library("gridExtra")
library("dplyr")
library("data.table")
library("tidyr")

#import clinical data as a table
clinical <- read.table("data_clinical_patient.txt", header=T, sep="\t")

#convert to dataframe
df <- data.frame(clinical)

#choose relevant columns for survival analysis
clinical_OS <- data.frame(df$PATIENT_ID, df$OS_STATUS, df$OS_MONTHS)

#new headers
colnames(clinical_OS) <-c ("PATIENT_ID", "status", "time")

#convert OS_STATUS to staus 0/1
clinical_OS[clinical_OS == "0:LIVING"] <- "0"
clinical_OS[clinical_OS == "1:DECEASED"] <- "1"

#remove rows where data is missing (NA)
clinical_OS <-clinical_OS %>% drop_na()

#output file for later use
write.csv(clinical_OS, "clinical_OS.csv")

#make list of patient ID and sample ID
ID <-read.table("data_clinical_sample.txt", header=T, sep="\t")
df <-data.frame(ID)
IDS <-data.frame(df$PATIENT_ID, df$SAMPLE_ID)
colnames(IDS) <-c("PATIENT_ID", "SAMPLE_ID")

#import CNA data as a table
CNA <-read.table("data_cna.txt", header=F, sep="\t")
df2 <-data.frame(CNA)

#convert to dataframe and transpose
df_t <- transpose(df2)

#rename the headers using the gene names
names(df_t) <- as.matrix(df_t[1, ])
df_t <- df_t[-1, ]
df_t <- df_t[-1, ]

names(df_t)[names(df_t) == 'Hugo_Symbol'] <- 'SAMPLE_ID'
CNV_patients<-merge(df_t, IDS, by="SAMPLE_ID")
CNV_clinical<-merge(CNV_patients, clinical_OS, by="PATIENT_ID")

#Run a test using one gene
test<-data.frame(CNV_clinical$PATIENT_ID, CNV_clinical$PTEN, CNV_clinical$time, CNV_clinical$status)

#new headers
colnames(test) <-c ("PATIENT_ID", "PTEN", "time", "status")
test <- transform(
  test,PTEN = as.numeric(PTEN))
test <- transform(
  test,status = as.numeric(status))

test$PTEN[test$PTEN == "0"] <- "Diploid"
test$PTEN[test$PTEN == "1"] <- "Gain"
test$PTEN[test$PTEN == "-1"] <- "Shallow Loss"
test$PTEN[test$PTEN == "-2"] <- "Deep Loss"

# Create survival object
surv_object <- Surv(time = test$time, event = test$status)

# coxPH analysis with Forest Plot
fit.coxph <- coxph(surv_object ~ PTEN, 
                   data = test)

summary(fit.coxph)

ggforest(fit.coxph, data = test)

# KM analysis
fit1 <- survfit(surv_object ~ PTEN, data = test)

ggsurvplot(fit1, data = test, 
           pval = TRUE, 
           conf.int = TRUE,
           risk.table.col = "strata",
           ggtheme = theme_bw() ,
           pval.method=TRUE,  
           risk.table = TRUE, 
           palette = "jco", 
           tables.theme = theme_cleantable())
