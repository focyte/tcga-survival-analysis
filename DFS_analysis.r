# Author: John @ Focyte
# Created: 22/04/2022
# Version: 1.2
# Description: This script performs survival analysis based on copy number estimation for a specified gene.
#              It organizes patients into groups and generates Kaplan-Meier plots and Cox proportional hazards model results.

# Load required packages
library(survival)
library(survminer)
library(stringr)
library(gridExtra)
library(dplyr)
library(data.table)
library(tidyr)

# Function to perform survival analysis
perform_survival_analysis <- function(clinical_data_file, sample_data_file, cna_data_file, gene_of_interest, output_file) {
  
  # Argument checks
  if (!file.exists(clinical_data_file) || !file.exists(sample_data_file) || !file.exists(cna_data_file)) {
    stop("Error: One or more input files do not exist.")
  }
  
  if (!is.character(gene_of_interest) || length(gene_of_interest) != 1) {
    stop("Error: 'gene_of_interest' must be a single character string.")
  }
  
  if (!is.character(output_file) || length(output_file) != 1) {
    stop("Error: 'output_file' must be a single character string.")
  }
  
  # Read clinical data
  df <- read.table(clinical_data_file, header = TRUE, sep = "\t")
  clinical_DFS <- data.frame(df$PATIENT_ID, df$DFS_STATUS, df$DFS_MONTHS)
  colnames(clinical_DFS) <- c("PATIENT_ID", "status", "time")
  clinical_DFS[clinical_DFS == "0:DiseaseFree"] <- "0"
  clinical_DFS[clinical_DFS == "1:Recurred/Progressed"] <- "1"
  clinical_DFS <- clinical_DFS %>% drop_na()
  
  # Read sample data
  ID <- read.table(sample_data_file, header = TRUE, sep = "\t")
  df <- data.frame(ID)
  IDS <- data.frame(df$PATIENT_ID, df$SAMPLE_ID)
  colnames(IDS) <- c("PATIENT_ID", "SAMPLE_ID")
  
  # Read CNA data
  CNA <- read.table(cna_data_file, header = FALSE, sep = "\t")
  df2 <- data.frame(CNA)
  df_t <- transpose(df2)
  names(df_t) <- as.matrix(df_t[1, ])
  df_t <- df_t[-1, ]
  df_t <- df_t[-1, ]
  
  # Merge CNA and Survival
  names(df_t)[names(df_t) == 'Hugo_Symbol'] <- 'SAMPLE_ID'
  CNV_patients <- merge(df_t, IDS, by = "SAMPLE_ID")
  CNV_clinical <- merge(CNV_patients, clinical_DFS, by = "PATIENT_ID")
  
  # Run a test using the specified gene
  gene <- CNV_clinical[[gene_of_interest]]
  test <- data.frame(CNV_clinical$PATIENT_ID, gene, CNV_clinical$time, CNV_clinical$status)
  colnames(test) <- c("PATIENT_ID", "gene", "time", "status")
  test <- transform(test, gene = as.numeric(gene))
  test <- transform(test, status = as.numeric(status))
  
  # Create gene copy number groups
  test$gene[test$gene == "0"] <- "1.Diploid"
  test$gene[test$gene == "1"] <- "4.Gain"
  test$gene[test$gene == "-1"] <- "2.Loss - Shallow"
  test$gene[test$gene == "-2"] <- "3.Loss - Deep"
  
  # Survival Object
  surv_object <- Surv(time = test$time, event = test$status)
  
  # Coxph analysis
  fit.coxph <- coxph(surv_object ~ gene, data = test)
  summary(fit.coxph)
  
  # Generate a Forest plot of the coxph analysis
  ggforest(fit.coxph, data = test)
  
  # KM Plot
  fit1 <- survfit(surv_object ~ gene, data = test)
  ggsurvplot(fit1, data = test, 
             title = paste("Effect of", gene_of_interest, "CNV on DFS in Prostate Cancer"),
             pval = TRUE, 
             conf.int = TRUE,
             risk.table.col = "strata",
             legend = "right",
             legend.title = paste(gene_of_interest, "CNV"),
             legend.labs = c("Diploid", "Shallow Loss", "Deep Loss", "Gain"),
             ggtheme = theme_bw(),
             axes.offset = TRUE,
             pval.method = TRUE,  
             risk.table = TRUE,
             palette = "jco", 
             tables.theme = theme_cleantable())
  
  # Save results to the output file
  write.table(summary(fit.coxph), file = output_file, append = TRUE, quote = FALSE, row.names = FALSE)
}
