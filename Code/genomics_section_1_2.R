library(carData)
library(car)
library(survival)
library(ggplot2)
library(survminer)
library(dplyr)  

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Q1<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
setwd("~/Desktop/Functional_genomics/Functional_genomics_report")
# Load DNAm age
load('LUAD-7.Rda')

# Creation of a dataframe of estimated DNAm age
DANm_barcode <- names(dm.age.subset)
DNAm_age <- as.vector(dm.age.subset)
DNAm <- data.frame(Barcode = DANm_barcode, DNAm_age = DNAm_age)

# Add patient number column and tumor types

barcode_split<- strsplit(as.character(DNAm$Barcode), "[.]")
DNAm$patient_number  <- sapply(barcode_split,"[",3) # Add the patient number column 
DNAm$tumor_type  <- sapply(barcode_split,"[",4) # T for tumor, N for normal
DNAm$tumor_type <- gsub('0..', 'T', DNAm$tumor_type)
DNAm$tumor_type <- gsub('1..', 'N', DNAm$tumor_type)


# Load clinical data

clinic_data <- read.table('/Users/chloe/Desktop/Functional_genomics/Functional_genomics_report/LUAD.Clinical_Pick_Tier1/LUAD.clin.merged.picked.txt', header = FALSE, row.names = 1, sep = '\t')
clinic_data <- as.data.frame(t(clinic_data))

# Add patient number column 
barcode_split<- strsplit(as.character(clinic_data$`Hybridization REF`), "[-]")
clinic_data$patient_number <- sapply(barcode_split,"[",3) 

# ------------creation of final tumor set and normal set ------------
merge <-merge(DNAm, clinic_data, by ="patient_number") # merge the clinical data with the DNA methylation age
#merge <- df.merge[!is.na(df.merge$years_to_birth),] # Deleting row for patient that we don't have the chronological age
merge$years_to_birth <- as.numeric(as.character(merge$years_to_birth)) # transform to numeric
tumor <- subset(merge, tumor_type=="T")
normal <- subset(merge, tumor_type=="N")
sum(!is.na(normal$years_to_birth))

# Q1.1 Spearman correlation for cancers
cor.test(tumor$years_to_birth, tumor$DNAm_age, method = "spearman") # weak positive correlation (0.209) + significant
linear_model <- lm(c(tumor$DNAm_age) ~ c(tumor$years_to_birth))

summary(linear_model)

plot(tumor$years_to_birth, tumor$DNAm_age, main="DNAm age vs chronological age for cancer tissues",
     xlab="Chronological age", ylab="DNAm age", pch=1, col = "black")
abline(linear_model, col="red")
abline(c(0,1),col="blue") 
legend("topleft", legend=c("Regression line", "Identity line", "Data"),
       col=c("red", "blue", "black"), lty= c(1,1, NA), pch = c(NA,NA, 1), cex=0.8)

table(tumor$years_to_birth)
# Q1.2 Spearman correlation for normal tissues
cor.test(normal$DNAm_age, normal$years_to_birth, method ="spearman") #strong positive correlation (0.757) + significant 
plot(normal$years_to_birth, normal$DNAm_age, main=" DNAm age vs chronological age for normal tissues",
     xlab="Chronological age", ylab="DNAm age", pch=1, col = "black")
linear_model <- lm(c(normal$DNAm_age) ~ c(normal$years_to_birth))
abline(linear_model, col="red")
abline(c(0,1),col="blue") 
legend("topleft", legend=c("Regression line", "Identity line", "Data"),
       col=c("red", "blue", "black"), lty= c(1,1, NA), pch = c(NA,NA, 1), cex=0.8)

# Q1.3 Spearman correlation DNAm age of tumors and their patient matched normal tissue
# rename disctinct colum to merge by patient number 
names(tumor)[names(tumor) == "DNAm_age"] <- "DNAm_age_tumor"
names(normal)[names(normal) == "DNAm_age"] <- "DNAm_age_normal"
df.patients <- merge(tumor, normal, by ="patient_number")
cor.test(df.patients$DNAm_age_tumor, df.patients$DNAm_age_normal, method = "spearman") # p-value > 0.05 (not significant), on ne peut pas rejeter l'hypothèse null)

# Q1.4

#------one option to calculate acceleration (not choose)------
#plot(df.patients$DNAm_age_tumor/df.patients$DNAm_age_normal, df.patients$years_to_birth.x, main="Tumor age accelaration vs chronological age",
#     xlab="DNAm age tumor/DNAm age normal ", ylab="Chronological age", pch=19, col = "red") # few data only 58 observations
#cor.test(df.patients$DNAm_age_tumor/df.patients$DNAm_age_normal, df.patients$years_to_birth.x, method = "spearman") # p-value > 0.05 (not significant)

# ------ Second option (not choose)------
#plot(tumor$DNAm_age_tumor/tumor$years_to_birth, tumor$years_to_birth, main="Tumor age accelaration vs chronological age",
#     xlab="DNAm age tumor/chronological age ", ylab="Chronological age", pch=19, col = "red") # more valuable => 241 observations
#cor.test(tumor$DNAm_age_tumor/tumor$years_to_birth, tumor$years_to_birth, method = "spearman") # Significant, weak correlation

#------Best option -----------

#Creation of a linear model to predict DNAm_age
linear_model <- lm(c(normal$DNAm_age) ~ c(normal$years_to_birth))
plot(normal$years_to_birth, normal$DNAm_age_normal, main="Chronological age vs DNAm age for normal tissue",
     xlab="Chronological age", ylab="DNAm age normal", pch=1, col = "black")
abline(linear_model, col="red")
abline(c(0,1),col="blue") 
legend("topleft", legend=c("Regression line", "Identity line", "Data"),
       col=c("red", "blue", "black"), lty= c(1,1, NA), pch = c(NA,NA, 1), cex=0.8)
summary(linear_model)

tumor$predicted_DNAm_age_normal <- c(coef(linear_model)[[2]]* tumor$years_to_birth + coef(linear_model)[[1]])
tumor$DNAm_age_acceleration <- c(tumor$DNAm_age_tumor/tumor$predicted_DNAm_age_normal)
plot(tumor$years_to_birth, tumor$DNAm_age_acceleration , main="Tumor age accelaration vs chronological age",
     xlab="Chronological age", ylab="Acceleration (DNAm age tumor/predicted DNAm age normal)", pch=1, col = "red")
cor.test(tumor$DNAm_age_acceleration, tumor$years_to_birth, method = "spearman")
abline(h=1,col="blue") 
legend("topleft", legend= c("y = 1","Data"),
       col= c("blue", "red"), lty=c(1, NA), pch = c(NA, 1), cex=0.8)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Q2.1<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Continuous variable preprocess

tumor$year_of_tobacco_smoking_onset<- as.numeric(as.character(tumor$year_of_tobacco_smoking_onset))
tumor$date_of_initial_pathologic_diagnosis<- as.numeric(as.character(tumor$date_of_initial_pathologic_diagnosis))
tumor$smoking_period <- c(tumor$date_of_initial_pathologic_diagnosis - tumor$year_of_tobacco_smoking_onset)
tumor$number_pack_years_smoked<- as.numeric(as.character(tumor$number_pack_years_smoked))
tumor$total_number_packs_smoked <- c(tumor$smoking_period * tumor$number_pack_years_smoked)
tumor$karnofsky_performance_score <- as.numeric(as.character(tumor$karnofsky_performance_score))

# Spearman correlation for continuous variable

cor.test(tumor$DNAm_age_tumor, tumor$smoking_period, method = "spearman")
cor.test(tumor$DNAm_age_acceleration, tumor$smoking_period, method = "spearman"	)


cor.test(tumor$DNAm_age_tumor, tumor$number_pack_years_smoked, method = "spearman")
cor.test(tumor$DNAm_age_acceleration, tumor$number_pack_years_smoked, method = "spearman")

cor.test(tumor$DNAm_age_tumor, tumor$total_number_packs_smoked, method = "spearman")
cor.test(tumor$total_number_packs_smoked, tumor$DNAm_age_acceleration, method = "spearman") 

cor.test(tumor$DNAm_age_tumor, tumor$karnofsky_performance_score, method = "spearman")
cor.test(tumor$DNAm_age_acceleration, tumor$karnofsky_performance_score, method = "spearman")

# Categorical variable preprocess

table(tumor$pathology_T_stage_reorganised)
tumor$pathology_T_stage_reorganised <-c(tumor$pathology_T_stage)
tumor$pathology_T_stage_reorganised <- gsub('t1.', 't1', tumor$pathology_T_stage_reorganised)
tumor$pathology_T_stage_reorganised <- gsub('t2.', 't2', tumor$pathology_T_stage_reorganised)

tumor_acceleration <- tumor[!is.na(tumor$years_to_birth),] 
table(tumor_acceleration$histological_type)
colSums(!is.na(tumor_acceleration))

table(tumor$pathology_M_stage)
tumor$pathology_M_stage_reorgarnised <-c(tumor$pathology_M_stage)
tumor$pathology_M_stage_reorgarnised <-gsub('m1.', 'm1', tumor$pathology_M_stage_reorgarnised)

# Categorical variable calculation

# Pathologic stage 
#DNAm age
myAOV = aov(DNAm_age_tumor ~ pathologic_stage, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_tumor ~ pathologic_stage, data = tumor)
kruskal.test(DNAm_age_tumor ~ pathologic_stage, data = tumor)
#Acceleration 
myAOV = aov(DNAm_age_acceleration ~ pathologic_stage, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_acceleration ~ pathologic_stage, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_acceleration ~ pathologic_stage, data = tumor)

# Pathology T stage 

#DNAm age
myAOV = aov(DNAm_age_tumor ~ pathology_T_stage_reorganised, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_tumor ~ pathology_T_stage_reorganised, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_tumor ~ pathology_T_stage_reorganised, data = tumor)
#Acceleration 
myAOV = aov(DNAm_age_acceleration ~ pathology_T_stage_reorganised, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_acceleration ~ pathology_T_stage_reorganised, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_acceleration ~ pathology_T_stage_reorganised, data = tumor)

# Pathology N stage 
#DNAm age
myAOV = aov(DNAm_age_tumor ~ pathology_N_stage , data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_tumor ~ pathology_N_stage, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_tumor ~ pathology_N_stage, data = tumor)
#Acceleration 
myAOV = aov(DNAm_age_acceleration ~ pathology_N_stage, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_acceleration ~ pathology_N_stage, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_acceleration ~ pathology_N_stage, data = tumor)

# Pathology M stage 
#DNAm age
myAOV = aov(DNAm_age_tumor ~ pathology_M_stage_reorgarnised , data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_tumor ~ pathology_M_stage_reorgarnised, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_tumor ~ pathology_M_stage_reorgarnised, data = tumor)
#Acceleration 
myAOV = aov(DNAm_age_acceleration ~ pathology_M_stage_reorgarnised, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_acceleration ~ pathology_M_stage_reorgarnised, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_acceleration ~ pathology_M_stage_reorgarnised, data = tumor)

# gender
#DNAm age
myAOV = aov(DNAm_age_tumor ~ gender , data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_tumor ~ gender, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_tumor ~ gender, data = tumor)
#Acceleration 
myAOV = aov(DNAm_age_acceleration ~ gender, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_acceleration ~ gender, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_acceleration ~ gender, data = tumor)

# Radiation Therapy
#DNAm age
myAOV = aov(DNAm_age_tumor ~ radiation_therapy , data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_tumor ~ radiation_therapy, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_tumor ~ radiation_therapy, data = tumor)
#Acceleration 
myAOV = aov(DNAm_age_acceleration ~ radiation_therapy, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_acceleration ~ radiation_therapy, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_acceleration ~ radiation_therapy, data = tumor)

# Histological type 
#DNAm age
myAOV = aov(DNAm_age_tumor ~ histological_type , data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_tumor ~ histological_type, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_tumor ~ histological_type, data = tumor)
#Acceleration 
myAOV = aov(DNAm_age_acceleration ~ histological_type, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_acceleration ~ histological_type, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_acceleration ~ histological_type, data = tumor)

# Residual tumor
#DNAm age
myAOV = aov(DNAm_age_tumor ~ residual_tumor , data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_tumor ~ residual_tumor, data = tumor) #for validation => should be > 0.05 at level alpha 5%
summary(myAOV)
#Acceleration 
myAOV = aov(DNAm_age_acceleration ~ residual_tumor, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_acceleration ~ residual_tumor, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_acceleration ~ residual_tumor, data = tumor)

# Race
#DNAm age
myAOV = aov(DNAm_age_tumor ~ race , data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_tumor ~ race, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_tumor ~ race, data = tumor)
#Acceleration 
myAOV = aov(DNAm_age_acceleration ~ race, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_acceleration ~ race, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_acceleration ~ race, data = tumor)

# Ethnicity
#DNAm age
myAOV = aov(DNAm_age_tumor ~ ethnicity , data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_tumor ~ ethnicity, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_tumor ~ ethnicity, data = tumor)
#Acceleration 
myAOV = aov(DNAm_age_acceleration ~ ethnicity, data = tumor )
shapiro.test(residuals(myAOV)) # for validation => should be > 0.05 at level alpha 5%
leveneTest(DNAm_age_acceleration ~ ethnicity, data = tumor) #for validation => should be > 0.05 at level alpha 5%
kruskal.test(DNAm_age_acceleration ~ ethnicity, data = tumor)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Q2.2<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#Survival analysis

tumor$vital_status <- as.numeric(as.character(tumor[, "vital_status"]))# transform to numeric
tumor$days_to_death <- as.numeric(as.character(tumor[,"days_to_death"])) 
tumor$days_to_last_followup <- as.numeric(as.character(tumor[,"days_to_last_followup"]))
event <- tumor$days_to_death
event[is.na(event)] <- tumor$days_to_last_followup[is.na(event)]
head(tumor)

#add kaplan-meirer curves
tumor <- within(tumor, quartile <- as.integer(cut(DNAm_age_acceleration, quantile(DNAm_age_acceleration, probs=0:4/4, na.rm=TRUE), include.lowest=TRUE)))
res_KM <- survfit(Surv(event, vital_status) ~ quartile, data = tumor)
summary(res_KM)
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(res_KM,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF", "#E74500", "#B800E7"))

# Corx analysis

res.cox_DNAm_age <- coxph(Surv(event, vital_status) ~ DNAm_age_tumor, data = tumor)
res.cox_DNAm_acc <- coxph(Surv(event, vital_status) ~ DNAm_age_acceleration, data = tumor)

summary(res.cox_DNAm_age)
summary(res.cox_DNAm_acc)