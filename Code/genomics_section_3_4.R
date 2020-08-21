library(limma)
library(edgeR)
library(msigdbr)
library(EnhancedVolcano)
library(maftools)
require(data.table)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Q3 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
DNAm$sample_type <- gsub('0..', 'T', DNAm$tumor_type)
DNAm$sample_type <- gsub('1..', 'N', DNAm$sample_type)
DNAm$tumor_type <- gsub('[a-zA-Z]', "", DNAm$tumor_type)

# Load clinical data
clinic_data <- read.table('/Users/chloe/Desktop/Functional_genomics/Functional_genomics_report/LUAD.Clinical_Pick_Tier1/LUAD.clin.merged.picked.txt', header = FALSE, row.names = 1, sep = '\t')
clinic_data <- as.data.frame(t(clinic_data))

# Add patient number column 
barcode_split<- strsplit(as.character(clinic_data$`Hybridization REF`), "[-]")
clinic_data$patient_number <- sapply(barcode_split,"[",3) 

# ------------creation of tumor set and normal set ------------
merge <-merge(DNAm, clinic_data, by ="patient_number") # merge the clinical data with the DNA methylation age
merge$years_to_birth <- as.numeric(as.character(merge$years_to_birth)) # transform to numeric
# delete patient for which we don't know the years to birth
merge <- merge[!is.na(merge$years_to_birth),]

tumor <- subset(merge, sample_type=="T")
normal <- subset(merge, sample_type=="N")

# rename disctinct colum 
names(tumor)[names(tumor) == "DNAm_age"] <- "DNAm_age_tumor"
names(normal)[names(normal) == "DNAm_age"] <- "DNAm_age_normal"

#Creation of a linear model to predict DNAm_age
linear_model <- lm(c(normal$DNAm_age) ~ c(normal$years_to_birth))
tumor$predicted_DNAm_age_normal <- c(coef(linear_model)[[2]]* tumor$years_to_birth + coef(linear_model)[[1]])
normal$predicted_DNAm_age_normal <- c(coef(linear_model)[[2]]* normal$years_to_birth + coef(linear_model)[[1]])
tumor$DNAm_age_acceleration <- c(tumor$DNAm_age_tumor/tumor$predicted_DNAm_age_normal)
normal$DNAm_age_acceleration <- c(normal$DNAm_age_normal/normal$predicted_DNAm_age_normal)

#Load rna-seq
rna <- read.table("/Users/chloe/Desktop/Functional_genomics/Functional_genomics_report/gdac.broadinstitute.org_LUAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/LUAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt", row.names=1, sep='\t')
rna <- rna[,(rna[2,] =="raw_count")] # deleted column scaled_estimate & transcript_id
#column name
names(rna) <- as.matrix(rna[1, ])
rna <- rna[-1, ]
rna <- rna[-1,]
rna[] <- lapply(rna, function(x) type.convert(as.numeric(as.character(x))))
#gene id clean
gene_list <-strsplit(as.character(rownames(rna) ), "[|]")
rownames(rna) <- c(sapply(gene_list, "[[", 2))

#get participant and sample type
dataframe <- as.data.frame(do.call(rbind,strsplit(colnames(rna), ".", fixed=TRUE)))
colnames(dataframe) <- c("Barcode")
head(dataframe)
barcode_list <-strsplit(as.character(dataframe$Barcode), "[-]")
dataframe <- cbind(dataframe, patient_number = c(sapply(barcode_list, "[[", 3)))
dataframe <- cbind(dataframe, tumor_type = c(sapply(barcode_list, "[[", 4)))
dataframe <- cbind(dataframe, sample_type = c(sapply(barcode_list, "[[", 4)))
dataframe$tumor_type <- gsub('[a-zA-Z]', "", dataframe$tumor_type)
levels(dataframe$sample_type)<- c("T", "T","T", "N", "N")
dataframe$Barcode <-gsub("-", ".", dataframe$Barcode)
#----------------

df_participant_tumor <- dataframe[(dataframe$sample_type != "N"),] #get all tumor participant rnaseq 
df_participant_normal <- dataframe[(dataframe$sample_type == "N"),] #get all normal participant for which we got the rnaseq

# Merging by patient_number AND tumor type (since different tumor type coexist)
tumor_initial <-tumor
normal_initial <- normal
tumor <- merge(tumor, df_participant_tumor, by= c("patient_number","tumor_type", "sample_type"), all= FALSE)
normal <- merge(normal, df_participant_normal, by= c("patient_number","tumor_type", "sample_type"), all= FALSE) 

#Creation of two subset of mRNA-seq data, one for healthy sample, one for tumor samples.
colnames(rna) <-gsub("-", ".", colnames(rna))
rna_tumor <- rna[,(dataframe$sample_type != "N")]
rna_normal <- rna[,(dataframe$sample_type == "N")]

rna_tumor <- rna_tumor[tumor$Barcode.y ]
rna_normal <-rna_normal[normal$Barcode.y]

#-------------------------------------------------------------------------------

# Creation of a design matrix 

design_tumor <- model.matrix(~DNAm_age_acceleration + years_to_birth, data = tumor)
design_normal <- model.matrix( ~DNAm_age_acceleration +years_to_birth, data = normal)

# DGE object
dge_tumor <- DGEList(counts= rna_tumor)
dge_normal<- DGEList(counts= rna_normal)

# Filterout low count genes

keep<-filterByExpr(dge_tumor,design = design_tumor)
dge_tumor<-dge_tumor[keep,]

keep<-filterByExpr(dge_normal,design = design_normal)
dge_normal<-dge_normal[keep,]

# normalization
dge_tumor <- calcNormFactors(dge_tumor)
dge_normal <- calcNormFactors(dge_normal)

# Voom transform data 
v_normal <- voom(dge_normal, design_normal, plot=TRUE)
v_tumor <- voom(dge_tumor, design_tumor, plot=TRUE)

# Loading the C2:CP dataset

msigdbr = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
msigdbr_cleaned <-(tapply(msigdbr$entrez_gene, msigdbr$gs_id, as.character,simplify =FALSE))  #list of character vectors, each vector containing the gene identifiers for a set of genes.

index_tumor <- ids2indices(msigdbr_cleaned, row.names(v_tumor), remove.empty=TRUE)
index_normal <-ids2indices(msigdbr_cleaned, row.names(v_normal), remove.empty=TRUE)

# Camera function

camera(v_tumor, index_tumor, design = design_tumor, contrast = "DNAm_age_acceleration" )
camera(v_normal, index_normal, design = design_normal, contrast = "DNAm_age_acceleration" )

#Find DE genes tumor
fit <- lmFit(v_tumor, design_tumor)
fit <- contrasts.fit(fit, coefficients = 2)
fit <- eBayes(fit)
plotSA(fit, main="Final model: Mean-variance trend")
table <- topTable(fit, sort.by = "p", n =15)
length(which(table$adj.P.Val < 0.05))  #4503 DE genes
volcanoplot(fit, highlight = 8,names = rownames(rna_tumor) )

summary(decideTests(fit))
colnames(table)
EnhancedVolcano(table,
                lab = rownames(table),
                x = 'logFC',
                y = 'adj.P.Val',
                xlim = c(-5, 5))

#Find DE genes normal
fit <- lmFit(v_normal, design_normal)
fit <- contrasts.fit(fit, coefficients = 2)
fit <- eBayes(fit)
plotSA(fit, main="Final model: Mean-variance trend")
topTable(fit, sort.by = "p")

#-----------------------------BONUS1---------------------------------
library(plyr)
# Load mutation files

maf_files <- list.files("~/Desktop/Functional_genomics/Functional_genomics_report/LUAD_mutation_packagers_call", pattern="\\.maf.txt$", full.names=TRUE)
#maf_merged <- merge_mafs(maf_files)
MAF <- ldply(maf_files, read.delim)
number_somatic_tumor <- setNames(data.frame(table(MAF$Tumor_Sample_Barcode)), c("Barcode", "n_somatic_mutation"))

# Split barcode into sample type, tumor type and patient number
barcode_split<- strsplit(as.character(number_somatic_tumor$Barcode), "[-]")
number_somatic_tumor$patient_number  <- sapply(barcode_split,"[",3) # Add the patient number column 
number_somatic_tumor$tumor_type  <- sapply(barcode_split,"[",4) # T for tumor, N for normal
number_somatic_tumor$sample_type <- gsub('0..', 'T', number_somatic_tumor$tumor_type)
number_somatic_tumor$sample_type <- gsub('1..', 'N', number_somatic_tumor$sample_type)
number_somatic_tumor$tumor_type <- gsub('[a-zA-Z]', "", number_somatic_tumor$tumor_type)

# Merge clinical data df (with DNAm_acceleration) with df containing the number of somatic tumor by patient and tumor type.
tumor_mutations <- merge(tumor_initial, number_somatic_tumor, by= c("patient_number","tumor_type"), all= FALSE) #addd column with number of somatic mutations

# Spearman correlation
cor.test(tumor_mutations$n_somatic_mutation , tumor_mutations$DNAm_age_acceleration, method = "spearman")

#-----------------------------BONUS2-----------------------------------

#Is there a relation between DNAm age and the number of DNA copy number breakpoints?

CNV <- read.table("/Users/chloe/Desktop/Functional_genomics/Functional_genomics_report/gdac.broadinstitute.org_LUAD.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

# clean
barcode_split<- strsplit(as.character(CNV$Sample), "[-]")
CNV$patient_number  <- sapply(barcode_split,"[",3) # Add the patient number column 
CNV$tumor_type  <- sapply(barcode_split,"[",4) # T for tumor, N for normal
CNV$sample_type <- gsub('0..', 'T', CNV$tumor_type)
CNV$sample_type <- gsub('1..', 'N', CNV$sample_type)
CNV$tumor_type <- gsub('[a-zA-Z]', "", CNV$tumor_type)

CNV_tumor <- subset(CNV, sample_type=="T")
CNV_normal <- subset(CNV, sample_type=="N")

# Add label 0 for neutral (-0.3< neutral <0.3), 1 for amplfication and -1 for deletion 
CNV_tumor$label <- ifelse(CNV_tumor$Segment_Mean >-0.3 & CNV_tumor$Segment_Mean <0.3, 0,NA)
CNV_tumor$label <- ifelse(CNV_tumor$Segment_Mean >=0.3, 1,CNV_tumor$label)
CNV_tumor$label <- ifelse(CNV_tumor$Segment_Mean <= -0.3, -1,CNV_tumor$label)

# Number of amplification for a patient + number of deletion of a patient

CNV_tumor_up <- subset(CNV_tumor, label == 1)
CNV_tumor_down <- subset(CNV_tumor, label == -1)

CNV_tumor_patient_up <-setNames(data.frame(table(CNV_tumor_up$Sample)), c("Barcode", "number_up"))
CNV_tumor_patient_down <- setNames(data.frame(table(CNV_tumor_down$Sample)), c("Barcode", "number_down"))

CNV_tumor_patient <- merge(CNV_tumor_patient_up, CNV_tumor_patient_down, by= "Barcode", all = TRUE)
barcode_split<- strsplit(as.character(CNV_tumor_patient$Barcode), "[-]")
CNV_tumor_patient$patient_number  <- sapply(barcode_split,"[",3) # Add the patient number column 

# Merge with DNAm age tumor by patient
CNV_tumor_patient_merged <- merge(tumor, CNV_tumor_patient, by= c("patient_number"), all= FALSE)

# Calculate correlation between DNAm_acc and number_up (number of amplification)
cor.test(CNV_tumor_patient_merged$number_up , CNV_tumor_patient_merged$DNAm_age_acceleration, method = "spearman")
plot( CNV_tumor_patient_merged$number_up, CNV_tumor_patient_merged$DNAm_age_acceleration, main="DNAm age acceleration vs number of segments with amplification for cancer tissues",
     xlab="Number of segments with amplification", ylab="DNAm age acceleration", pch=1, col = "black")


# Calculate correlation between DNAm_acc and number_down (number of delection)
cor.test(CNV_tumor_patient_merged$number_down , CNV_tumor_patient_merged$DNAm_age_acceleration, method = "spearman")
plot(CNV_tumor_patient_merged$number_down, CNV_tumor_patient_merged$DNAm_age_acceleration, main="DNAm age acceleration vs number of segments with deletion for cancer tissues",
     xlab="Number of segments with deletion", ylab="DNAm age acceleration", pch=1, col = "black")


# Calculate correlation for between DNAm_acc and number_down and number_up (number of delection/amplification)
cor.test(CNV_tumor_patient_merged$number_down + CNV_tumor_patient_merged$number_up, CNV_tumor_patient_merged$DNAm_age_acceleration, method = "spearman")
plot(CNV_tumor_patient_merged$number_down + CNV_tumor_patient_merged$number_up, CNV_tumor_patient_merged$DNAm_age_acceleration, main="DNAm age acceleration vs number of segments with deletions or amplifications for cancer tissues",
     xlab="Number of segments with deletions or amplifications", ylab="DNAm age acceleration", pch=1, col = "black")

