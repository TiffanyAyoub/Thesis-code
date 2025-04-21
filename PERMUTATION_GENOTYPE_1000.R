#### PERMUTATIONS FOR IL6KO MICE ####

## LIBRARY ##
library(lme4)
library(dplyr)
library(glue)

# LOAD DATA
load("/home/tayoub/MEDIAN_FOR_IL6DATA.RData")
dt.struc <- read.csv('/projects/jyeung/DSURQEE_60micron_invivo/DSURQEE_R_mapping.csv')
colnames(MEDIAN_FOR_IL6DATA)[15:197] <- as.character(dt.struc$Structure)

df_FINAL_fil <- MEDIAN_FOR_IL6DATA %>%
  filter(!grepl("G11", ID))



FINAL_filtered <- df_FINAL_fil %>%
  filter(!(ID == "G6_3" & Age == 28),  
         !(ID == "G13_8" & Age == 63),  
         !(ID == "G25_1" & Age == 14),
         !(ID == "G10_1" & Age == 42))


## ENSURE GROUP IS FACTOR ##
FINAL_filtered$Genotype <- as.factor(FINAL_filtered$Genotype)

## EXTRACT 183 BRAIN STRUCTURES ## 
structure_list <- dt.struc$Structure

#### LINEAR MODEL 1 ####
#List the Litter size for each Home cage
FINAL_filtered$Litter_size[ !duplicated(FINAL_filtered$Home_cage) ]
#Check how many unique litters
length(FINAL_filtered$Litter_size[ !duplicated(FINAL_filtered$Home_cage) ])
#Mean
mean(FINAL_filtered$Litter_size[ !duplicated(FINAL_filtered$Home_cage) ])
#SD
sd(FINAL_filtered$Litter_size[!duplicated(FINAL_filtered$Home_cage)])

mean_litter_size <- mean(FINAL_filtered$Litter_size[!duplicated(FINAL_filtered$Home_cage)])
sd_litter_size <- sd(FINAL_filtered$Litter_size[!duplicated(FINAL_filtered$Home_cage)])

# Create the z-score for Litter_size
FINAL_filtered <- FINAL_filtered %>%
  mutate(z_litter_size = (Litter_size - mean_litter_size) / sd_litter_size)


#### LINEAR MODEL ####
age_factor <- as.character(FINAL_filtered$Age)
sex_flag <- (FINAL_filtered$Sex == "Male") * 1

# Treatment flags for each time point
M_treatment_flag_28 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 28)
F_treatment_flag_28 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 28)
M_treatment_flag_42 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 42)
F_treatment_flag_42 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 42)
M_treatment_flag_63 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 63)
F_treatment_flag_63 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 63)

# Genotype flags for each time point
M_genotype_flag_14 <- (FINAL_filtered$Genotype == "KO") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 14)
F_genotype_flag_14 <- (FINAL_filtered$Genotype == "KO") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 14)
M_genotype_flag_28 <- (FINAL_filtered$Genotype == "KO") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 28)
F_genotype_flag_28 <- (FINAL_filtered$Genotype == "KO") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 28)
M_genotype_flag_42 <- (FINAL_filtered$Genotype == "KO") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 42)
F_genotype_flag_42 <- (FINAL_filtered$Genotype == "KO") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 42)
M_genotype_flag_63 <- (FINAL_filtered$Genotype == "KO") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 63)
F_genotype_flag_63 <- (FINAL_filtered$Genotype == "KO") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 63)

# make empty data frame
struct <- dt.struc$Structure
form <- glue("`{struct}` ~ -1 + age_factor + age_factor:sex_flag + 
            age_factor:z_litter_size +
            M_treatment_flag_28 + F_treatment_flag_28 +
            M_treatment_flag_42 + F_treatment_flag_42 + 
            M_treatment_flag_63 + F_treatment_flag_63 + 
            M_genotype_flag_14 + F_genotype_flag_14 +
            M_genotype_flag_28 + F_genotype_flag_28 +
            M_genotype_flag_42 + F_genotype_flag_42 + 
            M_genotype_flag_63 + F_genotype_flag_63 +
            M_treatment_flag_28:M_genotype_flag_28 +
            F_treatment_flag_28:F_genotype_flag_28 +
            M_treatment_flag_42:M_genotype_flag_42 +
            F_treatment_flag_42:F_genotype_flag_42 +
            M_treatment_flag_63:M_genotype_flag_63 +
            F_treatment_flag_63:F_genotype_flag_63 +
            (1|ID) + (1|Home_cage)")
clm <- lmer(as.formula(form), data = FINAL_filtered)
structure_list <- dt.struc$Structure


slm <- summary(clm)
allcols <- as.vector(outer(rownames(slm$coefficients), colnames(slm$coefficients), paste, sep="."))
allres <- data.frame(Structure=structure_list)
allres[,allcols] <- NA

allres$Permutation <- NA


## NUMBER OF PERMUTATIONS ##
n_permutation <- 1000
perm_results_df <- data.frame()  # Empty rn

## EXTRACTING ID, GENOTYPE & SEX ## 
# Ensures each mouse ID is assigned a single genotype and sex #
df_ID_gen <- distinct(FINAL_filtered, ID, Genotype, Sex, .keep_all = FALSE)

## REMOVE GENOTYPE FROM MAIN DATA FRAME ##
# Creates a version of the dataset without genotype b/c it will be replaced by the permuted genotype below #
df_no_gen <- select(FINAL_filtered, -Genotype)

## LOOP THROUGH PERMUTATIONS BELOW 
set.seed(123)
for (p in 1:n_permutation) {
  
  df_dummy <- allres
  
  # Splitting the subsetted data by sex #
  #This will preserve the sex-specific distribution of genotypes when doing the permutation#
  df_ID_F <- subset(df_ID_gen, Sex=="Female")
  df_ID_M <- subset(df_ID_gen, Sex=="Male")
  
  # Shuffle Genotype by ID#
  #This ensures that the genotype assignment remains balanced within each sex# 
  df_shuffle_F <- transform(df_ID_F, Genotype = sample(df_ID_F$Genotype))
  df_shuffle_M <- transform(df_ID_M, Genotype = sample(df_ID_M$Genotype))
  #Combining them back together 
  df_shuffle <- rbind(df_shuffle_F, df_shuffle_M)
  #Extract only the ID and genotype column#
  #The df_shuffle will be merged back into the data frame so we don't need sex since its in the og data frame#
  df_shuffle <- df_shuffle[c("ID", "Genotype")]
  #Merge the shuffled genotype column back into the original data frame#
  #df_perm will now be identical to original data frame except the genotype column is different (permuated)#
  #left_join() ensures all occurrences of an ID recieve the same shuffled genotype# 
  df_perm <- left_join(df_no_gen, df_shuffle, by = "ID")
  
  # Recompute ALL Treatment & Genotype flags after shuffling 
  df_perm$age_factor <- as.character(df_perm$Age)
  df_perm$sex_flag <- (df_perm$Sex == "Male") * 1
  
  df_perm$M_treatment_flag_28 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 28)
  df_perm$F_treatment_flag_28 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 28)
  df_perm$M_treatment_flag_42 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 42)
  df_perm$F_treatment_flag_42 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 42)
  df_perm$M_treatment_flag_63 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 63)
  df_perm$F_treatment_flag_63 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 63)
  
  df_perm$M_genotype_flag_14 <- (df_perm$Genotype == "KO") * (df_perm$Sex == "Male") * (df_perm$Age == 14)
  df_perm$F_genotype_flag_14 <- (df_perm$Genotype == "KO") * (df_perm$Sex == "Female") * (df_perm$Age == 14)
  df_perm$M_genotype_flag_28 <- (df_perm$Genotype == "KO") * (df_perm$Sex == "Male") * (df_perm$Age == 28)
  df_perm$F_genotype_flag_28 <- (df_perm$Genotype == "KO") * (df_perm$Sex == "Female") * (df_perm$Age == 28)
  df_perm$M_genotype_flag_42 <- (df_perm$Genotype == "KO") * (df_perm$Sex == "Male") * (df_perm$Age == 42)
  df_perm$F_genotype_flag_42 <- (df_perm$Genotype == "KO") * (df_perm$Sex == "Female") * (df_perm$Age == 42)
  df_perm$M_genotype_flag_63 <- (df_perm$Genotype == "KO") * (df_perm$Sex == "Male") * (df_perm$Age == 63)
  df_perm$F_genotype_flag_63 <- (df_perm$Genotype == "KO") * (df_perm$Sex == "Female") * (df_perm$Age == 63)
  
  ## STORE RESULTS FOR THE PERMUTATIONS
  
  ## LOOP THROUGH ALL 183 STRUCTURES
  for (j in 1:length(structure_list)) {
    cstruct <- as.character(df_dummy$Structure[j])
    
    form <- glue("`{cstruct}` ~ -1 + factor(Age) + age_factor:sex_flag + 
                  age_factor:z_litter_size +
                  M_treatment_flag_28 + F_treatment_flag_28 +
                  M_treatment_flag_42 + F_treatment_flag_42 + 
                  M_treatment_flag_63 + F_treatment_flag_63 + 
                  M_genotype_flag_14 + F_genotype_flag_14 +
                  M_genotype_flag_28 + F_genotype_flag_28 +
                  M_genotype_flag_42 + F_genotype_flag_42 + 
                  M_genotype_flag_63 + F_genotype_flag_63 +
                  M_treatment_flag_28:M_genotype_flag_28 +
                  F_treatment_flag_28:F_genotype_flag_28 +
                  M_treatment_flag_42:M_genotype_flag_42 +
                  F_treatment_flag_42:F_genotype_flag_42 +
                  M_treatment_flag_63:M_genotype_flag_63 +
                  F_treatment_flag_63:F_genotype_flag_63 +
                  (1|ID) + (1|Home_cage)")
    
    # Fit the mixed model with permuted Genotype
    clm_perm <- lmer(as.formula(form), data = df_perm)
    
    slm <- summary(clm_perm)
    #Extract and store model coefficients in df_dummy for the current brain structure (j)#
    #Keeps track of the estimated effects for each brain structure#
    df_dummy[j,allcols]<-as.vector(slm$coefficients)
    #Assigns the current permutation number to all rows in df_dummy
    #This helps track which permutation the results belong to when combining multiple permutations below#
    df_dummy$Permutation <- p
  }
  
  
  # Combine all results into one big data frame
  perm_results_df <- rbind(perm_results_df, df_dummy)
}

save(perm_results_df, file = "/projects/tayoub/PERMUTATION_RESULTS_1000.RData")
write.csv(perm_results_df, "/projects/tayoub/PERMUTATION_RESULTS_1000.csv", row.names = FALSE)



##############################################################################################
##############################################################################################
##############################################################################################

load("/projects/tayoub/PERMUTATION_RESULTS_1000.RData")

library(dplyr)
library(ggplot2)


# Normalize estimates for each age and sex
perm_results_df <- perm_results_df %>%
  mutate(
    # FEMALE: Normalize interaction by female baseline at each age
    f_estimates_14_p = (`F_genotype_flag_14.Estimate` / `age_factor14.Estimate`) * 100,
    f_estimates_28_p = (`F_genotype_flag_28.Estimate` / `age_factor28.Estimate`) * 100,
    f_estimates_42_p = (`F_genotype_flag_42.Estimate` / `age_factor42.Estimate`) * 100,
    f_estimates_63_p = (`F_genotype_flag_63.Estimate` / `age_factor63.Estimate`) * 100,
    
    # MALE: Compute new baseline first, then normalize
    Male_WT_Saline_Estimate_14_p = `age_factor14:sex_flag.Estimate` + `age_factor14.Estimate`,
    Male_WT_Saline_Estimate_28_p = `age_factor28:sex_flag.Estimate` + `age_factor28.Estimate`,
    Male_WT_Saline_Estimate_42_p = `age_factor42:sex_flag.Estimate` + `age_factor42.Estimate`,
    Male_WT_Saline_Estimate_63_p = `age_factor63:sex_flag.Estimate` + `age_factor63.Estimate`,
    
    m_estimates_14_p = (`M_genotype_flag_14.Estimate` / Male_WT_Saline_Estimate_14_p) * 100,
    m_estimates_28_p = (`M_genotype_flag_28.Estimate` / Male_WT_Saline_Estimate_28_p) * 100,
    m_estimates_42_p = (`M_genotype_flag_42.Estimate` / Male_WT_Saline_Estimate_42_p) * 100,
    m_estimates_63_p = (`M_genotype_flag_63.Estimate` / Male_WT_Saline_Estimate_63_p) * 100
  )

# Compute median for each permutation
perm_medians <- perm_results_df %>%
  group_by(Permutation) %>%
  summarise(
    Median_M_Gen_14_p = median(m_estimates_14_p, na.rm = TRUE),
    Median_M_Gen_28_p = median(m_estimates_28_p, na.rm = TRUE),
    Median_M_Gen_42_p = median(m_estimates_42_p, na.rm = TRUE),
    Median_M_Gen_63_p = median(m_estimates_63_p, na.rm = TRUE),
    
    Median_F_Gen_14_p = median(f_estimates_14_p, na.rm = TRUE),
    Median_F_Gen_28_p = median(f_estimates_28_p, na.rm = TRUE),
    Median_F_Gen_42_p = median(f_estimates_42_p, na.rm = TRUE),
    Median_F_Gen_63_p = median(f_estimates_63_p, na.rm = TRUE)
  )

# Print the median results
print(perm_medians)


load("/projects/tayoub/ALLRES.RData")

# baseline estimates for females 
age_factor_14 <- allres$`age_factor14.Estimate`
age_factor_28 <- allres$`age_factor28.Estimate`
age_factor_42 <- allres$`age_factor42.Estimate`
age_factor_63 <- allres$`age_factor63.Estimate`

# Male baselines (WT+Saline) for each age
Male_WT_Saline_Estimate_14 <- age_factor_14 + allres$`age_factor14:sex_flag.Estimate`
Male_WT_Saline_Estimate_28 <- age_factor_28 + allres$`age_factor28:sex_flag.Estimate`
Male_WT_Saline_Estimate_42 <- age_factor_42 + allres$`age_factor42:sex_flag.Estimate`
Male_WT_Saline_Estimate_63 <- age_factor_63 + allres$`age_factor63:sex_flag.Estimate`

# Normalize estimates for each term

#female first
f_estimates_14 <- (allres$`F_genotype_flag_14.Estimate` / age_factor_14) * 100
f_estimates_28 <- (allres$`F_genotype_flag_28.Estimate` / age_factor_28) * 100
f_estimates_42 <- (allres$`F_genotype_flag_42.Estimate` / age_factor_42) * 100
f_estimates_63 <- (allres$`F_genotype_flag_63.Estimate` / age_factor_63) * 100

#now male
m_estimates_14 <- (allres$`M_genotype_flag_14.Estimate` / Male_WT_Saline_Estimate_14) * 100
m_estimates_28 <- (allres$`M_genotype_flag_28.Estimate` / Male_WT_Saline_Estimate_28) * 100
m_estimates_42 <- (allres$`M_genotype_flag_42.Estimate` / Male_WT_Saline_Estimate_42) * 100
m_estimates_63 <- (allres$`M_genotype_flag_63.Estimate` / Male_WT_Saline_Estimate_63) * 100


observed_F_med_14 <- median(f_estimates_14)
observed_F_med_28 <- median(f_estimates_28)
observed_F_med_42 <- median(f_estimates_42)
observed_F_med_63 <- median(f_estimates_63)

observed_M_med_14 <- median(m_estimates_14)
observed_M_med_28 <- median(m_estimates_28)
observed_M_med_42 <- median(m_estimates_42)
observed_M_med_63 <- median(m_estimates_63)

## HISTOGRAMS ##

## FEMALES FIRST ##

hist(perm_medians$Median_F_Gen_14_p, breaks=20, main="Histogram of Median_F_Gen_14_p",
     xlab="Median_F_Gen_14_p", col="white", border="black")
abline(v = observed_F_med_14, col="red", lwd=2, lty=2)  # Add red dashed line for observed median

hist(perm_medians$Median_F_Gen_28_p, breaks=20, main="Histogram of Median_F_Gen_28_p",
     xlab="Median_F_Gen_28_p", col="white", border="black")
abline(v = observed_F_med_28, col="red", lwd=2, lty=2)  # Add red dashed line for observed median

hist(perm_medians$Median_F_Gen_42_p, breaks=20, main="Histogram of Median_F_Gen_42_p",
     xlab="Median_F_Gen_42_p", col="white", border="black")
abline(v = observed_F_med_42, col="red", lwd=2, lty=2)  # Add red dashed line for observed median

hist(perm_medians$Median_F_Gen_63_p, breaks=20, main="Histogram of Median_F_Gen_63_p",
     xlab="Median_F_Gen_63_p", col="white", border="black")
abline(v = observed_F_med_63, col="red", lwd=2, lty=2)  # Add red dashed line for observed median


## MALES ## 

hist(perm_medians$Median_M_Gen_14_p, breaks=20, main="Histogram of Median_M_Gen_14_p",
     xlab="Median_M_Gen_14_p", col="white", border="black")
abline(v = observed_M_med_14, col="red", lwd=2, lty=2)  # Add red dashed line for observed median

hist(perm_medians$Median_M_Gen_28_p, breaks=20, main="Histogram of Median_M_Gen_28_p",
     xlab="Median_M_Gen_28_p", col="white", border="black")
abline(v = observed_M_med_28, col="red", lwd=2, lty=2)  # Add red dashed line for observed median

hist(perm_medians$Median_M_Gen_42_p, breaks=20, main="Histogram of Median_M_Gen_42_p",
     xlab="Median_M_Gen_42_p", col="white", border="black")
abline(v = observed_M_med_42, col="red", lwd=2, lty=2)  # Add red dashed line for observed median

hist(perm_medians$Median_M_Gen_63_p, breaks=20, main="Histogram of Median_M_Gen_63_p",
     xlab="Median_M_Gen_63_p", col="white", border="black")
abline(v = observed_M_med_63, col="red", lwd=2, lty=2)  # Add red dashed line for observed median


#### TESTING ####

#two-tailed empirical p-value

## FEMALES ##

# P14 #
N <- length(perm_medians$Median_F_Gen_14_p)
j <- sum(observed_F_med_14 > perm_medians$Median_F_Gen_14_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0.772 

# P28 #
N <- length(perm_medians$Median_F_Gen_28_p)
j <- sum(observed_F_med_28 > perm_medians$Median_F_Gen_28_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0.264

# P42 #
N <- length(perm_medians$Median_F_Gen_42_p)
j <- sum(observed_F_med_42 > perm_medians$Median_F_Gen_42_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0.106

# P63 #
N <- length(perm_medians$Median_F_Gen_63_p)
j <- sum(observed_F_med_63 > perm_medians$Median_F_Gen_63_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0.066

## MALES ##

# P14 #
N <- length(perm_medians$Median_M_Gen_14_p)
j <- sum(observed_M_med_14 > perm_medians$Median_M_Gen_14_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0.01 

# P28 #
N <- length(perm_medians$Median_M_Gen_28_p)
j <- sum(observed_M_med_28 > perm_medians$Median_M_Gen_28_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0.116

# P42 #
N <- length(perm_medians$Median_M_Gen_42_p)
j <- sum(observed_M_med_42 > perm_medians$Median_M_Gen_42_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0.182

# P63 #
N <- length(perm_medians$Median_M_Gen_63_p)
j <- sum(observed_M_med_63 > perm_medians$Median_M_Gen_63_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0.4
