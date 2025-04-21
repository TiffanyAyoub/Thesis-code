## PERMUTATION (1000) FOR TREATMENT EFFECT ##

# Get data frame
load("/home/tayoub/MEDIAN_FOR_IL6DATA.RData")
#dt.struc <- read.csv('/hpf/largeprojects/MICe/jyeung/DSURQEE_60micron_invivo/DSURQEE_R_mapping.csv')
dt.struc <- read.csv('/projects/jyeung/DSURQEE_60micron_invivo/DSURQEE_R_mapping.csv')

colnames(MEDIAN_FOR_IL6DATA)[15:197] <- as.character(dt.struc$Structure)

#Filter
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


## SPLIT MALE AND FEMALE USING LINEAR MODEL 2
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


age_factor <- as.character(FINAL_filtered$Age)
sex_flag <- (FINAL_filtered$Sex == "Male") * 1

M_KO_flag_14 <- (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 14) * (FINAL_filtered$Genotype == "KO")
M_WT_flag_14 <- (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 14) * (FINAL_filtered$Genotype != "KO")
F_KO_flag_14 <- (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 14) * (FINAL_filtered$Genotype == "KO")
F_WT_flag_14 <- (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 14) * (FINAL_filtered$Genotype != "KO")

M_MTX_KO_flag_28 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 28) * (FINAL_filtered$Genotype == "KO")
M_MTX_WT_flag_28 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 28) * (FINAL_filtered$Genotype != "KO")
F_MTX_KO_flag_28 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 28) * (FINAL_filtered$Genotype == "KO")
F_MTX_WT_flag_28 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 28) * (FINAL_filtered$Genotype != "KO")
M_Saline_KO_flag_28 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 28) * (FINAL_filtered$Genotype == "KO")
M_Saline_WT_flag_28 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 28) * (FINAL_filtered$Genotype != "KO")
F_Saline_KO_flag_28 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 28) * (FINAL_filtered$Genotype == "KO")
F_Saline_WT_flag_28 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 28) * (FINAL_filtered$Genotype != "KO")

M_MTX_KO_flag_42 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 42) * (FINAL_filtered$Genotype == "KO")
M_MTX_WT_flag_42 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 42) * (FINAL_filtered$Genotype != "KO")
F_MTX_KO_flag_42 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 42) * (FINAL_filtered$Genotype == "KO")
F_MTX_WT_flag_42 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 42) * (FINAL_filtered$Genotype != "KO")
M_Saline_KO_flag_42 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 42) * (FINAL_filtered$Genotype == "KO")
M_Saline_WT_flag_42 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 42) * (FINAL_filtered$Genotype != "KO")
F_Saline_KO_flag_42 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 42) * (FINAL_filtered$Genotype == "KO")
F_Saline_WT_flag_42 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 42) * (FINAL_filtered$Genotype != "KO")

M_MTX_KO_flag_63 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 63) * (FINAL_filtered$Genotype == "KO")
M_MTX_WT_flag_63 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 63) * (FINAL_filtered$Genotype != "KO")
F_MTX_KO_flag_63 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 63) * (FINAL_filtered$Genotype == "KO")
F_MTX_WT_flag_63 <- (FINAL_filtered$Treatment == "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 63) * (FINAL_filtered$Genotype != "KO")
M_Saline_KO_flag_63 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 63) * (FINAL_filtered$Genotype == "KO")
M_Saline_WT_flag_63 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Male") * (FINAL_filtered$Age == 63) * (FINAL_filtered$Genotype != "KO")
F_Saline_KO_flag_63 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 63) * (FINAL_filtered$Genotype == "KO")
F_Saline_WT_flag_63 <- (FINAL_filtered$Treatment != "MTX") * (FINAL_filtered$Sex == "Female") * (FINAL_filtered$Age == 63) * (FINAL_filtered$Genotype != "KO")


# make empty data frame
struct <- dt.struc$Structure
form <- glue("`{struct}` ~ -1 + age_factor + age_factor:sex_flag + 
            age_factor:z_litter_size +
            M_KO_flag_14 + 
            F_KO_flag_14 + 
            M_MTX_KO_flag_28 + M_MTX_WT_flag_28 +
            F_MTX_KO_flag_28 + F_MTX_WT_flag_28 + 
            M_Saline_KO_flag_28 + 
            F_Saline_KO_flag_28 + 
            M_MTX_KO_flag_42 + M_MTX_WT_flag_42 + 
            F_MTX_KO_flag_42 + F_MTX_WT_flag_42 +
            M_Saline_KO_flag_42 + 
            F_Saline_KO_flag_42 + 
            M_MTX_KO_flag_63 + M_MTX_WT_flag_63 +
            F_MTX_KO_flag_63 + F_MTX_WT_flag_63 +
            M_Saline_KO_flag_63 + 
            F_Saline_KO_flag_63 + 
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

## EXTRACTING ID & GENOTYPE ## 
df_ID_treat <- distinct(FINAL_filtered, ID, Treatment, Sex, .keep_all = FALSE)
## REMOVE GENOTYPE FROM MAIN DATA FRAME ##
df_no_treat <- select(FINAL_filtered, -Treatment)

## LOOP THROUGH PERMUTATIONS BELOW 
set.seed(123)
for (p in 1:n_permutation) {
  
  df_dummy <- allres
  
  df_ID_F <- subset(df_ID_treat, Sex=="Female")
  df_ID_M <- subset(df_ID_treat, Sex=="Male")
  
  # Shuffle Genotype by ID (ensure same ID gets same Genotype across timepoints)
  df_shuffle_F <- transform(df_ID_F, Treatment = sample(df_ID_F$Treatment))
  df_shuffle_M <- transform(df_ID_M, Treatment = sample(df_ID_M$Treatment))
  df_shuffle <- rbind(df_shuffle_F, df_shuffle_M)
  df_shuffle <- df_shuffle[c("ID", "Treatment")]
  df_perm <- left_join(df_no_treat, df_shuffle, by = "ID")
  
  age_factor <- as.character(df_perm$Age)
  sex_flag <- (df_perm$Sex == "Male") * 1
  
  M_KO_flag_14 <- (df_perm$Sex == "Male") * (df_perm$Age == 14) * (df_perm$Genotype == "KO")
  M_WT_flag_14 <- (df_perm$Sex == "Male") * (df_perm$Age == 14) * (df_perm$Genotype != "KO")
  F_KO_flag_14 <- (df_perm$Sex == "Female") * (df_perm$Age == 14) * (df_perm$Genotype == "KO")
  F_WT_flag_14 <- (df_perm$Sex == "Female") * (df_perm$Age == 14) * (df_perm$Genotype != "KO")
  
  M_MTX_KO_flag_28 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 28) * (df_perm$Genotype == "KO")
  M_MTX_WT_flag_28 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 28) * (df_perm$Genotype != "KO")
  F_MTX_KO_flag_28 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 28) * (df_perm$Genotype == "KO")
  F_MTX_WT_flag_28 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 28) * (df_perm$Genotype != "KO")
  M_Saline_KO_flag_28 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 28) * (df_perm$Genotype == "KO")
  M_Saline_WT_flag_28 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 28) * (df_perm$Genotype != "KO")
  F_Saline_KO_flag_28 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 28) * (df_perm$Genotype == "KO")
  F_Saline_WT_flag_28 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 28) * (df_perm$Genotype != "KO")
  
  M_MTX_KO_flag_42 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 42) * (df_perm$Genotype == "KO")
  M_MTX_WT_flag_42 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 42) * (df_perm$Genotype != "KO")
  F_MTX_KO_flag_42 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 42) * (df_perm$Genotype == "KO")
  F_MTX_WT_flag_42 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 42) * (df_perm$Genotype != "KO")
  M_Saline_KO_flag_42 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 42) * (df_perm$Genotype == "KO")
  M_Saline_WT_flag_42 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 42) * (df_perm$Genotype != "KO")
  F_Saline_KO_flag_42 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 42) * (df_perm$Genotype == "KO")
  F_Saline_WT_flag_42 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 42) * (df_perm$Genotype != "KO")
  
  M_MTX_KO_flag_63 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 63) * (df_perm$Genotype == "KO")
  M_MTX_WT_flag_63 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 63) * (df_perm$Genotype != "KO")
  F_MTX_KO_flag_63 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 63) * (df_perm$Genotype == "KO")
  F_MTX_WT_flag_63 <- (df_perm$Treatment == "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 63) * (df_perm$Genotype != "KO")
  M_Saline_KO_flag_63 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 63) * (df_perm$Genotype == "KO")
  M_Saline_WT_flag_63 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Male") * (df_perm$Age == 63) * (df_perm$Genotype != "KO")
  F_Saline_KO_flag_63 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 63) * (df_perm$Genotype == "KO")
  F_Saline_WT_flag_63 <- (df_perm$Treatment != "MTX") * (df_perm$Sex == "Female") * (df_perm$Age == 63) * (df_perm$Genotype != "KO")
  
  ## STORE RESULTS FOR THE PERMUTATIONS
  
  ## LOOP THROUGH ALL 183 STRUCTURES
  for (j in 1:length(structure_list)) {
    cstruct <- as.character(df_dummy$Structure[j])
    
    form <- glue("`{cstruct}` ~ -1 + age_factor + age_factor:sex_flag + 
            age_factor:z_litter_size +
            M_KO_flag_14 + 
            F_KO_flag_14 + 
            M_MTX_KO_flag_28 + M_MTX_WT_flag_28 +
            F_MTX_KO_flag_28 + F_MTX_WT_flag_28 + 
            M_Saline_KO_flag_28 + 
            F_Saline_KO_flag_28 + 
            M_MTX_KO_flag_42 + M_MTX_WT_flag_42 + 
            F_MTX_KO_flag_42 + F_MTX_WT_flag_42 +
            M_Saline_KO_flag_42 + 
            F_Saline_KO_flag_42 + 
            M_MTX_KO_flag_63 + M_MTX_WT_flag_63 +
            F_MTX_KO_flag_63 + F_MTX_WT_flag_63 +
            M_Saline_KO_flag_63 + 
            F_Saline_KO_flag_63 + 
            (1|ID) + (1|Home_cage)")
    
    
    # Fit the mixed model with permuted Genotype
    clm_perm <- lmer(as.formula(form), data = df_perm)
    
    slm <- summary(clm_perm)
    df_dummy[j,allcols]<-as.vector(slm$coefficients)
    df_dummy$Permutation <- p
  }
  
  # Combine all results into one big data frame
  perm_results_df <- rbind(perm_results_df, df_dummy)
}

save(perm_results_df, file = "/projects/tayoub/PERMUTATION_RESULTS_TREATMENT_1000.RData")
write.csv(perm_results_df, "/projects/tayoub/PERMUTATION_RESULTS_TREATMENT_1000.csv", row.names = FALSE)



################################################################################
################################################################################
################################################################################

load("/projects/tayoub/PERMUTATION_RESULTS_TREATMENT_1000.RData")

perm_results_df <- perm_results_df %>%
  mutate(
    # WT FEMALE (Using age_factor28 because females are baseline)
    F_WT_MTX_28_p = (`F_MTX_WT_flag_28.Estimate` / `age_factor28.Estimate`) * 100,
    
    # WT MALE: Compute new baseline first, then normalize
    Male_WT_Saline_Estimate_28_p = `age_factor28:sex_flag.Estimate` + `age_factor28.Estimate`,
    M_WT_MTX_28_p = (`M_MTX_WT_flag_28.Estimate` / Male_WT_Saline_Estimate_28_p) * 100,
    
    # KO FEMALE: Compute new baseline first, then normalize
    F_KO_Saline_28_Estimate_p = `F_Saline_KO_flag_28.Estimate` + `age_factor28.Estimate`,
    F_KO_MTX_28_p = (`F_MTX_KO_flag_28.Estimate` / F_KO_Saline_28_Estimate_p) * 100,
    
    # KO MALE: Compute new baseline first, then normalize
    M_KO_Saline_28_Estimate_p = `M_Saline_KO_flag_28.Estimate` + Male_WT_Saline_Estimate_28_p,
    M_KO_MTX_28_p = (`M_MTX_KO_flag_28.Estimate` / M_KO_Saline_28_Estimate_p) * 100
  )

# Compute median for each permutation
perm_medians <- perm_results_df %>%
  group_by(Permutation) %>%  
  summarise(
    Median_F_WT_MTX = median(F_WT_MTX_28_p),
    Median_F_KO_MTX = median(F_KO_MTX_28_p),
    Median_M_WT_MTX = median(M_WT_MTX_28_p),
    Median_M_KO_MTX = median(M_KO_MTX_28_p)
  )

# Display result
print(perm_medians)


hist(perm_medians$Median_F_WT_MTX, breaks=20)
hist(perm_medians$Median_F_KO_MTX, breaks=20)
hist(perm_medians$Median_M_WT_MTX, breaks=20)
hist(perm_medians$Median_M_KO_MTX, breaks=20)

load("/projects/tayoub/ALLRES_2.RData")

#For WT Female (using age_fatcor28 b/c females are basline)
F_WT_MTX_28 <- (allres$F_MTX_WT_flag_28.Estimate / allres$age_factor28.Estimate) * 100

#For WT Male 
Male_WT_Saline_Estimate_28 <- allres$'age_factor28:sex_flag.Estimate' + allres$'age_factor28.Estimate'
M_WT_MTX_28 <- (allres$M_MTX_WT_flag_28.Estimate / Male_WT_Saline_Estimate_28) * 100


#For KO Female
F_KO_Saline_28_Estimate <- allres$'F_Saline_KO_flag_28.Estimate' + allres$'age_factor28.Estimate'
F_KO_MTX_28 <- (allres$F_MTX_KO_flag_28.Estimate / F_KO_Saline_28_Estimate) * 100

#For KO Male
M_KO_Saline_28_Estimate <- allres$'M_Saline_KO_flag_28.Estimate' + Male_WT_Saline_Estimate_28
M_KO_MTX_28 <- (allres$M_MTX_KO_flag_28.Estimate / M_KO_Saline_28_Estimate) * 100

observed_F_WT_MTX_median <- median(F_WT_MTX_28)
observed_F_KO_MTX_median <- median(F_KO_MTX_28)
observed_M_WT_MTX_median <- median(M_WT_MTX_28)
observed_M_KO_MTX_median <- median(M_KO_MTX_28)

print(observed_F_WT_MTX_median)
print(observed_F_KO_MTX_median)
print(observed_M_WT_MTX_median)
print(observed_M_KO_MTX_median)

# Histograms


hist(perm_medians$Median_F_WT_MTX, breaks=20, main="Histogram of Median_F_WT_MTX",
     xlab="Median_F_WT_MTX", col="white", border="black")
abline(v = observed_F_WT_MTX_median, col="red", lwd=2, lty=2)  # Add red dashed line for observed median


hist(perm_medians$Median_F_KO_MTX, breaks=20, main="Histogram of Median_F_KO_MTX",
     xlab="Median_F_KO_MTX", col="white", border="black")
abline(v = observed_F_KO_MTX_median, col="red", lwd=2, lty=2)  # Add red dashed line for observed median


hist(perm_medians$Median_M_WT_MTX, breaks=20, main="Histogram of Median_M_WT_MTX",
     xlab="Median_M_WT_MTX", col="white", border="black")
abline(v = observed_M_WT_MTX_median, col="red", lwd=2, lty=2)  # Add red dashed line for observed median


hist(perm_medians$Median_M_KO_MTX, breaks=20, main="Histogram of Median_M_KO_MTX",
     xlab="Median_M_KO_MTX", col="white", border="black")
abline(v = observed_M_KO_MTX_median, col="red", lwd=2, lty=2)  # Add red dashed line for observed median


##### TESTING #####

N <- length(perm_medians$Median_F_WT_MTX)
j <- sum(observed_F_WT_MTX_median > perm_medians$Median_F_WT_MTX)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0.078 

N <- length(perm_medians$Median_F_KO_MTX)
j <- sum(observed_F_KO_MTX_median > perm_medians$Median_F_KO_MTX)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0 

N <- length(perm_medians$Median_M_WT_MTX)
j <- sum(observed_M_WT_MTX_median > perm_medians$Median_M_WT_MTX)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0.006

N <- length(perm_medians$Median_M_KO_MTX)
j <- sum(observed_M_KO_MTX_median > perm_medians$Median_M_KO_MTX)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
mypval
#0

###########################################################################

#Comparing one histogram against another using wilcoxon signed rank test 

load("/projects/tayoub/ALLRES_2.RData")

#For WT Female (using age_fatcor28 b/c females are basline)
F_WT_MTX_28 <- (allres$F_MTX_WT_flag_28.Estimate / allres$age_factor28.Estimate) * 100

#For WT Male 
Male_WT_Saline_Estimate_28 <- allres$'age_factor28:sex_flag.Estimate' + allres$'age_factor28.Estimate'
M_WT_MTX_28 <- (allres$M_MTX_WT_flag_28.Estimate / Male_WT_Saline_Estimate_28) * 100


#For KO Female
F_KO_Saline_28_Estimate <- allres$'F_Saline_KO_flag_28.Estimate' + allres$'age_factor28.Estimate'
F_KO_MTX_28 <- (allres$F_MTX_KO_flag_28.Estimate / F_KO_Saline_28_Estimate) * 100

#For KO Male
M_KO_Saline_28_Estimate <- allres$'M_Saline_KO_flag_28.Estimate' + Male_WT_Saline_Estimate_28
M_KO_MTX_28 <- (allres$M_MTX_KO_flag_28.Estimate / M_KO_Saline_28_Estimate) * 100


######## TEST ###########

# Compare F KO vs F WT  
wilcox_female_MTX_comparison <- wilcox.test(F_WT_MTX_28, F_KO_MTX_28, 
                                           alternative = "two.sided")
print(wilcox_female_MTX_comparison)
# Significant #
# p-value = 0.0009795 #

# Compare M KO vs M WT
wilcox_male_MTX_comparison <- wilcox.test(M_WT_MTX_28, M_KO_MTX_28, 
                                            alternative = "two.sided")
print(wilcox_male_MTX_comparison)
# Significant #
# p-value = 6.635e-05 #

