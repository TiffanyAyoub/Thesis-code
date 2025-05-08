###### THESIS R CODE ######
###### TIFFANY AYOUB ######

#########################################################
## LOAD PACKAGES                                       ##
#########################################################

library(tidyverse)
library(glue)
library(lme4)
library(lmerTest)
library(ggplot2)
library(gridExtra)


#########################################################
## LOAD DATA                                           ##
#########################################################

## LOAD: "FINAL_DATA_FRAME.RData ## THIS IS THE CSV FILE WITH BRAIN VOLUMES FROM THE STUDY ##
load("FINAL_DATA_FRAME.RData")
df<-FINAL_DATA_FRAME
## LOAD: DSURQEE_R_mapping.csv ## THIS IS THE CSV FILE WITH LABELS FOR THE BRAIN ATLAS WHICH CONTAINS 183 BRAIN STRUCTURES ##
dt.struc<-read.csv2("DSURQEE_R_mapping.csv",sep=",")
## ENSURE MATCHED STRUCTURE NAMES IN ATLAS TO STRUCTURE NAMES IN DATA ##
dt.struc$Structure <- colnames(df)[which(colnames(df) == "amygdala"):(ncol(df))]


#########################################################
## SUMMARY VOLUME CALCULATION                          ##
#########################################################

## FILTER ATLAS TO GET ONLY THE STRUCTURES LISTED UNDER WHITE MATTER TISSUE TYPE ##
## SUM ALL WHITE MATTER STRUCTURES ##
WM_list <- subset(dt.struc, tissue.type == "WM")[,1] 
WM_list <- c(WM_list, "anterior lobule white matter") # ADDING STRUCTURE MANUALLY THAT IS WM #
# WHITE MATTER LIST CONSISTS OF 40 STRUCTURES #
df$WM <- rowSums(df[,as.character(WM_list)], na.rm = TRUE)

## FILTER ATLAS TO GET ONLY THE STRUCTURES LISTED UNDER GREY MATTER TISSUE TYPE ##
## SUM ALL GREY MATTER STRUCTURES ##
GM_list <- subset(dt.struc, tissue.type == "GM")[,1]
GM_list <- c(GM_list, "midbrain", "medulla", "lateral septum", "medial septum") # ADDING STRUCTURES MANUALLY -- THESE ARE MIXED TISSUE MADE UP MOSTLY OF GM #
# GREY MATTER LIST CONSISTS OF 138 STRUCTURES #
df$GM <- rowSums(df[,GM_list], na.rm = TRUE)

## FILTER ATLAS TO GET ONLY THE STRUCTURES LISTED UNDER HIPPOCAMPUS TISSUE TYPE ##
## SUM ALL HIPPOCAMPUS STRUCTURES ##
HP_list <- subset(dt.struc, hierarchy == "Hippocampal region")[,1]
# HIPPOCAMPUS HAS 15 STRUCTUES #
df$HP <- rowSums(df[,HP_list], na.rm = TRUE)

## FILTER ATLAS TO GET ONLY THE STRUCTURES LISTED UNDER PARIETAL REGION TISSUE TYPE ##
## SUM ALL PR STRUCTURES ##
PR_list <- subset(dt.struc, hierarchy == "Parietal region")[,1]
# PR HAS 13 STRUCTUES #
df$PR <- rowSums(df[,PR_list], na.rm = TRUE)


#########################################################
## LINEAR MODEL SETUP                                  ##
#########################################################

## FIRST SETUP CONVENIENT FLAGS FOR LINEAR MODEL ##
df$age_factor <- factor(as.character(df$Age),levels=c("14","28","42","63"))
df$sex_flag <- (df$Sex == "Male") * 1

## TREATMENT FLAGS FOR EACH TIME POINT ##
df$M_treatment_flag_28 <- (df$Treatment == "MTX") * (df$Sex == "Male") * (df$Age == 28)
df$F_treatment_flag_28 <- (df$Treatment == "MTX") * (df$Sex == "Female") * (df$Age == 28)
df$M_treatment_flag_42 <- (df$Treatment == "MTX") * (df$Sex == "Male") * (df$Age == 42)
df$F_treatment_flag_42 <- (df$Treatment == "MTX") * (df$Sex == "Female") * (df$Age == 42)
df$M_treatment_flag_63 <- (df$Treatment == "MTX") * (df$Sex == "Male") * (df$Age == 63)
df$F_treatment_flag_63 <- (df$Treatment == "MTX") * (df$Sex == "Female") * (df$Age == 63)

## GENOTYPE FLAGS FOR EACH TIME POINT ##
df$M_genotype_flag_14 <- (df$Genotype == "KO") * (df$Sex == "Male") * (df$Age == 14)
df$F_genotype_flag_14 <- (df$Genotype == "KO") * (df$Sex == "Female") * (df$Age == 14)
df$M_genotype_flag_28 <- (df$Genotype == "KO") * (df$Sex == "Male") * (df$Age == 28)
df$F_genotype_flag_28 <- (df$Genotype == "KO") * (df$Sex == "Female") * (df$Age == 28)
df$M_genotype_flag_42 <- (df$Genotype == "KO") * (df$Sex == "Male") * (df$Age == 42)
df$F_genotype_flag_42 <- (df$Genotype == "KO") * (df$Sex == "Female") * (df$Age == 42)
df$M_genotype_flag_63 <- (df$Genotype == "KO") * (df$Sex == "Male") * (df$Age == 63)
df$F_genotype_flag_63 <- (df$Genotype == "KO") * (df$Sex == "Female") * (df$Age == 63)

## CALCULATE MEAN AND SD OF LITTER SIZE (ONLY ONCE PER HOME CAGE TO AVOID COUNTING THE SAME CAGE MULTIPLE TIMES)
mean_litter_size <- mean(df$Litter_size[!duplicated(df$Home_cage)])
sd_litter_size <- sd(df$Litter_size[!duplicated(df$Home_cage)])

## CREATE Z-SCORE FOR LITTER SIZE TO STANDARDIZE IT ##
df <- df %>% mutate(z_litter_size = (Litter_size - mean_litter_size) / sd_litter_size)

## LINEAR MIXED EFFECTS MODEL FORMULA (RHS) ##
formula_rhs <- " -1 + age_factor + age_factor:sex_flag + 
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
                 (1|ID) + (1|Home_cage)"


###################################################################
## GENERIC FUNCTION TO MODEL AND GENERATE DIFFERENT PLOTS        ##
###################################################################

myplotfunc<-function(structname,
                     formula_rhs,
                     GroupFilter=c("WT-Saline", "WT-MTX", "KO-Saline", "KO-MTX"),
                     PltType="longitudinal", #or, e.g., "AgeP28"
                     PltTitle="Structure Plot"){
  cform <- as.formula( paste0("`",structname,"`"," ~ ",formula_rhs) )
  clm <- lmer(cform, data = df)
  slm <- summary(clm)
  dfmut <- df %>% mutate(Group = paste(Genotype, Treatment, sep = "-"))
  ## EXTRACT RANDOM EFFECTS FOR ID AND HOME_CAGE ## 
  adjusted_vol <- df[,structname] - ranef(clm)$ID[as.character(df$ID), "(Intercept)"] - ranef(clm)$Home_cage[as.character(df$Home_cage), "(Intercept)"]
  dfmut <- dfmut %>% mutate(adjusted_vol = adjusted_vol)
  ## FILTER TO INCLUDE ONLY WT-SALINE AND WT-MTX ## 
  dfmut <- dfmut %>% filter(Group %in% GroupFilter)
  dfmut$Group <- factor(dfmut$Group,levels=GroupFilter) #input GroupFilter determines plot order
  if (grepl("^AgeP",PltType)){
    dfmut <- dfmut %>% filter(as.character(age_factor)==sub("AgeP","",PltType))
  } else if (grepl("^AgeGroupP",PltType)){
    #this is really the sampe as "AgeP" but formatted slightly differently
    dfmut <- dfmut %>% filter(as.character(age_factor)==sub("AgeGroupP","",PltType))
  }
  ## CALCULATE SUMMARY STATS ##
  summary_df <- dfmut %>%
    group_by(Age, Group, Sex) %>%
    summarise(
      mean_vol = mean(adjusted_vol, na.rm = TRUE),
      se_vol = sd(adjusted_vol, na.rm = TRUE) / sqrt(length(adjusted_vol)),
      ci_low = mean_vol - qt(0.975, df = length(adjusted_vol) - 1) * se_vol,
      ci_high = mean_vol + qt(0.975, df = length(adjusted_vol) - 1) * se_vol,
      .groups = 'drop'
    )
  ## DEFINE COLOURS ##
  my_colours_long <- c("WT-Saline" = "blue","WT-MTX" = "red","KO-Saline" = "black","KO-MTX" = "grey")
  my_colours_long <- my_colours_long[names(my_colours_long) %in% GroupFilter]
  ## DEFINE GROUP LABELS ##
  group_labels <- c("WT-Saline" = bquote(italic("Il6")^"+/+" ~ "Saline"),
                    "WT-MTX"    = bquote(italic("Il6")^"+/+" ~ "MTX"),
                    "KO-Saline" = bquote(italic("Il6")^"-/-" ~ "Saline"),
                    "KO-MTX"    = bquote(italic("Il6")^"-/-" ~ "MTX") )
  group_labels <- group_labels[names(group_labels) %in% GroupFilter]
  ## PLOT ##
  if (PltType=="longitudinal"){
    cplt <- ggplot(summary_df, aes(x = as.character(Age), y = mean_vol, group = Group, color = Group)) + 
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                    width = 0.3, position = position_dodge(width = 0.8), size = 1) +
      geom_point(size = 4, position = position_dodge(width = 0.8)) +
      geom_jitter(data = dfmut, aes(y = adjusted_vol, color = Group),
                  alpha = 0.4, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)) +
      geom_vline(xintercept = 1.5, linetype = "dashed", color = "black", size = 0.5) +
      facet_wrap(~Sex) +
      labs(
        title = PltTitle,
        x = "Age (Days)", 
        y = expression("Volume" ~ "(" * mm^3 * ")"),
        color = "Groups"
      ) +
      scale_color_manual(values = my_colours_long, labels = group_labels) +
      theme_minimal(base_size = 18) +
      theme(
        text = element_text(family = "Arial", size = 18),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, size = 0.7)
      )
  } else if (grepl("^AgeP",PltType)){
    cplt <- ggplot(summary_df, aes(x = as.character(Age), y = mean_vol, group = Group, color = Group)) + 
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                    width = 0.3, position = position_dodge(width = 0.8), size = 1) +  # Error bars for CI
      geom_point(size = 5, position = position_dodge(width = 0.8)) +  # Larger points for emphasis
      geom_jitter(data = dfmut, aes(y = adjusted_vol, color = Group),  
                  alpha = 0.4, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)) +  
      facet_wrap(~Sex, labeller = as_labeller(c(Female = "Female", Male = "Male")), strip.position = "bottom") +  
      labs(
        title = PltTitle,
        x = NULL,
        y = expression("Volume" ~ "(" * mm^3 * ")"),
        color = "Groups"
      ) +
      scale_color_manual(values = my_colours_long, labels = group_labels) +
      theme_minimal(base_size = 18) +  # Match base font size
      theme(
        text = element_text(family = "Arial", size = 18),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 16),
        strip.placement = "outside",
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, size = 0.7)
      )
  } else if (grepl("^AgeGroupP",PltType)){
    cplt <- ggplot(summary_df, aes(x = Group, y = mean_vol, color = Group)) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                    width = 0.3, size = 1, position = position_dodge(width = 0.8)) +
      geom_point(size = 5, position = position_dodge(width = 0.8)) +
      geom_jitter(data = dfmut, aes(y = adjusted_vol, color = Group),
                  alpha = 0.4, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)) +
      facet_wrap(~Sex, strip.position = "bottom") +
      labs(
        title = PltTitle,
        x = NULL,
        y = expression("Volume" ~ "(" * mm^3 * ")")
      ) +
      scale_color_manual(values = my_colours_long) +
      scale_x_discrete(labels = group_labels) +
      theme_minimal(base_size = 14) +
      theme(
        text = element_text(family = "Arial", size = 14),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "plain", size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.line.x = element_line(color = "black", size = 0.5)
      )
  }
  return(list("slm"=slm,"plt"=cplt))
}


#########################################################
## SUMMARY LONGITUDINAL VOLUME PLOTS (WT MTX)          ##
#########################################################

## GREY MATTER ##
gmresult<-myplotfunc("GM",formula_rhs,GroupFilter=c("WT-Saline", "WT-MTX"),PltType="longitudinal",PltTitle="Grey Matter")
slm_gm <- gmresult[["slm"]]
# SIGNIFICANT RESULTS FOR WT MICE =
# M_TREATMENT_FLAG_28 = P < 0.001 
# F_TREATMENT_FLAG_28 = P < 0.001 
# M_TREATMENT_FLAG_42 = P < 0.05
cairo_pdf("Figure3B.pdf",width=8.0,height=10.0)
print(gmresult[["plt"]]) 
dev.off()


## WHITE MATTER ##
wmresult<-myplotfunc("WM",formula_rhs,GroupFilter=c("WT-Saline", "WT-MTX"),PltType="longitudinal",PltTitle="White Matter")
slm_wm <- wmresult[["slm"]]
# M_TREATMENT_FLAG_28 = P < 0.001 
# F_TREATMENT_FLAG_28 = P < 0.01 
cairo_pdf("Figure3C.pdf",width=8.0,height=10.0)
print(wmresult[["plt"]]) 
dev.off()


## HIPPOCAMPUS ##
hpresult<-myplotfunc("HP",formula_rhs,GroupFilter=c("WT-Saline", "WT-MTX"),PltType="longitudinal",PltTitle="Hippocampus")
slm_hp <- hpresult[["slm"]]
# M_TREATMENT_FLAG_28 = P < 0.01 
# F_TREATMENT_FLAG_28 = P < 0.05
cairo_pdf("Figure4A.pdf",width=8.0,height=10.0)
print(hpresult[["plt"]]) 
dev.off()


## THALAMUS ##
thalresult<-myplotfunc("thalamus",formula_rhs,GroupFilter=c("WT-Saline", "WT-MTX"),PltType="longitudinal",PltTitle="Thalamus")
slm_thal <- thalresult[["slm"]]
cairo_pdf("Figure4B.pdf",width=8.0,height=10.0)
print(thalresult[["plt"]]) 
dev.off()



#########################################################
## SAMPLE P28 VOLUME PLOTS (WT MTX)                    ##
#########################################################

## PRIMARY VISUAL CORTEX: BINOCULAR AREA ##
pvc_b_result<-myplotfunc("Primary visual cortex: binocular area",formula_rhs,GroupFilter=c("WT-Saline", "WT-MTX"),
                         PltType="AgeP28",PltTitle="Primary Visual Cortex: Binocular Area (P28)")
slm_pvc_b <- pvc_b_result[["slm"]]
cairo_pdf("Figure5B.pdf",width=10.0,height=5.5)
print(pvc_b_result[["plt"]]) 
dev.off()


## PRIMARY VISUAL CORTEX: MONOCULAR AREA ##
pvc_m_result<-myplotfunc("Primary visual cortex: monocular area",formula_rhs,GroupFilter=c("WT-Saline", "WT-MTX"),
                         PltType="AgeP28",PltTitle="Primary Visual Cortex: Monocular Area (P28)")
slm_pvc_m <- pvc_m_result[["slm"]]
cairo_pdf("Figure5C.pdf",width=10.0,height=5.5)
print(pvc_m_result[["plt"]]) 
dev.off()


## SECONDARY VISUAL CORTEX: LATERAL AREA ##
svc_result<-myplotfunc("Secondary visual cortex: lateral area",formula_rhs,GroupFilter=c("WT-Saline", "WT-MTX"),
                       PltType="AgeP28",PltTitle="Secondary Visual Cortex: Lateral Area (P28)")
slm_svc <- svc_result[["slm"]]
cairo_pdf("Figure5D.pdf",width=10.0,height=5.5)
print(svc_result[["plt"]]) 
dev.off()


## PRIMARY SOMATOSENSORY CORTEX: HINDLIMB REGION ##
psc_hindlimb_result<-myplotfunc("Primary somatosensory cortex: hindlimb region",formula_rhs,GroupFilter=c("WT-Saline", "WT-MTX"),
                                PltType="AgeP28",PltTitle="Primary Somatosensory Cortex: Hindlimb Region (P28)")
slm_psc_hindlimb <- psc_hindlimb_result[["slm"]]
cairo_pdf("Figure5E.pdf",width=10.0,height=5.5)
print(psc_hindlimb_result[["plt"]]) 
dev.off()



#########################################################
## SAMPLE P28 VOLUME PLOTS (WT + KO MTX)               ##
#########################################################

## Hippocampus ##
hpresult<-myplotfunc("HP",formula_rhs,GroupFilter=c("WT-Saline", "WT-MTX", "KO-Saline", "KO-MTX"),
                     PltType="AgeGroupP28",PltTitle="Hippocampus (P28)")
slm_hp <- hpresult[["slm"]]
cairo_pdf("Figure8Atop.pdf",width=12.0,height=5.5)
print(hpresult[["plt"]]) 
dev.off()

## Parietal Region ##
prresult<-myplotfunc("PR",formula_rhs,GroupFilter=c("WT-Saline", "WT-MTX", "KO-Saline", "KO-MTX"),
                     PltType="AgeGroupP28",PltTitle="Parietal Region (P28)")
slm_pr <- prresult[["slm"]]
cairo_pdf("Figure8Abottom.pdf",width=12.0,height=5.5)
print(prresult[["plt"]]) 
dev.off()

#########################################################################
## RUN LINEAR MIXED EFFECTS MODELS FOR EACH STRUCTURE AND SAVE RESULTS ##
## CODE USED TO GENERATE ALLRES.RData ##
#########################################################################

if (file.exists("ALLRES.RData")){
  load("ALLRES.RData")
} else {
  ## GENERATE TEMPLATE OUTPUT FROM LMER
  cform <- as.formula( paste0("GM"," ~ ",formula_rhs) )  
  clm <- lmer(cform, data = df)
  slm <- summary(clm)
  allcols <- as.vector(outer(rownames(slm$coefficients), colnames(slm$coefficients), paste, sep="."))
  structure_list<-dt.struc$Structure
  allres <- data.frame(Structure=structure_list)
  allres[,allcols] <- NA
  ## LOOP THROUGH ALL STRUCTURES AND RUN LMER ## 
  for (i in 1:length(structure_list)){
    cstruct <- as.character(allres$Structure[i])
    cform <- as.formula( paste0("`",cstruct,"`"," ~ ",formula_rhs) )  
    clm <- lmer(cform, data = df)
    slm <- summary(clm)
    allres[i,allcols]<-as.vector(slm$coefficients)
    ## SAVE ALLRES DATA ## 
    save(allres, file = "ALLRES.RData")}
  
}

################################################################################################
## STRUCTUREWISE TREATMENT EFFECT IN WT MTX MICE                                              ##
## THIS IS FOR FIGURE 5A TO COMPARE MTX WT MICE (P28-P63) TO WT SALINE MICE                   ##
## USED TO CREATE THE BRAIN MAP IN FIGURE 5A                                                  ##
################################################################################################

## LOAD DATAFRAME IF RUNNING INTERACTIVELY FROM HERE: load("ALLRES.RData")

## RUN FDR ACROSS WT MTX MICE FROM P28-P63 ##
allres_p <- p.adjust(c(
  allres$`F_treatment_flag_28.Pr(>|t|)`, allres$`M_treatment_flag_28.Pr(>|t|)`,
  allres$`F_treatment_flag_42.Pr(>|t|)`, allres$`M_treatment_flag_42.Pr(>|t|)`,
  allres$`F_treatment_flag_63.Pr(>|t|)`, allres$`M_treatment_flag_63.Pr(>|t|)`
), "fdr")

allres$`F_treatment_flag_28.FDR` <- allres_p[1:nrow(allres)]
allres$`M_treatment_flag_28.FDR` <- allres_p[(nrow(allres) + 1):(2 * nrow(allres))]
allres$`F_treatment_flag_42.FDR` <- allres_p[(2 * nrow(allres) + 1):(3 * nrow(allres))]
allres$`M_treatment_flag_42.FDR` <- allres_p[(3 * nrow(allres) + 1):(4 * nrow(allres))]
allres$`F_treatment_flag_63.FDR` <- allres_p[(4 * nrow(allres) + 1):(5 * nrow(allres))]
allres$`M_treatment_flag_63.FDR` <- allres_p[(5 * nrow(allres) + 1):(6 * nrow(allres))]

## FOR INTERPRETABILITY, COMPUTE PERCENT CHANGE IN STRUCTURE VOLS AT EACH TIME POINT ##

## FOR WT FEMALE: USING AGE_FACTOR28 B/C F ARE BASELINES ## 
F_WT_MTX_28_pc <- allres$F_treatment_flag_28.Estimate/allres$age_factor28.Estimate*100
F_WT_MTX_42_pc <- allres$F_treatment_flag_42.Estimate/allres$age_factor42.Estimate*100
F_WT_MTX_63_pc <- allres$F_treatment_flag_63.Estimate/allres$age_factor63.Estimate*100

## FOR WT MALE: ADDING SEX FLAG TO AGE FACTOR TO GET MALE BASELINE MALE ##
Male_WT_Saline_Estimate_28 <- allres$'age_factor28:sex_flag.Estimate' + allres$'age_factor28.Estimate'
M_WT_MTX_28_pc <- allres$M_treatment_flag_28.Estimate/Male_WT_Saline_Estimate_28*100
Male_WT_Saline_Estimate_42 <- allres$'age_factor42:sex_flag.Estimate' + allres$'age_factor42.Estimate'
M_WT_MTX_42_pc <- allres$M_treatment_flag_42.Estimate/Male_WT_Saline_Estimate_42*100
Male_WT_Saline_Estimate_63 <- allres$'age_factor63:sex_flag.Estimate' + allres$'age_factor63.Estimate'
M_WT_MTX_63_pc <- allres$M_treatment_flag_63.Estimate/Male_WT_Saline_Estimate_63*100

## SUMMARIZE # OF AFFECTED STRUCTURES ##
print( paste0("MTX VOL CHANGE (FEMALE, P28, FDR<0.1): ",as.character(sum( allres$`F_treatment_flag_28.FDR` < 0.1 ))," structures") ) #34
print( paste0("MTX VOL CHANGE (MALE, P28, FDR<0.1): ",as.character(sum( allres$`M_treatment_flag_28.FDR` < 0.1 ))," structures") )   #59
print( paste0("MTX VOL CHANGE (FEMALE, P42, FDR<0.1): ",as.character(sum( allres$`F_treatment_flag_42.FDR` < 0.1 ))," structures") ) #4
print( paste0("MTX VOL CHANGE (MALE, P42, FDR<0.1): ",as.character(sum( allres$`M_treatment_flag_42.FDR` < 0.1 ))," structures") )   #8
print( paste0("MTX VOL CHANGE (FEMALE, P63, FDR<0.1): ",as.character(sum( allres$`F_treatment_flag_63.FDR` < 0.1 ))," structures") ) #0
print( paste0("MTX VOL CHANGE (MALE, P63, FDR<0.1): ",as.character(sum( allres$`M_treatment_flag_63.FDR` < 0.1 ))," structures") )   #0


#############################################################################################
## STRUCTUREWISE GENOTYPE EFFECT IN KO SALINE MICE                                         ##
## THIS IS FOR FIGURE 6A TO COMPARE KO SALINE MICE (P14-P63) TO WT SALINE MICE             ##
## USED TO CREATE THE BRAIN MAP IN FIGURE 6A                                               ##
#############################################################################################

## LOAD DATAFRAME IF RUNNING INTERACTIVELY FROM HERE: load("ALLRES.RData")

## RUN FDR across all genotype/age/sex groups ##
allres_p <- p.adjust(c(
  allres$`F_genotype_flag_14.Pr(>|t|)`, allres$`F_genotype_flag_28.Pr(>|t|)`,
  allres$`F_genotype_flag_42.Pr(>|t|)`, allres$`F_genotype_flag_63.Pr(>|t|)`,
  allres$`M_genotype_flag_14.Pr(>|t|)`, allres$`M_genotype_flag_28.Pr(>|t|)`,
  allres$`M_genotype_flag_42.Pr(>|t|)`, allres$`M_genotype_flag_63.Pr(>|t|)`
), method = "fdr")

allres$`F_genotype_flag_14.FDR` <- allres_p[1:nrow(allres)]
allres$`F_genotype_flag_28.FDR` <- allres_p[(nrow(allres) + 1):(2 * nrow(allres))]
allres$`F_genotype_flag_42.FDR` <- allres_p[(2 * nrow(allres) + 1):(3 * nrow(allres))]
allres$`F_genotype_flag_63.FDR` <- allres_p[(3 * nrow(allres) + 1):(4 * nrow(allres))]
allres$`M_genotype_flag_14.FDR` <- allres_p[(4 * nrow(allres) + 1):(5 * nrow(allres))]
allres$`M_genotype_flag_28.FDR` <- allres_p[(5 * nrow(allres) + 1):(6 * nrow(allres))]
allres$`M_genotype_flag_42.FDR` <- allres_p[(6 * nrow(allres) + 1):(7 * nrow(allres))]
allres$`M_genotype_flag_63.FDR` <- allres_p[(7 * nrow(allres) + 1):(8 * nrow(allres))]

## FOR INTERPRETABILITY, COMPUTE PERCENT CHANGE IN STRUCTURE VOLS AT EACH TIME POINT ##

## % CHANGE FOR FEMALE KO MICE ##
F_KO_14_pc <- allres$F_genotype_flag_14.Estimate / allres$age_factor14.Estimate * 100
F_KO_28_pc <- allres$F_genotype_flag_28.Estimate / allres$age_factor28.Estimate * 100
F_KO_42_pc <- allres$F_genotype_flag_42.Estimate / allres$age_factor42.Estimate * 100
F_KO_63_pc <- allres$F_genotype_flag_63.Estimate / allres$age_factor63.Estimate * 100

## FOR WT MALE: ADDING SEX FLAG TO AGE FACTOR TO GET MALE BASELINE MALE ##
Male_WT_Saline_Estimate_14 <- allres$`age_factor14:sex_flag.Estimate` + allres$age_factor14.Estimate
M_KO_14_pc <- allres$M_genotype_flag_14.Estimate / Male_WT_Saline_Estimate_14 * 100
Male_WT_Saline_Estimate_28 <- allres$`age_factor28:sex_flag.Estimate` + allres$age_factor28.Estimate
M_KO_28_pc <- allres$M_genotype_flag_28.Estimate / Male_WT_Saline_Estimate_28 * 100
Male_WT_Saline_Estimate_42 <- allres$`age_factor42:sex_flag.Estimate` + allres$age_factor42.Estimate
M_KO_42_pc <- allres$M_genotype_flag_42.Estimate / Male_WT_Saline_Estimate_42 * 100
Male_WT_Saline_Estimate_63 <- allres$`age_factor63:sex_flag.Estimate` + allres$age_factor63.Estimate
M_KO_63_pc <- allres$M_genotype_flag_63.Estimate / Male_WT_Saline_Estimate_63 * 100

## SUMMARIZE # OF AFFECTED STRUCTURES ##
print( paste0("IL6KO SALINE VOL CHANGE (FEMALE, P14, FDR<0.1): ",as.character(sum( allres$`F_genotype_flag_14.FDR` < 0.1 ))," structures") ) #0
print( paste0("IL6KO SALINE VOL CHANGE (FEMALE, P28, FDR<0.1): ",as.character(sum( allres$`F_genotype_flag_28.FDR` < 0.1 ))," structures") ) #1
print( paste0("IL6KO SALINE VOL CHANGE (FEMALE, P42, FDR<0.1): ",as.character(sum( allres$`F_genotype_flag_42.FDR` < 0.1 ))," structures") ) #0
print( paste0("IL6KO SALINE VOL CHANGE (FEMALE, P63, FDR<0.1): ",as.character(sum( allres$`F_genotype_flag_63.FDR` < 0.1 ))," structures") ) #1
print( paste0("IL6KO SALINE VOL CHANGE (MALE, P14, FDR<0.1): ",as.character(sum( allres$`M_genotype_flag_14.FDR` < 0.1 ))," structures") ) #5
print( paste0("IL6KO SALINE VOL CHANGE (MALE, P28, FDR<0.1): ",as.character(sum( allres$`M_genotype_flag_28.FDR` < 0.1 ))," structures") ) #0
print( paste0("IL6KO SALINE VOL CHANGE (MALE, P42, FDR<0.1): ",as.character(sum( allres$`M_genotype_flag_42.FDR` < 0.1 ))," structures") ) #0
print( paste0("IL6KO SALINE VOL CHANGE (MALE, P63, FDR<0.1): ",as.character(sum( allres$`M_genotype_flag_63.FDR` < 0.1 ))," structures") ) #0



##################################################################################
## GENERATING GENOTYPE EFFECT HISTOGRAM PLOTS                                   ##
## THIS IS FOR FIGURE 6B TO COMPARE KO SALINE MICE (P14-P63) TO WT SALINE MICE  ##
##################################################################################

## LOAD DATAFRAME IF RUNNING INTERACTIVELY FROM HERE: load("ALLRES.RData")

## BASELINE ESTIMATES FOR FEMALES ## 
age_factor_14 <- allres$`age_factor14.Estimate`
age_factor_28 <- allres$`age_factor28.Estimate`
age_factor_42 <- allres$`age_factor42.Estimate`
age_factor_63 <- allres$`age_factor63.Estimate`

## MALE BASLEINES (WT+SALINE) FOR EACH AGE ## FEMALES ARE BASLINE SO HAD TO ADD SEX:FLAG TO GET MALE BASLINE ##
Male_WT_Saline_Estimate_14 <- age_factor_14 + allres$`age_factor14:sex_flag.Estimate`
Male_WT_Saline_Estimate_28 <- age_factor_28 + allres$`age_factor28:sex_flag.Estimate`
Male_WT_Saline_Estimate_42 <- age_factor_42 + allres$`age_factor42:sex_flag.Estimate`
Male_WT_Saline_Estimate_63 <- age_factor_63 + allres$`age_factor63:sex_flag.Estimate`

## PERCENT CHANGE FOR EACH TERM ##

## FEMALE ##
f_estimates_14 <- (allres$`F_genotype_flag_14.Estimate` / age_factor_14) * 100
f_estimates_28 <- (allres$`F_genotype_flag_28.Estimate` / age_factor_28) * 100
f_estimates_42 <- (allres$`F_genotype_flag_42.Estimate` / age_factor_42) * 100
f_estimates_63 <- (allres$`F_genotype_flag_63.Estimate` / age_factor_63) * 100

## MALE ##
m_estimates_14 <- (allres$`M_genotype_flag_14.Estimate` / Male_WT_Saline_Estimate_14) * 100
m_estimates_28 <- (allres$`M_genotype_flag_28.Estimate` / Male_WT_Saline_Estimate_28) * 100
m_estimates_42 <- (allres$`M_genotype_flag_42.Estimate` / Male_WT_Saline_Estimate_42) * 100
m_estimates_63 <- (allres$`M_genotype_flag_63.Estimate` / Male_WT_Saline_Estimate_63) * 100

# Combine female and male data with Genotype labels and Sex column
f_data <- data.frame(
  Estimate = c(f_estimates_14, f_estimates_28, f_estimates_42, f_estimates_63),
  Genotype = factor(rep(c("P14 KO", "P28 KO", "P42 KO", "P63 KO"), each = length(f_estimates_14))),
  Sex = "Female"
)

m_data <- data.frame(
  Estimate = c(m_estimates_14, m_estimates_28, m_estimates_42, m_estimates_63),
  Genotype = factor(rep(c("P14 KO", "P28 KO", "P42 KO", "P63 KO"), each = length(m_estimates_14))),
  Sex = "Male"
)

combined_data <- rbind(f_data, m_data)

# Set x-axis limits for consistent alignment of zero
x_min <- min(combined_data$Estimate) - 1
x_max <- max(combined_data$Estimate) + 1

# Define colors for each genotype
genotype_colors <- c("P14 KO" = "#B2DF8A", "P28 KO" = "#A6CEE3", "P42 KO" = "#33A02C", "P63 KO" = "#1F78B4")

cplt<-ggplot(combined_data, aes(x = Estimate, fill = Genotype)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.8, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(values = genotype_colors) +
  labs(
    x = "% Volume Change (Relative to Age- and Sex-Matched Controls)",
    y = "Number of Brain Structures"
  ) +
  xlim(x_min, x_max) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = 12, family = "Arial"),
    axis.title.y = element_text(size = 12, family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold", family = "Arial"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(1, "lines")
  ) +
  facet_grid(Sex ~ Genotype, scales = "free_y")

cairo_pdf("Figure6B.pdf",width=10.0,height=5.5)
print(cplt) 
dev.off()



########################################################################################
## GENERATE 1000 PERMUTATIONS TO TEST GENOTYPE EFFECTS                                ##
## THIS SCRIPT SHUFFLES GENOTYPE WITHIN SEX AND RE-RUNS MODELS FOR EACH STRUCTURE     ##
## OUTPUT: PERMUTATION_RESULTS_1000.RData FOR STATISTICAL COMPARISON (FIGURE 6)       ##
#########################################################################################


if (file.exists("PERMUTATION_RESULTS_1000.RData")){
  load("PERMUTATION_RESULTS_1000.RData") 
} else {
  ## SET UP EMPTY RESULTS DATA FRAME ##
  allres$Permutation <- NA
  
  ## PERMUTATION PARAMETERS ##
  n_permutation <- 1000
  perm_results_df <- data.frame()  
  
  ## GET UNIQUE ID-GENOTYPE-SEX MAPPING ##
  df_ID_gen <- distinct(df, ID, Genotype, Sex, .keep_all = FALSE)
  
  ## REMOVE GENOTYPE FOR PERMUTATION ##
  df_no_gen <- select(df, -Genotype)
  
  ## BEGIN PERMUTATION LOOP ##
  set.seed(123)
  for (p in 1:n_permutation) {
    
    df_dummy <- allres
    
    ## SPLIT BY SEX ##
    df_ID_F <- subset(df_ID_gen, Sex == "Female")
    df_ID_M <- subset(df_ID_gen, Sex == "Male")
    
    ## SHUFFLE GENOTYPE ##
    df_shuffle_F <- transform(df_ID_F, Genotype = sample(df_ID_F$Genotype))
    df_shuffle_M <- transform(df_ID_M, Genotype = sample(df_ID_M$Genotype))
    df_shuffle <- rbind(df_shuffle_F, df_shuffle_M)
    df_shuffle <- df_shuffle[c("ID", "Genotype")]
    
    ## MERGE BACK SHUFFLED GENOTYPE ##
    df_perm <- left_join(df_no_gen, df_shuffle, by = "ID")
    
    ## REDEFINE FLAGS WITH SHUFFLED GENOTYPE ##
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
    
    ## LOOP THROUGH ALL STRUCTURES ##
    for (j in 1:length(structure_list)) {
      cstruct <- as.character(structure_list[j])
      cform <- as.formula( paste0("`",cstruct,"`"," ~ ",formula_rhs) )  
      clm_perm <- lmer(cform, data = df_perm)
      slm <- summary(clm_perm)    
      df_dummy[j, allcols] <- as.vector(slm$coefficients)
      df_dummy$Permutation <- p
    }
    
    ## COMBINE ALL RESULTS INTO ONE BIG DATA FRAME ##
    perm_results_df <- rbind(perm_results_df, df_dummy)
  }
  
  ## SAVE RESULTS ##
  save(perm_results_df, file = "PERMUTATION_RESULTS_1000.RData")
  
}


####################################################################################################
## COMPARE OBSERVED GENOTYPE EFFECTS TO PERMUTATION DISTRIBUTIONS USING EMPIRICAL P-VALUES       ##
## NORMALIZE ESTIMATES BY SEX AND AGE, COMPUTE MEDIANS, AND CALCULATE TWO-TAILED P-VALUES        ##
## DATA SOURCES: PERMUTATION_RESULTS_1000.RData & ALLRES.RData                                   ##
####################################################################################################

## NORMALIZE ESTIMATES FOR EACH AGE AND SEX ## 
perm_results_df <- perm_results_df %>%
  mutate(
    ## FEMALE ##
    f_estimates_14_p = (`F_genotype_flag_14.Estimate` / `age_factor14.Estimate`) * 100,
    f_estimates_28_p = (`F_genotype_flag_28.Estimate` / `age_factor28.Estimate`) * 100,
    f_estimates_42_p = (`F_genotype_flag_42.Estimate` / `age_factor42.Estimate`) * 100,
    f_estimates_63_p = (`F_genotype_flag_63.Estimate` / `age_factor63.Estimate`) * 100,
    
    ## MALE: COMPUTE NEW BASELINE FIRST, THEN NORMALIZE ## 
    Male_WT_Saline_Estimate_14_p = `age_factor14:sex_flag.Estimate` + `age_factor14.Estimate`,
    Male_WT_Saline_Estimate_28_p = `age_factor28:sex_flag.Estimate` + `age_factor28.Estimate`,
    Male_WT_Saline_Estimate_42_p = `age_factor42:sex_flag.Estimate` + `age_factor42.Estimate`,
    Male_WT_Saline_Estimate_63_p = `age_factor63:sex_flag.Estimate` + `age_factor63.Estimate`,
    
    m_estimates_14_p = (`M_genotype_flag_14.Estimate` / Male_WT_Saline_Estimate_14_p) * 100,
    m_estimates_28_p = (`M_genotype_flag_28.Estimate` / Male_WT_Saline_Estimate_28_p) * 100,
    m_estimates_42_p = (`M_genotype_flag_42.Estimate` / Male_WT_Saline_Estimate_42_p) * 100,
    m_estimates_63_p = (`M_genotype_flag_63.Estimate` / Male_WT_Saline_Estimate_63_p) * 100
  )

## COMPUTE MEDIAN FOR EACH PERMUTATION ## 
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

# print(perm_medians)

observed_F_med_14 <- median(f_estimates_14)
observed_F_med_28 <- median(f_estimates_28)
observed_F_med_42 <- median(f_estimates_42)
observed_F_med_63 <- median(f_estimates_63)

observed_M_med_14 <- median(m_estimates_14)
observed_M_med_28 <- median(m_estimates_28)
observed_M_med_42 <- median(m_estimates_42)
observed_M_med_63 <- median(m_estimates_63)

## STATISTICAL TESTING USING TWO-TAILED EMPERICAL P-VALUE ##

## FEMALES ##

# P14 #
N <- length(perm_medians$Median_F_Gen_14_p)    #should match n_permutation
j <- sum(observed_F_med_14 > perm_medians$Median_F_Gen_14_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median IL6KO Effect across all structures (FEMALE, P14): ",as.character(observed_F_med_14),"p=",as.character(mypval)) )  #0.808 

# P28 #
N <- length(perm_medians$Median_F_Gen_28_p)
j <- sum(observed_F_med_28 > perm_medians$Median_F_Gen_28_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median IL6KO Effect across all structures (FEMALE, P28): ",as.character(observed_F_med_28),"p=",as.character(mypval)) )  #0.44

# P42 #
N <- length(perm_medians$Median_F_Gen_42_p)
j <- sum(observed_F_med_42 > perm_medians$Median_F_Gen_42_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median IL6KO Effect across all structures (FEMALE, P42): ",as.character(observed_F_med_42),"p=",as.character(mypval)) )  #0.214

# P63 #
N <- length(perm_medians$Median_F_Gen_63_p)
j <- sum(observed_F_med_63 > perm_medians$Median_F_Gen_63_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median IL6KO Effect across all structures (FEMALE, P63): ",as.character(observed_F_med_63),"p=",as.character(mypval)) )  #0.2

## MALES ##

# P14 #
N <- length(perm_medians$Median_M_Gen_14_p)
j <- sum(observed_M_med_14 > perm_medians$Median_M_Gen_14_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median IL6KO Effect across all structures (MALE, P14): ",as.character(observed_M_med_14),"p=",as.character(mypval)) )  #0.082 

# P28 #
N <- length(perm_medians$Median_M_Gen_28_p)
j <- sum(observed_M_med_28 > perm_medians$Median_M_Gen_28_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median IL6KO Effect across all structures (MALE, P28): ",as.character(observed_M_med_28),"p=",as.character(mypval)) )  #0.274 

# P42 #
N <- length(perm_medians$Median_M_Gen_42_p)
j <- sum(observed_M_med_42 > perm_medians$Median_M_Gen_42_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median IL6KO Effect across all structures (MALE, P42): ",as.character(observed_M_med_42),"p=",as.character(mypval)) )  #0.385

# P63 #
N <- length(perm_medians$Median_M_Gen_63_p)
j <- sum(observed_M_med_63 > perm_medians$Median_M_Gen_63_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median IL6KO Effect across all structures (MALE, P63): ",as.character(observed_M_med_63),"p=",as.character(mypval)) )  #0.614

########################################################################################
## USE A SECOND LINEAR MODEL FOR CONVENIENCE TO PROBE IL6KO-MTX VS IL6KO-SALINE
## SEPARATED BY AGE AND SEX
##
## RUN LINEAR MIXED EFFECTS MODEL #2 FOR EACH STRUCTURE AND SAVE RESULTS ##
## CODE USED TO GENERATE ALLRES_2.RData ##
#########################################################################################

M_KO_flag_14 <- (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 14) * (FINAL_DATA_FRAME$Genotype == "KO")
M_WT_flag_14 <- (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 14) * (FINAL_DATA_FRAME$Genotype != "KO")
F_KO_flag_14 <- (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 14) * (FINAL_DATA_FRAME$Genotype == "KO")
F_WT_flag_14 <- (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 14) * (FINAL_DATA_FRAME$Genotype != "KO")

M_MTX_KO_flag_28 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 28) * (FINAL_DATA_FRAME$Genotype == "KO")
M_MTX_WT_flag_28 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 28) * (FINAL_DATA_FRAME$Genotype != "KO")
F_MTX_KO_flag_28 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 28) * (FINAL_DATA_FRAME$Genotype == "KO")
F_MTX_WT_flag_28 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 28) * (FINAL_DATA_FRAME$Genotype != "KO")
M_Saline_KO_flag_28 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 28) * (FINAL_DATA_FRAME$Genotype == "KO")
M_Saline_WT_flag_28 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 28) * (FINAL_DATA_FRAME$Genotype != "KO")
F_Saline_KO_flag_28 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 28) * (FINAL_DATA_FRAME$Genotype == "KO")
F_Saline_WT_flag_28 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 28) * (FINAL_DATA_FRAME$Genotype != "KO")

M_MTX_KO_flag_42 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 42) * (FINAL_DATA_FRAME$Genotype == "KO")
M_MTX_WT_flag_42 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 42) * (FINAL_DATA_FRAME$Genotype != "KO")
F_MTX_KO_flag_42 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 42) * (FINAL_DATA_FRAME$Genotype == "KO")
F_MTX_WT_flag_42 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 42) * (FINAL_DATA_FRAME$Genotype != "KO")
M_Saline_KO_flag_42 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 42) * (FINAL_DATA_FRAME$Genotype == "KO")
M_Saline_WT_flag_42 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 42) * (FINAL_DATA_FRAME$Genotype != "KO")
F_Saline_KO_flag_42 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 42) * (FINAL_DATA_FRAME$Genotype == "KO")
F_Saline_WT_flag_42 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 42) * (FINAL_DATA_FRAME$Genotype != "KO")

M_MTX_KO_flag_63 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 63) * (FINAL_DATA_FRAME$Genotype == "KO")
M_MTX_WT_flag_63 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 63) * (FINAL_DATA_FRAME$Genotype != "KO")
F_MTX_KO_flag_63 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 63) * (FINAL_DATA_FRAME$Genotype == "KO")
F_MTX_WT_flag_63 <- (FINAL_DATA_FRAME$Treatment == "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 63) * (FINAL_DATA_FRAME$Genotype != "KO")
M_Saline_KO_flag_63 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 63) * (FINAL_DATA_FRAME$Genotype == "KO")
M_Saline_WT_flag_63 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Male") * (FINAL_DATA_FRAME$Age == 63) * (FINAL_DATA_FRAME$Genotype != "KO")
F_Saline_KO_flag_63 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 63) * (FINAL_DATA_FRAME$Genotype == "KO")
F_Saline_WT_flag_63 <- (FINAL_DATA_FRAME$Treatment != "MTX") * (FINAL_DATA_FRAME$Sex == "Female") * (FINAL_DATA_FRAME$Age == 63) * (FINAL_DATA_FRAME$Genotype != "KO")

formula_rhs_2 <- " -1 + age_factor + age_factor:sex_flag + 
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
                    (1|ID) + (1|Home_cage)"


if (file.exists("ALLRES_2.RData")){
  load("ALLRES_2.RData")
} else {
  ## GENERATE TEMPLATE OUTPUT FROM LMER
  cform <- as.formula( paste0("GM"," ~ ",formula_rhs_2) )  
  clm <- lmer(cform, data = df)
  slm <- summary(clm)
  allcols <- as.vector(outer(rownames(slm$coefficients), colnames(slm$coefficients), paste, sep="."))
  allres_2 <- data.frame(Structure=structure_list)
  allres_2[,allcols] <- NA
  ## LOOP THROUGH ALL STRUCTURES AND RUN LMER ## 
  for (i in 1:length(structure_list)){
    cstruct <- as.character(allres_2$Structure[i])
    cform <- as.formula( paste0("`",cstruct,"`"," ~ ",formula_rhs_2) )  
    clm <- lmer(cform, data = df)
    slm <- summary(clm)
    allres_2[i,allcols]<-as.vector(slm$coefficients)}
  ## SAVE ALLRES_2.RData ##
  save(allres_2, file = "ALLRES_2.RData")
}

#################################################################################################################
## STRUCTUREWISE TREATMENT EFFECT IN MTX MICE AT P28                                                           ##
## THIS IS FOR FIGURE 7A TO COMPARE 1) WT MTX MICE TO WT SALINE 2) KO MTX TO KO SALINE (SPLIT BY SEX AT P28)   ##
## USED TO CREATE THE BRAIN MAP IN FIGURE 7A                                                                   ##
#################################################################################################################

## LOAD DATAFRAME IF RUNNING INTERACTIVELY FROM HERE: load("ALLRES_2.RData")

## RUN FDR correction across MTX comparisons at P28 ##
allres_p <- p.adjust(c(
  allres_2$`F_MTX_WT_flag_28.Pr(>|t|)`,
  allres_2$`M_MTX_WT_flag_28.Pr(>|t|)`,
  allres_2$`F_MTX_KO_flag_28.Pr(>|t|)`,
  allres_2$`M_MTX_KO_flag_28.Pr(>|t|)`
), method = "fdr")

allres_2$`F_MTX_WT_flag_28.FDR` <- allres_p[1:nrow(allres_2)]
allres_2$`M_MTX_WT_flag_28.FDR` <- allres_p[(nrow(allres_2) + 1):(2 * nrow(allres_2))]
allres_2$`F_MTX_KO_flag_28.FDR` <- allres_p[(2 * nrow(allres_2) + 1):(3 * nrow(allres_2))]
allres_2$`M_MTX_KO_flag_28.FDR` <- allres_p[(3 * nrow(allres_2) + 1):(4 * nrow(allres_2))]

## FOR INTERPRETABILITY, COMPUTE PERCENT CHANGE IN STRUCTURE VOLS AT P28 ##

## FOR WT FEMALE ##
F_WT_MTX_28_pc <- allres_2$F_MTX_WT_flag_28.Estimate / allres_2$age_factor28.Estimate * 100

## FOR WT MALE ##
Male_WT_Saline_Estimate_28 <- allres_2$`age_factor28:sex_flag.Estimate` + allres_2$age_factor28.Estimate
M_WT_MTX_28_pc <- allres_2$M_MTX_WT_flag_28.Estimate / Male_WT_Saline_Estimate_28 * 100

## FOR KO FEMALE ##
F_KO_Saline_28_Estimate <- allres_2$`F_Saline_KO_flag_28.Estimate` + allres_2$age_factor28.Estimate
F_KO_MTX_28_pc <- allres_2$F_MTX_KO_flag_28.Estimate / F_KO_Saline_28_Estimate * 100

## FOR KO MALE ##
M_KO_Saline_28_Estimate <- allres_2$`M_Saline_KO_flag_28.Estimate` + Male_WT_Saline_Estimate_28
M_KO_MTX_28_pc <- allres_2$M_MTX_KO_flag_28.Estimate / M_KO_Saline_28_Estimate * 100

## SUMMARIZE # OF AFFECTED STRUCTURES ##
print(paste0("MTX EFFECT (FEMALE WT, P28, FDR<0.1): ", sum(allres_2$`F_MTX_WT_flag_28.FDR` < 0.1), " structures")) #56
print(paste0("MTX EFFECT (MALE WT, P28, FDR<0.1): ", sum(allres_2$`M_MTX_WT_flag_28.FDR` < 0.1), " structures")) #79
print(paste0("MTX EFFECT (FEMALE KO, P28, FDR<0.1): ", sum(allres_2$`F_MTX_KO_flag_28.FDR` < 0.1), " structures")) #75
print(paste0("MTX EFFECT (MALE KO, P28, FDR<0.1): ", sum(allres_2$`M_MTX_KO_flag_28.FDR` < 0.1), " structures")) #44

######################################################################################
## GENERATING TREATMENT EFFECT HISTOGRAM PLOTS                                      ##
## THIS IS FOR FIGURE 7B TO COMPARE MTX-TREATED MICE TO SALINE-TREATED MICE AT P28  ##
## COMPARING WT MTX TO WT SALINE & KO MTX TO KO SALINE (SEPARATED BY SEX)           ##
######################################################################################

## LOAD DATAFRAME IF RUNNING INTERACTIVELY FROM HERE: load("ALLRES_2.RData")

## WT FEMALE ##
F_WT_MTX_28_pc <- allres_2$F_MTX_WT_flag_28.Estimate / allres_2$age_factor28.Estimate * 100

## WT MALE ##
Male_WT_Saline_Estimate_28 <- allres_2$'age_factor28:sex_flag.Estimate' + allres_2$'age_factor28.Estimate'
M_WT_MTX_28_pc <- allres_2$M_MTX_WT_flag_28.Estimate / Male_WT_Saline_Estimate_28 * 100

## KO FEMALE ##
F_KO_Saline_28_Estimate <- allres_2$'F_Saline_KO_flag_28.Estimate' + allres_2$'age_factor28.Estimate'
F_KO_MTX_28_pc <- allres_2$F_MTX_KO_flag_28.Estimate / F_KO_Saline_28_Estimate * 100

## KO MALE ##
M_KO_Saline_28_Estimate <- allres_2$'M_Saline_KO_flag_28.Estimate' + Male_WT_Saline_Estimate_28
M_KO_MTX_28_pc <- allres_2$M_MTX_KO_flag_28.Estimate / M_KO_Saline_28_Estimate * 100

## CREATE COMBINED_DATA WITH SEX AND GENOTYPE ##
combined_data <- data.frame(
  Estimate = c(F_WT_MTX_28_pc, M_WT_MTX_28_pc, F_KO_MTX_28_pc, M_KO_MTX_28_pc),
  Sex = rep(c("Female", "Male", "Female", "Male"), each = nrow(allres_2)),
  Genotype = rep(c("WT", "WT", "KO", "KO"), each = nrow(allres_2))
)

## SET AXIS LIMITS FOR CONSISTENT X-AXIS ##
x_min <- min(combined_data$Estimate, na.rm = TRUE)
x_max <- max(combined_data$Estimate, na.rm = TRUE)

## SET CUSTOM COLOURS ##
sex_colors <- c("Female" = "#66CC99", "Male" = "#FF9966")

## ENSURE CORRECT ORDERING OF FACTOR LEVELS ##
combined_data$Sex <- factor(combined_data$Sex, levels = c("Female", "Male"))
combined_data$Genotype <- factor(combined_data$Genotype, levels = c("WT", "KO"))

## PLOT ##
cplt<-ggplot(combined_data, aes(x = Estimate, fill = Sex)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.8, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(values = sex_colors) +
  labs(
    x = "% Volume Change (Relative to P28 Sex- and Genotype-Matched Saline Controls)",
    y = "Number of Brain Structures"
  ) +
  xlim(x_min, x_max) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14, family = "Arial"),
    axis.title.y = element_text(size = 14, family = "Arial"),
    axis.text = element_text(size = 14, family = "Arial"),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold", family = "Arial"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(1, "lines")
  ) +
  facet_grid(Genotype ~ Sex, scales = "free_y")

cairo_pdf("Figure7B.pdf",width=10.0,height=5.5)
print(cplt) 
dev.off()

###############################################################################################
## GENERATE 1000 PERMUTATIONS TO TEST TREATMENT EFFECTS                                      ##
## THIS SCRIPT SHUFFLES TREATMENT WITHIN SEX AND RE-RUNS MODELS FOR EACH STRUCTURE           ##
## OUTPUT: PERMUTATION_RESULTS_TREATMENT_1000.RData FOR STATISTICAL COMPARISON (FIGURE 7)    ##
###############################################################################################

if (file.exists("PERMUTATION_RESULTS_TREATMENT_1000.RData")){
  load("PERMUTATION_RESULTS_TREATMENT_1000.RData") 
} else {
  ## SET UP EMPTY RESULTS DATA FRAME ##
  allres_2$Permutation <- NA
  
  ## PERMUTATION PARAMETERS ##
  n_permutation <- 1000
  perm_treat_df <- data.frame() 
  
  ## GET UNIQUE ID-GENOTYPE-SEX MAPPING ##
  df_ID_treat <- distinct(df, ID, Treatment, Sex, .keep_all = FALSE)
  
  ## REMOVE TREATMENT FOR PERMUTATION ##
  df_no_treat <- select(df, -Treatment)
  
  ## BEGIN PERMUTATION LOOP ## 
  set.seed(123)
  for (p in 1:n_permutation) {
    
    df_dummy_treat <- allres_2
    
    ## SPLIT BY SEX ##
    df_ID_F_TREAT <- subset(df_ID_treat, Sex=="Female")
    df_ID_M_TREAT <- subset(df_ID_treat, Sex=="Male")
    
    ## SHUFFLE TREATMENT ##
    df_shuffle_F_TREAT <- transform(df_ID_F_TREAT, Treatment = sample(df_ID_F_TREAT$Treatment))
    df_shuffle_M_TREAT <- transform(df_ID_M_TREAT, Treatment = sample(df_ID_M_TREAT$Treatment))
    df_shuffle_TREAT <- rbind(df_shuffle_F_TREAT, df_shuffle_M_TREAT)
    df_shuffle_TREAT <- df_shuffle_TREAT[c("ID", "Treatment")]
    
    ## MERGE BACK SHUFFLED TREATMENT ##
    df_perm_treat <- left_join(df_no_treat, df_shuffle_TREAT, by = "ID")
    
    ## REDEFINE FLAGS WITH SHUFFLED TREATMENT ##
    age_factor <- as.character(df_perm_treat$Age)
    sex_flag <- (df_perm_treat$Sex == "Male") * 1
    
    M_KO_flag_14 <- (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 14) * (df_perm_treat$Genotype == "KO")
    M_WT_flag_14 <- (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 14) * (df_perm_treat$Genotype != "KO")
    F_KO_flag_14 <- (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 14) * (df_perm_treat$Genotype == "KO")
    F_WT_flag_14 <- (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 14) * (df_perm_treat$Genotype != "KO")
    
    M_MTX_KO_flag_28 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 28) * (df_perm_treat$Genotype == "KO")
    M_MTX_WT_flag_28 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 28) * (df_perm_treat$Genotype != "KO")
    F_MTX_KO_flag_28 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 28) * (df_perm_treat$Genotype == "KO")
    F_MTX_WT_flag_28 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 28) * (df_perm_treat$Genotype != "KO")
    M_Saline_KO_flag_28 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 28) * (df_perm_treat$Genotype == "KO")
    M_Saline_WT_flag_28 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 28) * (df_perm_treat$Genotype != "KO")
    F_Saline_KO_flag_28 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 28) * (df_perm_treat$Genotype == "KO")
    F_Saline_WT_flag_28 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 28) * (df_perm_treat$Genotype != "KO")
    
    M_MTX_KO_flag_42 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 42) * (df_perm_treat$Genotype == "KO")
    M_MTX_WT_flag_42 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 42) * (df_perm_treat$Genotype != "KO")
    F_MTX_KO_flag_42 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 42) * (df_perm_treat$Genotype == "KO")
    F_MTX_WT_flag_42 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 42) * (df_perm_treat$Genotype != "KO")
    M_Saline_KO_flag_42 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 42) * (df_perm_treat$Genotype == "KO")
    M_Saline_WT_flag_42 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 42) * (df_perm_treat$Genotype != "KO")
    F_Saline_KO_flag_42 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 42) * (df_perm_treat$Genotype == "KO")
    F_Saline_WT_flag_42 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 42) * (df_perm_treat$Genotype != "KO")
    
    M_MTX_KO_flag_63 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 63) * (df_perm_treat$Genotype == "KO")
    M_MTX_WT_flag_63 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 63) * (df_perm_treat$Genotype != "KO")
    F_MTX_KO_flag_63 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 63) * (df_perm_treat$Genotype == "KO")
    F_MTX_WT_flag_63 <- (df_perm_treat$Treatment == "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 63) * (df_perm_treat$Genotype != "KO")
    M_Saline_KO_flag_63 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 63) * (df_perm_treat$Genotype == "KO")
    M_Saline_WT_flag_63 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Male") * (df_perm_treat$Age == 63) * (df_perm_treat$Genotype != "KO")
    F_Saline_KO_flag_63 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 63) * (df_perm_treat$Genotype == "KO")
    F_Saline_WT_flag_63 <- (df_perm_treat$Treatment != "MTX") * (df_perm_treat$Sex == "Female") * (df_perm_treat$Age == 63) * (df_perm_treat$Genotype != "KO")
    
    ## LOOP THROUGH ALL STRUCTURES ##
    for (j in 1:length(structure_list)) {
      cstruct <- as.character(structure_list[j])
      cform <- as.formula( paste0("`",cstruct,"`"," ~ ",formula_rhs_2) )  
      clm_perm <- lmer(cform, data = df_perm_treat)
      slm <- summary(clm_perm)    
      df_dummy_treat[j, allcols] <- as.vector(slm$coefficients)
      df_dummy_treat$Permutation <- p
    }
    
    ## COMBINE ALL RESULTS INTO ONE BIG DATA FRAME ##
    perm_treat_df <- rbind(perm_treat_df, df_dummy_treat)
  }
  
  ## SAVE RESULTS ##
  save(perm_treat_df, file = "PERMUTATION_RESULTS_TREATMENT_1000.RData") 
  
} 

####################################################################################################
## COMPARE TREATMENT EFFECTS TO PERMUTATION DISTRIBUTIONS USING EMPIRICAL P-VALUES               ##
## NORMALIZE ESTIMATES BY SEX AND AGE, COMPUTE MEDIANS, AND CALCULATE TWO-TAILED P-VALUES        ##
## DATA SOURCES: PERMUTATION_TREATMENT_1000.RData & ALLRES_2.RData                               ##
####################################################################################################

perm_treat_df <- perm_treat_df %>%
  mutate(
    ## WT FEMALE ##
    F_WT_MTX_28_p = (`F_MTX_WT_flag_28.Estimate` / `age_factor28.Estimate`) * 100,
    
    ## WT MALE ##
    Male_WT_Saline_Estimate_28_p = `age_factor28:sex_flag.Estimate` + `age_factor28.Estimate`,
    M_WT_MTX_28_p = (`M_MTX_WT_flag_28.Estimate` / Male_WT_Saline_Estimate_28_p) * 100,
    
    ## KO FEMALE ##
    F_KO_Saline_28_Estimate_p = `F_Saline_KO_flag_28.Estimate` + `age_factor28.Estimate`,
    F_KO_MTX_28_p = (`F_MTX_KO_flag_28.Estimate` / F_KO_Saline_28_Estimate_p) * 100,
    
    # KO MALE ##
    M_KO_Saline_28_Estimate_p = `M_Saline_KO_flag_28.Estimate` + Male_WT_Saline_Estimate_28_p,
    M_KO_MTX_28_p = (`M_MTX_KO_flag_28.Estimate` / M_KO_Saline_28_Estimate_p) * 100
  )

## COMPUTE MEDIANS FOR PERMUTATIONS ##
perm_medians <- perm_treat_df %>%
  group_by(Permutation) %>%  
  summarise(
    Median_F_WT_MTX = median(F_WT_MTX_28_p),
    Median_F_KO_MTX = median(F_KO_MTX_28_p),
    Median_M_WT_MTX = median(M_WT_MTX_28_p),
    Median_M_KO_MTX = median(M_KO_MTX_28_p)
  )

observed_F_WT_MTX_median <- median(F_WT_MTX_28)
observed_F_KO_MTX_median <- median(F_KO_MTX_28)
observed_M_WT_MTX_median <- median(M_WT_MTX_28)
observed_M_KO_MTX_median <- median(M_KO_MTX_28)

## STATISTICAL TESTING USING TWO-TAILED EMPERICAL P-VALUE ##

## FEMALES ##

## WT ##
N <- length(perm_medians$Median_F_WT_MTX)
j <- sum(observed_F_WT_MTX_median > perm_medians$Median_F_WT_MTX)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median WT-MTX Effect across all structures (FEMALE, P28): ",as.character(observed_F_WT_MTX_median),"p=",as.character(mypval)) )  #0.078  

## KO ##
N <- length(perm_medians$Median_F_KO_MTX)
j <- sum(observed_F_KO_MTX_median > perm_medians$Median_F_KO_MTX)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median KO-MTX Effect across all structures (FEMALE, P28): ",as.character(observed_F_KO_MTX_median),"p=",as.character(mypval)) )  #0  

## MALES ##

## WT ##
N <- length(perm_medians$Median_M_WT_MTX)
j <- sum(observed_M_WT_MTX_median > perm_medians$Median_M_WT_MTX)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median WT-MTX Effect across all structures (MALE, P28): ",as.character(observed_M_WT_MTX_median),"p=",as.character(mypval)) )  #0.006  

## KO ##
N <- length(perm_medians$Median_M_KO_MTX)
j <- sum(observed_M_KO_MTX_median > perm_medians$Median_M_KO_MTX)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median KO-MTX Effect across all structures (MALE, P28): ",as.character(observed_M_KO_MTX_median),"p=",as.character(mypval)) )  #0  


####################################################################################################
## STATISTICAL COMPARISON OF KO VS WT UNDER MTX AT P28 USING WILCOXON RANK-SUM TEST              ##
## GROUPS COMPARED: FEMALE KO vs WT AND MALE KO vs WT                                            ##
## OUTPUT: TWO-SIDED P-VALUES INDICATING SIGNIFICANT DISTRIBUTION DIFFERENCES                    ##
####################################################################################################

# Compare F KO vs F WT  
wilcox_female_MTX_comparison <- wilcox.test(F_WT_MTX_28, F_KO_MTX_28, 
                                            alternative = "two.sided")
print(wilcox_female_MTX_comparison)
# Significant #
# p-value = 0.0009795 #
# THE TWO HISTOGRAMS ARE SIGNIFICANTLY DIFFERENCE THAN ONE ANOTHER #

# Compare M KO vs M WT
wilcox_male_MTX_comparison <- wilcox.test(M_WT_MTX_28, M_KO_MTX_28, 
                                          alternative = "two.sided")
print(wilcox_male_MTX_comparison)
# Significant #
# p-value = 6.635e-05 #
# THE TWO HISTOGRAMS ARE SIGNIFICANTLY DIFFERENCE THAN ONE ANOTHER #

#################################################################################################################
## STRUCTUREWISE GENOTYPE × TREATMENT INTERACTION EFFECTS AT P28                                               ##
## FOR SUPPLEMENTARY FIGURE 1: SIGNIFICANT EFFECTS ONLY (P < 0.05), SPLIT BY SEX                               ##
## USED TO CREATE BRAIN MAP FOR SUPPLEMENTARY FIGURE                                                           ##
#################################################################################################################

## FOR INTERPRETABILITY, COMPUTE PERCENT CHANGE IN STRUCTURE VOLS AT EACH TIME POINT ##

## FEMALE ##
F_28_pc <- allres$`F_treatment_flag_28:F_genotype_flag_28.Estimate` / allres$age_factor28.Estimate * 100

## MALE ##
Male_WT_Saline_Estimate_28 <- allres$'age_factor28:sex_flag.Estimate' + allres$'age_factor28.Estimate'
M_28_pc <- allres$`M_treatment_flag_28:M_genotype_flag_28.Estimate` / Male_WT_Saline_Estimate_28 * 100

## SUMMARIZE # OF AFFECTED STRUCTURES ##
print( paste0("IL6KO MTX VOL CHANGE (FEMALE, P28, FDR<0.1): ",as.character(sum( allres$`F_treatment_flag_28:F_genotype_flag_28.Pr(>|t|)` < 0.05 ))," structures") ) #19
print( paste0("IL6KO MTX VOL CHANGE (MALE, P28, FDR<0.1): ",as.character(sum( allres$`M_treatment_flag_28:M_genotype_flag_28.Pr(>|t|)` < 0.05 ))," structures") ) #5

######################################################################################################
## GENERATE HISTOGRAM TO GENERATE GENOTYPE × TREATMENT INTERACTION EFFECTS                          ##
## OUTPUT: PERMUTATION_RESULTS_INT_1000.RData FOR STATISTICAL COMPARISON (FIGURE 8B)                ##
######################################################################################################

# NORMALIZE ESTIMATES #

## FEMALE INTERACTION EFFECT P28 ##
F_INT_28_pc <- (allres$`F_treatment_flag_28:F_genotype_flag_28.Estimate` / allres$`age_factor28.Estimate`) * 100

## MALE INTERACTION EFFECT P28 ##
Male_WT_Saline_Estimate_28 <- allres$`age_factor28:sex_flag.Estimate` + allres$`age_factor28.Estimate`
M_INT_28_pc <- (allres$`M_treatment_flag_28:M_genotype_flag_28.Estimate` / Male_WT_Saline_Estimate_28) * 100

# CREATE DATA FRAMES FOR EACH INTERACTION ESTIMATE (ADDING COLUMNS FOR SEX AND INTERACTION)
m_data_INT <- data.frame(Estimate = M_INT_28_pc, Sex = "Male", Interaction = "M_treatment_flag_28:M_genotype_flag_28")
f_data_INT <- data.frame(Estimate = F_INT_28_pc, Sex = "Female", Interaction = "F_treatment_flag_28:F_genotype_flag_28")

combined_data <- rbind(f_data_INT, m_data_INT)

## RELABEL FACET TITLES ##
combined_data$Interaction <- factor(
  combined_data$Interaction,
  levels = c(
    "F_treatment_flag_28:F_genotype_flag_28",
    "M_treatment_flag_28:M_genotype_flag_28"
  ),
  labels = c(
    "Female IL6KO MTX Interaction",
    "Male IL6KO MTX Interaction"
  )
)

# DEFINE X-AXIS LIMITS #
x_min <- min(combined_data$Estimate, na.rm = TRUE)
x_max <- max(combined_data$Estimate, na.rm = TRUE)

# COLOURS #
sex_colors <- c("Male" = "#9370DB", "Female" = "#D68CA6")

## MAKE FEMALE APPEAR FIRST FOR PLOTTING ##
combined_data$Sex <- factor(combined_data$Sex, levels = c("Female", "Male"))

## PLOT ##
cplt<- ggplot(combined_data, aes(x = Estimate, fill = Sex)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.8, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_fill_manual(values = sex_colors) +
  labs(
    x = "% Volume Change (Relative to Interaction of Genotype and Treatment Differences)",
    y = "Number of Brain Structures"
  ) +
  xlim(x_min, x_max) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"),
    axis.title.x = element_text(size = 14, family = "Arial"),
    axis.title.y = element_text(size = 14, family = "Arial"),
    axis.text = element_text(size = 14, family = "Arial"),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold", family = "Arial"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~ Interaction, ncol = 2)

cairo_pdf("Figure8B.pdf",width=10.0,height=5.5)
print(cplt) 
dev.off()

######################################################################################################
## GENERATE 1000 PERMUTATIONS TO TEST GENOTYPE × TREATMENT INTERACTION EFFECTS                     ##
## THIS SCRIPT SHUFFLES GENOTYPE WITHIN METHOTREXATE-TREATED MICE (BY SEX) AND RE-RUNS MODELS      ##
## OUTPUT: PERMUTATION_RESULTS_INT_1000.RData FOR STATISTICAL COMPARISON (FIGURE 8)                ##
######################################################################################################

if (file.exists("PERMUTATION_RESULTS_INT_1000.RData")){
  load("PERMUTATION_RESULTS_INT_1000.RData") #load perm_results_df if already saved to file
} else {
  ## SET UP EMPTY RESULTS DATA FRAME ##
  allres$Permutation <- NA
  
  ## PERMUTATION PARAMETERS ##
  n_permutation <- 1000
  perm_results_df_2 <- data.frame()  
  
  ## FILTER FOR ONLY MTX MICE ## 
  FINAL_DATA_FRAME_MTX <- df %>%
    filter(Treatment == "MTX")
  
  ## FILTER FOR ONLY SALINE MICE ## 
  FINAL_DATA_FRAME_Saline <- df %>%
    filter(Treatment == "Saline")
  
  ## EXTRACTING ID & GENOTYPE ## 
  df_ID_gen_2 <- distinct(df, ID, Genotype, Sex, .keep_all = FALSE)
  
  ## REMOVE GENOTYPE FROM MAIN DATA FRAME ##
  df_no_gen_2 <- select(df, -Genotype)
  
  ## LOOP THROUGH PERMUTATIONS BELOW 
  set.seed(123)
  for (p in 1:n_permutation) {
    
    df_dummy_int <- allres
    
    df_ID_F_int <- subset(df_ID_gen_2, Sex=="Female")
    df_ID_M_int <- subset(df_ID_gen_2, Sex=="Male")
    
    ## SHUFFLE GENOTYOPE BY ID (ENSURE SAME ID GETS SAME GENOTYPE ACROSS TIMEPOINTS) 
    df_shuffle_F_int <- transform(df_ID_F_int, Genotype = sample(df_ID_F_int$Genotype))
    df_shuffle_M_int <- transform(df_ID_M_int, Genotype = sample(df_ID_M_int$Genotype))
    df_shuffle_int <- rbind(df_shuffle_F_int, df_shuffle_M_int)
    df_shuffle_int <- df_shuffle_int[c("ID", "Genotype")]
    
    ## MERGE BACK SHUFFLED TREATMENT ##
    df_perm_1 <- left_join(df_no_gen_2, df_shuffle_int, by = "ID")
    df_perm_2 <- rbind(df_perm_1, FINAL_DATA_FRAME_Saline)
    
    ## REDEFINE FLAGS WITH SHUFFLED TREATMENT ## 
    df_perm_2$age_factor <- as.character(df_perm_2$Age)
    df_perm_2$sex_flag <- (df_perm_2$Sex == "Male") * 1
    
    df_perm_2$M_treatment_flag_28 <- (df_perm_2$Treatment == "MTX") * (df_perm_2$Sex == "Male") * (df_perm_2$Age == 28)
    df_perm_2$F_treatment_flag_28 <- (df_perm_2$Treatment == "MTX") * (df_perm_2$Sex == "Female") * (df_perm_2$Age == 28)
    df_perm_2$M_treatment_flag_42 <- (df_perm_2$Treatment == "MTX") * (df_perm_2$Sex == "Male") * (df_perm_2$Age == 42)
    df_perm_2$F_treatment_flag_42 <- (df_perm_2$Treatment == "MTX") * (df_perm_2$Sex == "Female") * (df_perm_2$Age == 42)
    df_perm_2$M_treatment_flag_63 <- (df_perm_2$Treatment == "MTX") * (df_perm_2$Sex == "Male") * (df_perm_2$Age == 63)
    df_perm_2$F_treatment_flag_63 <- (df_perm_2$Treatment == "MTX") * (df_perm_2$Sex == "Female") * (df_perm_2$Age == 63)
    
    df_perm_2$M_genotype_flag_14 <- (df_perm_2$Genotype == "KO") * (df_perm_2$Sex == "Male") * (df_perm_2$Age == 14)
    df_perm_2$F_genotype_flag_14 <- (df_perm_2$Genotype == "KO") * (df_perm_2$Sex == "Female") * (df_perm_2$Age == 14)
    df_perm_2$M_genotype_flag_28 <- (df_perm_2$Genotype == "KO") * (df_perm_2$Sex == "Male") * (df_perm_2$Age == 28)
    df_perm_2$F_genotype_flag_28 <- (df_perm_2$Genotype == "KO") * (df_perm_2$Sex == "Female") * (df_perm_2$Age == 28)
    df_perm_2$M_genotype_flag_42 <- (df_perm_2$Genotype == "KO") * (df_perm_2$Sex == "Male") * (df_perm_2$Age == 42)
    df_perm_2$F_genotype_flag_42 <- (df_perm_2$Genotype == "KO") * (df_perm_2$Sex == "Female") * (df_perm_2$Age == 42)
    df_perm_2$M_genotype_flag_63 <- (df_perm_2$Genotype == "KO") * (df_perm_2$Sex == "Male") * (df_perm_2$Age == 63)
    df_perm_2$F_genotype_flag_63 <- (df_perm_2$Genotype == "KO") * (df_perm_2$Sex == "Female") * (df_perm_2$Age == 63)
    
    ## LOOP THROUGH ALL STRUCTURES ##
    for (j in 1:length(structure_list)) {
      cstruct <- as.character(structure_list[j])
      cform <- as.formula( paste0("`",cstruct,"`"," ~ ",formula_rhs) )  
      clm_perm <- lmer(cform, data = df_perm_2)
      slm <- summary(clm_perm)    
      df_dummy_int[j, allcols] <- as.vector(slm$coefficients)
      df_dummy_int$Permutation <- p
    }
    
    ## COMBINE ALL RESULTS INTO ONE BIG DATA FRAME ##
    perm_results_df_2 <- rbind(perm_results_df_2, df_dummy_int)
  }
  
  ## SAVE RESULTS ##
  save(perm_results_df_2, file = "PERMUTATION_RESULTS_INT_1000.RData") 
  
}  

####################################################################################################
## COMPARE OBSERVED GENOTYPE × TREATMENT INTERACTIONS TO PERMUTATION DISTRIBUTIONS AT P28        ##
## NORMALIZE INTERACTION ESTIMATES BY SEX, COMPUTE MEDIANS, AND CALCULATE EMPIRICAL P-VALUES     ##
## DATA SOURCES: PERMUTATION_RESULTS_INT_1000.RData & ALLRES.RData                               ##
####################################################################################################
perm_results_df_2 <- perm_results_df_2 %>%
  mutate(
    ## FEMALE ##
    F_INT_28_pc_p = (`F_treatment_flag_28:F_genotype_flag_28.Estimate` / `age_factor28.Estimate`) * 100,
    
    ## MALE ##
    Male_WT_Saline_Estimate_28_p = `age_factor28:sex_flag.Estimate` + `age_factor28.Estimate`,
    M_INT_28_pc_p = (`M_treatment_flag_28:M_genotype_flag_28.Estimate` / Male_WT_Saline_Estimate_28_p) * 100
  )

## COMPUTE MEDIANS FOR EACH PERMUTATION ##
perm_medians <- perm_results_df_2 %>%
  group_by(Permutation) %>%  
  summarise(
    Median_M_INT_28_pc_p = median(M_INT_28_pc_p, na.rm = TRUE),
    Median_F_INT_28_pc_p = median(F_INT_28_pc_p, na.rm = TRUE)
  )

observed_F_median <- median(F_INT_28_pc)
observed_M_median <- median(M_INT_28_pc)

## TWO-TAILED EMPERICAL P-VALUE ##

N <- length(perm_medians$Median_F_INT_28_pc_p)
j <- sum(observed_F_median> perm_medians$Median_F_INT_28_pc_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median KO-MTX Interaction Effect across all structures (FEMALE, P28): ",as.character(observed_F_median),"p=",as.character(mypval)) )  #0.454  


N <- length(perm_medians$Median_M_INT_28_pc_p)
j <- sum(observed_M_median> perm_medians$Median_M_INT_28_pc_p)
mypval <- ifelse(j>N/2,2*(N-j)/N,2*j/N)
print( paste0("Median KO-MTX Interaction Effect across all structures (MALE, P28): ",as.character(observed_M_median),"p=",as.character(mypval)) )  #0.512  


####################################################################################################
## STATISTICAL COMPARISON OF GENOTYPE × TREATMENT INTERACTION EFFECTS AT P28 USING WILCOX TEST   ##
## GROUPS COMPARED: FEMALE INTERACTION vs MALE INTERACTION (KO × MTX)                            ##
## OUTPUT: TWO-SIDED P-VALUE INDICATING HIGHLY SIGNIFICANT DIFFERENCE BETWEEN DISTRIBUTIONS      ##
####################################################################################################


wilcox_interaction_comparison <- wilcox.test(F_INT_28_pc, M_INT_28_pc, 
                                             alternative = "two.sided")
print(wilcox_interaction_comparison)
# Significant #
# p-value = 1.585e-15 #


####################################################################################################
## MOUSE WEIGHT TRAJECTORY PLOT                                                                   ##
####################################################################################################

## LOAD DATA ##
FINAL_WEIGHT <- read.csv("FINAL_WEIGHT.csv")

## CHANGE TO NUMERIC ##
FINAL_WEIGHT$Age <- as.numeric(as.character(FINAL_WEIGHT$Age))
FINAL_WEIGHT$Weight <- as.numeric(as.character(FINAL_WEIGHT$Weight))

## CALCULATE MEAN CONFIDENCE INTERVALS ##
weightc <- FINAL_WEIGHT %>%
  group_by(Group, Sex, Age) %>%
  summarise(
    mean = mean(Weight, na.rm = TRUE),
    se = sd(Weight, na.rm = TRUE) / sqrt(n()),
    ci_lower = mean - qt(0.975, df = n() - 1) * se,
    ci_upper = mean + qt(0.975, df = n() - 1) * se,
    N = n(),
    .groups = 'drop'
  )

## RE-ORDER GROUP LEVELS FOR LEGEND ##
weightc$Group <- factor(weightc$Group, levels = c("WT+Saline", "WT+MTX", "KO+Saline", "KO+MTX"))

## PLOTTING ##
cplt <- ggplot(weightc, aes(x = Age, y = mean, group = Group, color = Group)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Group),
              alpha = 0.1, linetype = "dotted", size = 0.5) +
  labs(title = "Mouse Weight",
       x = "Age (days)",
       y = "Weight (g)",
       color = "Group") +
  scale_color_manual(values = c("WT+Saline" = "blue",
                                "WT+MTX" = "red",
                                "KO+Saline" = "black",
                                "KO+MTX" = "grey")) +
  scale_fill_manual(values = c("WT+Saline" = "blue",
                               "WT+MTX" = "red",
                               "KO+Saline" = "black",
                               "KO+MTX" = "grey")) +
  facet_wrap(~Sex) +
  scale_x_continuous(trans = 'log10',
                     breaks = c(13, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 35, 41, 49, 56, 62)) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    panel.spacing = unit(2, "lines")
  )

## SAVE TO PDF ##
cairo_pdf("Figure2.pdf", width = 10.0, height = 5.5)
print(cplt)
dev.off()

####################################################################################################
## MOUSE WEIGHT STATISTICAL CALCULATIONS                                                          ##
####################################################################################################


####################################################################################################
## FEMALE WT SALINE BASELINE                                                                      ##
####################################################################################################

## SET BASELINE LEVELS ##
FINAL_WEIGHT$Group <- relevel(factor(FINAL_WEIGHT$Group), ref = "WT+Saline")
FINAL_WEIGHT$Age <- as.factor(FINAL_WEIGHT$Age)
FINAL_WEIGHT$Sex <- factor(FINAL_WEIGHT$Sex, levels = c("Female", "Male"))

## FIT THE MODEL ##
Linear_mixed_model <- lmer(
  Weight ~ Age * Group * Sex + (1 | ID) + (1 | Home_Cage),
  data = FINAL_WEIGHT
)

coefs <- summary(Linear_mixed_model)$coefficients
sig_effects_WT_F <- coefs[grep("^Age[0-9]+:GroupWT\\+MTX$", rownames(coefs)) , ]
sig_effects_WT_F <- sig_effects_WT_F[sig_effects_WT_F[, "Pr(>|t|)"] < 0.05, ]
print(sig_effects_WT_F[, "Pr(>|t|)"])
## WT MTX FEMALE MICE ARE SIGNIFICANTLY SMALLER THAN WT SALINE FEMALE MICE FROM P22 TO P27 ##

####################################################################################################
## FEMALE KO SALINE BASELINE                                                                      ##
####################################################################################################

## SET BASELINE LEVELS ##
FINAL_WEIGHT$Group <- relevel(factor(FINAL_WEIGHT$Group), ref = "KO+Saline")
FINAL_WEIGHT$Sex <- factor(FINAL_WEIGHT$Sex, levels = c("Female", "Male"))

## FIT THE MODEL ##
Linear_mixed_model <- lmer(
  Weight ~ Age * Group * Sex + (1 | ID) + (1 | Home_Cage),
  data = FINAL_WEIGHT
)

coefs <- summary(Linear_mixed_model)$coefficients
sig_effects_KO_F <- coefs[grep("^Age[0-9]+:GroupKO\\+MTX$", rownames(coefs)) , ]
sig_effects_KO_F <- sig_effects_KO_F[sig_effects_KO_F[, "Pr(>|t|)"] < 0.05, ]
print(sig_effects_KO_F[, "Pr(>|t|)"])
## KO MTX FEMALE MICE ARE SIGNIFICANTLY SMALLER THAN KO SALINE FEMALE MICE FROM P22 TO P35 ##

####################################################################################################
## MALE KO SALINE BASELINE                                                                      ##
####################################################################################################

## SET BASELINE LEVELS ##
FINAL_WEIGHT$Group <- relevel(factor(FINAL_WEIGHT$Group), ref = "KO+Saline")
FINAL_WEIGHT$Sex <- factor(FINAL_WEIGHT$Sex, levels = c("Male", "Female"))

## FIT THE MODEL ##
Linear_mixed_model <- lmer(
  Weight ~ Age * Group * Sex + (1 | ID) + (1 | Home_Cage),
  data = FINAL_WEIGHT
)

coefs <- summary(Linear_mixed_model)$coefficients
sig_effects_KO_M <- coefs[grep("^Age[0-9]+:GroupKO\\+MTX$", rownames(coefs)) , ]
sig_effects_KO_M <- sig_effects_KO_M[sig_effects_KO_M[, "Pr(>|t|)"] < 0.05, ]
print(sig_effects_KO_M[, "Pr(>|t|)"])
## KO MTX MALE MICE ARE SIGNIFICANTLY SMALLER THAN KO SALINE MALE MICE FROM P22 TO P49 ##

####################################################################################################
## MALE WT SALINE BASELINE                                                                        ##
####################################################################################################

## SET BASELINE LEVELS ##
FINAL_WEIGHT$Group <- relevel(factor(FINAL_WEIGHT$Group), ref = "WT+Saline")
FINAL_WEIGHT$Sex <- factor(FINAL_WEIGHT$Sex, levels = c("Male", "Female"))

## FIT THE MODEL ##
Linear_mixed_model <- lmer(
  Weight ~ Age * Group * Sex + (1 | ID) + (1 | Home_Cage),
  data = FINAL_WEIGHT
)

coefs <- summary(Linear_mixed_model)$coefficients
sig_effects_WT_M <- coefs[grep("^Age[0-9]+:GroupWT\\+MTX$", rownames(coefs)) , ]
sig_effects_WT_M <- sig_effects_WT_M[sig_effects_WT_M[, "Pr(>|t|)"] < 0.05, ]
print(sig_effects_WT_M[, "Pr(>|t|)"])
## WT MTX MALE MICE ARE SIGNIFICANTLY SMALLER THAN WT SALINE MALE MICE FROM P22 TO P35 ##


