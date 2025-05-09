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


