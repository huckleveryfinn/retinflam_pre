#
# retpsy_do.R
#
# created on Wed Jun 14 12:03:36 2023
# Finn Rabe, <finn dot rabe at bli dot uzh dot ch>
#-----------------------------------------------------------------------


## Load data and create timestamp for makefile
# setwd("retinflam/src")
source("retinflam_load.R")
file.create("../output/R/retinflam_do.Rout")

# run python script that create filtered dataframes and figures
# use_python("/usr/local/bin/python")
# source_python("./retinflam_do1.py")

df_norm <- read.csv("../output/data/df_normtest.csv")
df_norm <- subset(df_norm, select = c("Subfield", "Statistic", "p"))
df_norm$p <- format_p(df_norm$p, stars = )


# function to round p vals
p_round <- function(x, n = 2) {
    max(abs(round(x, n)), abs(signif(x, 1))) * sign(x)
}

# function converts p-values to asterisk
signif.num <- function(x) {
    symnum(x,
        corr = FALSE, na = FALSE, legend = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
        symbols = c("***", "**", "*", ".", "")
    )
}

# load filtered subfields and compute residuals
df_filt <- read.csv("../output/data/df_mlm.csv")
prsvar <- "PRSSZ"
ret_list <- c(
    "Inner_Inferior_left", "Inner_Inferior_right",
    "Outer_Inferior_left", "Outer_Inferior_right",
    "Inner_Nasal_left", "Inner_Nasal_right",
    "Outer_Nasal_left", "Outer_Nasal_right",
    "Inner_Superior_left", "Inner_Superior_right",
    "Outer_Superior_left", "Outer_Superior_right",
    "Inner_Temporal_left", "Inner_Temporal_right",
    "Outer_Temporal_left", "Outer_Temporal_right",
    "Central_left", "Central_right"
)
df_ret_res <- setNames(data.frame(matrix(ncol = length(ret_list), nrow = dim(df_filt)[1])), ret_list)
cov <- "Age + Age_squared + Sex + Smoking_status + BMI + Eye_disorders + Genotype_array + Townsend_index + Genetic.PC1+Genetic.PC2+ Genetic.PC3+ Genetic.PC4+ Genetic.PC5+Genetic.PC6+ Genetic.PC7+ Genetic.PC8+ Genetic.PC9+ Genetic.PC10"

# Loop over the new variable names
for (var_name in ret_list) {
    # Get eye label
    lat_lst <- as.list(strsplit(var_name, "_")[[1]])[-1]
    lat <- lat_lst[length(lat_lst)]
    octqc <- paste0("OCT_quality_", lat)
    maccen <- paste0("Macula_centered_", lat)

    # Create the formula string
    pc_str <- paste(var_name, "~", cov, "+", octqc, "+", maccen)
    # Convert the string to a formula
    pc_obj <- as.formula(pc_str)
    # Run the regression
    m_pc <- rlm(pc_obj, data = df_filt, psi = psi.huber)
    # Obtain residuals
    res <- scale(m_pc$resid)
    df_ret_res[var_name] <- res
}
write.csv(df_ret_res, "../output/data/retinflam_rlm_resid.csv", row.names = FALSE)

# run analysis with resid regressed out
# source_python("./retinflam_do.py")

# load filtered overall macula data frame and compute rlm
df_macula <- read.csv("../output/data/df_mlm_macula.csv")
paired_test_result <- t.test(df_macula$Macula_right, df_macula$Macula_left, paired = TRUE)
ttest_macula_stat <- round(paired_test_result$statistic[[1]], digits = 2)
ttest_macula_pval <- paired_test_result$p.value

## Associations between overall macular thickness and PRSSZ while controlling for confounds
m_macula_left <- rlm(Macula_left ~ PRSSZ + Age + Age_squared + Sex + Smoking_status + BMI + Eye_disorders + OCT_quality_left + Macula_centered_left + Genotype_array + Townsend_index + Genetic.PC1 + Genetic.PC2 + Genetic.PC3 + Genetic.PC4 + Genetic.PC5 + Genetic.PC6 + Genetic.PC7 + Genetic.PC8 + Genetic.PC9 + Genetic.PC10,
    data = df_macula, psi = psi.huber
)
pc_dd <- data.frame(summary(m_macula_left)$coefficients)
lmacula_b <- round(pc_dd$Value[2], digits = 2)
lmacula_se <- round(pc_dd$Std..Error[2], digits = 5)
lmacula_conflow <- round(confint.default(object = m_macula_left, parm = prsvar, level = 0.95)[1], digits = 2)
lmacula_confhigh <- round(confint.default(object = m_macula_left, parm = prsvar, level = 0.95)[2], digits = 2)
fpc <- f.robftest(m_macula_left, var = prsvar)
fpcFval <- round(fpc$statistic, digits = 2)
# lmacula_F <- map(fpcFval, 1)[[1]]
lmacula_F <- fpcFval[[1]]
lmacula_p <- round(fpc$p.value, digits = 5)

m_macula_right <- rlm(Macula_right ~ PRSSZ + Age + Age_squared + Sex + Smoking_status + BMI + Eye_disorders + OCT_quality_right + Macula_centered_right + Genotype_array + Townsend_index + Genetic.PC1 + Genetic.PC2 + Genetic.PC3 + Genetic.PC4 + Genetic.PC5 + Genetic.PC6 + Genetic.PC7 + Genetic.PC8 + Genetic.PC9 + Genetic.PC10,
    data = df_macula, psi = psi.huber
)
pc_dd <- data.frame(summary(m_macula_right)$coefficients)
rmacula_b <- round(pc_dd$Value[2], digits = 2)
rmacula_se <- round(pc_dd$Std..Error[2], digits = 5)
rmacula_conflow <- round(confint.default(object = m_macula_right, parm = prsvar, level = 0.95)[1], digits = 2)
rmacula_confhigh <- round(confint.default(object = m_macula_right, parm = prsvar, level = 0.95)[2], digits = 2)
fpc <- f.robftest(m_macula_right, var = prsvar)
fpcFval <- round(fpc$statistic, digits = 2)
# rmacula_F <- map(fpcFval, 1)[[1]]
rmacula_F <- fpcFval[[1]]
rmacula_p <- round(fpc$p.value, digits = 5)
tab_model(m_macula_left, m_macula_right, file = "../output/figures/TableS1_PCRLM_Macula.html", show.fstat = TRUE)
macularfwer <- p.adjust(c(lmacula_p, rmacula_p), method = "holm")


### Load filtered/standardized subfield's and pc data
df_filt <- read.csv("../output/data/df_mlm_pc.csv")

## Associations between polygenic risk for schizophrenia and PC1 and PC2
m_pc1 <- rlm(PC1 ~ PRSSZ,
    data = df_filt, psi = psi.huber
)
pc1_dd <- data.frame(summary(m_pc1)$coefficients)
pc1_b <- round(pc1_dd$Value[2], digits = 2)
pc1_se <- round(pc1_dd$Std..Error[2], digits = 5)
pc1_conflow <- round(confint.default(object = m_macula_left, parm = prsvar, level = 0.95)[1], digits = 2)
pc1_confhigh <- round(confint.default(object = m_macula_left, parm = prsvar, level = 0.95)[2], digits = 2)
fpc <- f.robftest(m_pc1, var = prsvar)
fpcFval <- round(fpc$statistic, digits = 2)
# pc1_F <- map(fpcFval, 1)[[1]]
pc1_F <- fpcFval[[1]]
pc1_p <- round(fpc$p.value, digits = 4)

m_pc2 <- rlm(PC2 ~ PRSSZ,
    data = df_filt, psi = psi.huber
)
pc2_dd <- data.frame(summary(m_pc2)$coefficients)
pc2_b <- round(pc2_dd$Value[2], digits = 2)
pc2_se <- round(pc2_dd$Std..Error[2], digits = 5)
pc2_conflow <- round(confint.default(object = m_macula_left, parm = prsvar, level = 0.95)[1], digits = 3)
pc2_confhigh <- round(confint.default(object = m_macula_left, parm = prsvar, level = 0.95)[2], digits = 3)
fpc <- f.robftest(m_pc2, var = prsvar)
fpcFval <- round(fpc$statistic, digits = 2)
# pc2_F <- map(fpcFval, 1)[[1]]
pc2_F <- fpcFval[[1]]
pc2_p <- round(fpc$p.value, digits = 4)
m_table <- tab_model(m_pc1, m_pc2, file = "../output/figures/Table1_PCRLM.html", show.fstat = TRUE)
#webshot("../output/figures/Table1_PCRLM.html", "../output/figures/Table1_PCRLM.png", vwidth = 580, vheight = 300,zoom = 1)
pcs_fwer <- p.adjust(c(pc1_p, pc2_p), method = "holm")


## Pathway-specific associations
df_prs <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("PRS", "PC", "Coef", "CI(low)", "CI(high)", "SE", "F", "p"))
prs_list <- c("PRSSZ", "PRSSZ_0_001", "PRSSZ_0_05", "PRSSZ_0_1", "PRSSZ_0_2", "PRSSZ_0_3", "PRSSZ_0_4", "PRSSZ_0_5", "PRSSZ_1") # ,"PRSSZ_bestfit")
prs_list_mod <- c("PRSSZ", "PRSSZ 0.001", "PRSSZ 0.05", "PRSSZ 0.1", "PRSSZ 0.2", "PRSSZ 0.3", "PRSSZ 0.4", "PRSSZ 0.5", "PRSSZ 1") # ,"PRSSZ_bestfit")
prs_list_mod <- rep(prs_list_mod, each = 2)
pc_list <- c("PC1", "PC2")
# cov <- "Age + Age_squared + Sex + Smoking_status + BMI + Eye_disorders + OCT_quality_left + OCT_quality_right + Genotype_array + Townsend_index + Genetic.PC1+Genetic.PC2+ Genetic.PC3+ Genetic.PC4+ Genetic.PC5+Genetic.PC6+ Genetic.PC7+ Genetic.PC8+ Genetic.PC9+ Genetic.PC10"
# Loop over the new variable names
i <- 0
pfdr_cols <- c()
for (var_name in prs_list) {
    pvals_pc <- c()
    for (pc in pc_list) {
        # Create the formula string
        pc_str <- paste(pc, "~", var_name)
        # Convert the string to a formula
        pc_obj <- as.formula(pc_str)
        # Run the regression
        m_pc <- rlm(pc_obj, data = df_filt, psi = psi.huber)
        n_model_left <- nobs(m_pc)

        # pc_intercept<- round(map(m_pc$coefficients["(Intercept)"], 1)[[1]],digits = 2)
        pc_intercept <- round(m_pc$coefficients["(Intercept)"][[1]], digits = 2)
        pc1bp <- bptest(m_pc)
        tab1_pcbppval <- round(map(pc1bp$p.value, 1)[[1]], digits = 2)
        # tab1_pcbpstat <- round(map(pc1bp$statistic, 1)[[1]],digits = 2)
        tab1_pcbpstat <- round(pc1bp$statistic[[1]], digits = 2)

        pc_dd <- data.frame(summary(m_pc)$coefficients)
        pc_b <- round(pc_dd$Value[2], digits = 2)
        pc_se <- round(pc_dd$Std..Error[2], digits = 5)
        lgc_conflow <- round(confint.default(object = m_pc, parm = var_name, level = 0.95)[1], digits = 3)
        lgc_confhigh <- round(confint.default(object = m_pc, parm = var_name, level = 0.95)[2], digits = 3)

        fpc <- f.robftest(m_pc, var = var_name)
        fpcFval <- round(fpc$statistic, digits = 2)
        # fpcFval <- map(fpcFval, 1)[[1]]
        fpcFval <- fpcFval[[1]]
        fpcpval <- round(fpc$p.value, digits = 4)

        pvals_pc <- append(pvals_pc, fpcpval)
        # data <- unlist(list(var_name,pc,pc_b,lgc_conflow,lgc_confhigh,pc_se,tab1_pcbpstat,tab1_pcbppval,fpcFval,fpcpval),recursive = FALSE)
        data <- unlist(list(var_name, pc, pc_b, lgc_conflow, lgc_confhigh, pc_se, fpcFval, fpcpval), recursive = FALSE)
        df_it <- as.data.frame(t(data))
        new <- rep(i, ncol(df_it))
        df_prs[nrow(df_it) + i, ] <- df_it
        i <- i + 1
    }
    pc_fdr <- p.adjust(pvals_pc, method = "holm")
    pfdr_cols <- append(pfdr_cols, pc_fdr)
}
# p_ast <- signif.num(pfdr_cols)
# pfdr_cols_ast <- paste(pfdr_cols,p_ast,sep="")
df_prs$PRS <- prs_list_mod
df_prs$Coef <- as.numeric(df_prs$Coef)
df_prs$F <- as.numeric(df_prs$F)
df_prs["pFWER"] <- pfdr_cols # pfdr_cols_ast

df_pathprs <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("Pathway", "Gesea Code", "PC", "Coef", "CI(low)", "CI(high)", "SE", "BPstat", "BPp", "F", "p"))
pathprs_list <- c(
    "NEUROINFLAM_PRS", "ACUTEINFLAM_PRS", "CHROINFLAM_PRS", "TGFB_PRS", "WNT_PRS", "CATENIN_PRS", "DOPPOSREG_PRS",
    # "MITO_PRS",
    "ABNOVAS_PRS", "CORART_PRS"
)
pc_list <- c("PC1", "PC2")
cpathprs_code <- c("M24927", "M6557", "M15140", "M18933", "M25305", "M17761", "M24111", "M43559", "M36658")
cpathprs_code <- rep(cpathprs_code, each = length(pc_list))
# Loop over the new variable names
i <- 0
pfdr_cols <- c()
for (var_name in pathprs_list) {
    pvals_pc <- c()
    for (pc in pc_list) {
        # Create the formula string
        pc_str <- paste(pc, "~", var_name)
        # Convert the string to a formula
        pc_obj <- as.formula(pc_str)
        # Run the regression
        m_pc <- rlm(pc_obj, data = df_filt, psi = psi.huber)

        # pc_intercept<- round(map(m_pc$coefficients["(Intercept)"], 1)[[1]],digits = 2)
        pc_intercept <- round(m_pc$coefficients["(Intercept)"][[1]], digits = 2)
        pc1bp <- bptest(m_pc)
        tab1_pcbppval <- round(map(pc1bp$p.value, 1)[[1]], digits = 2)
        # tab1_pcbpstat <- round(map(pc1bp$statistic, 1)[[1]],digits = 2)
        tab1_pcbpstat <- round(pc1bp$statistic[[1]], digits = 2)

        pc_dd <- data.frame(summary(m_pc)$coefficients)
        pc_b <- round(pc_dd$Value[2], digits = 2)
        pc_se <- round(pc_dd$Std..Error[2], digits = 5)
        lgc_conflow <- round(confint.default(object = m_pc, parm = var_name, level = 0.95)[1], digits = 3)
        lgc_confhigh <- round(confint.default(object = m_pc, parm = var_name, level = 0.95)[2], digits = 3)

        fpc <- f.robftest(m_pc, var = var_name)
        fpcFval <- round(fpc$statistic, digits = 2)
        # fpcFval <- map(fpcFval, 1)[[1]]
        fpcFval <- fpcFval[[1]]
        fpcpvalnum <- as.numeric(fpc$p.value)
        fpcpval <- round(fpcpvalnum, digits = 5)

        pvals_pc <- append(pvals_pc, fpcpval)
        data <- unlist(list(var_name, cpathprs_code[i + 1], pc, pc_b, lgc_conflow, lgc_confhigh, pc_se, tab1_pcbpstat, tab1_pcbppval, fpcFval, fpcpval), recursive = FALSE)
        df_it <- as.data.frame(t(data))
        new <- rep(i, ncol(df_it))
        df_pathprs[nrow(df_it) + i, ] <- df_it
        i <- i + 1
    }
    pc_fdr <- p.adjust(pvals_pc, method = "holm")
    pfdr_cols <- append(pfdr_cols, pc_fdr)
}
# p_ast <- signif.num(pfdr_cols)
# pfdr_cols_ast <- paste(pfdr_cols,p_ast,sep="")
df_pathprs$Coef <- as.numeric(df_pathprs$Coef)
df_pathprs$p <- as.numeric(df_pathprs$p)
df_pathprs$F <- as.numeric(df_pathprs$F)
df_pathprs["pFWER"] <- pfdr_cols # pfdr_cols_ast
df_pathprs["Groups"] <- c(
    "Inflammation", "Inflammation", "Inflammation", "Inflammation",
    "Inflammation", "Inflammation", "Inflammation", "Inflammation",
    "Developmental gene expression", "Developmental gene expression", "Developmental gene expression",
    "Developmental gene expression", "Developmental gene expression", "Developmental gene expression",
    #' Mitochondrial regulation','Mitochondrial regulation',
    "Microvasculature", "Microvasculature",
    "Microvasculature", "Microvasculature"
)
df_pathprs_ordered <- df_pathprs[order(df_pathprs$Pathway, df_pathprs$PC), ]

## add competitive p vals
df_pthcomp <- read.csv("../output/data/pathway_comp.csv")
df_pathcomp_reord <- df_pthcomp %>% arrange(factor(Pathway, levels = unique(df_pathprs$Pathway)))
df_pathcomp_reord["GeseaPathwayCode"] <- cpathprs_code
df_pathcomp_reord["Groups"] <- c(
    "Inflammation", "Inflammation", "Inflammation", "Inflammation",
    "Inflammation", "Inflammation", "Inflammation", "Inflammation",
    "Developmental", "Developmental", "Developmental",
    "Developmental", "Developmental", "Developmental",
    "Microvasculature", "Microvasculature",
    "Microvasculature", "Microvasculature"
)
df_pathcomp_reord <- subset(df_pathcomp_reord, select = c("Pathway", "GeseaPathwayCode", "Num_SNP", "PC", "Estimate", "SE", "selfcontained.p", "competitive.p", "Groups"))
# correct for multiple comparisons
df_pathcomp_reord <- df_pathcomp_reord %>%
    group_by(Pathway) %>%
    mutate(selfcon.pFWER = if_else(PC %in% c("PC1", "PC2"), p.adjust(selfcontained.p, method = "holm"), NA_real_)) %>%
    mutate(comp.pFWER = if_else(PC %in% c("PC1", "PC2"), p.adjust(competitive.p, method = "holm"), NA_real_)) %>%
    ungroup()
df_pathcomp_reord <- as.data.frame(df_pathcomp_reord)
df_pathcomp_reord["CIhigh"] <- round(df_pathcomp_reord$Estimate + (1.96 * df_pathcomp_reord$SE), digits = 3)
df_pathcomp_reord["CIlow"] <- round(df_pathcomp_reord$Estimate - (1.96 * df_pathcomp_reord$SE), digits = 3)
df_pathcomp_ord <- df_pathcomp_reord[order(df_pathcomp_reord$Pathway, df_pathcomp_reord$PC), ]
df_pathcomp_ord$Pathway <- gsub("_", " ", df_pathcomp_ord$Pathway)


## Best fit pathway PRS
df_pathprs_bfit <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Pathway", "PC", "F", "p"))
pathprs_bfit_list <- c(
    "NEUROINFLAM_Bestfit_PRS", "ACUTEINFLAM_Bestfit_PRS", "CHROINFLAM_Bestfit_PRS", "TGFB_Bestfit_PRS", "WNT_Bestfit_PRS", "CATENIN_Bestfit_PRS", "DOPPOSREG_Bestfit_PRS",
    # "MITO_Bestfit_PRS",
    "ABNOVAS_Bestfit_PRS", "CORART_Bestfit_PRS"
)
pc_list <- c("PC1", "PC2")
# Loop over the new variable names
i <- 0
pfdr_cols <- c()
for (var_name in pathprs_bfit_list) {
    pvals_pc <- c()
    for (pc in pc_list) {
        # Create the formula string
        pc_str <- paste(pc, "~", var_name)
        # Convert the string to a formula
        pc_obj <- as.formula(pc_str)
        # Run the regression
        m_pc <- rlm(pc_obj, data = df_filt, psi = psi.huber)

        fpc <- f.robftest(m_pc, var = var_name)
        fpcFval <- round(fpc$statistic, digits = 2)
        # fpcFval <- map(fpcFval, 1)[[1]]
        fpcFval <- fpcFval[[1]]
        fpcpval <- round(fpc$p.value, digits = 5)

        pvals_pc <- append(pvals_pc, fpcpval)
        data <- unlist(list(var_name, pc, fpcFval, fpcpval), recursive = FALSE)
        df_it <- as.data.frame(t(data))
        new <- rep(i, ncol(df_it))
        df_pathprs_bfit[nrow(df_it) + i, ] <- df_it
        i <- i + 1
    }
    pc_fdr <- p.adjust(pvals_pc, method = "holm")
    pfdr_cols <- append(pfdr_cols, pc_fdr)
}
df_pathprs_bfit["pFWER"] <- pfdr_cols
# df_pathprs['FbestFit'] <- as.numeric(df_pathprs_bfit$F)
# df_pathprs['pFDRbestFit'] <- as.numeric(df_pathprs_bfit$pFDR)


### Mediation model neuroinflam_prs --CRP--> PC1
# Path C
model_direct <- rlm(PC1 ~ NEUROINFLAM_PRS,
    data = df_filt, psi = psi.huber
)
pc_dd <- data.frame(summary(model_direct)$coefficients)
mdirect_b <- round(pc_dd$Value[2], digits = 2)
mdirect_se <- round(pc_dd$Std..Error[2], digits = 5)
mdirect_conflow <- round(confint.default(object = model_direct, parm = "NEUROINFLAM_PRS", level = 0.95)[1], digits = 3)
mdirect_confhigh <- round(confint.default(object = model_direct, parm = "NEUROINFLAM_PRS", level = 0.95)[2], digits = 3)
fpc <- f.robftest(model_direct, var = "NEUROINFLAM_PRS")
fpcFval <- round(fpc$statistic, digits = 2)
# crp_F <- map(fpcFval, 1)[[1]]
mdirect_F <- fpcFval[[1]]
mdirect_p <- round(fpc$p.value, digits = 4)

# Path A
model_mediator <- rlm(CRP ~ NEUROINFLAM_PRS + Age + Age_squared + Sex + Smoking_status + BMI + Eye_disorders + OCT_quality_left + Macula_centered_left + Genotype_array + Townsend_index + Genetic.PC1 + Genetic.PC2 + Genetic.PC3 + Genetic.PC4 + Genetic.PC5 + Genetic.PC6 + Genetic.PC7 + Genetic.PC8 + Genetic.PC9 + Genetic.PC10,
    data = df_filt
)
pc_dd <- data.frame(summary(model_mediator)$coefficients)
crp_b <- round(pc_dd$Value[2], digits = 2)
crp_se <- round(pc_dd$Std..Error[2], digits = 5)
crp_conflow <- round(confint.default(object = model_mediator, parm = "NEUROINFLAM_PRS", level = 0.95)[1], digits = 3)
crp_confhigh <- round(confint.default(object = model_mediator, parm = "NEUROINFLAM_PRS", level = 0.95)[2], digits = 3)
fpc <- f.robftest(model_mediator, var = "NEUROINFLAM_PRS")
fpcFval <- round(fpc$statistic, digits = 2)
# crp_F <- map(fpcFval, 1)[[1]]
crp_F <- fpcFval[[1]]
crp_p <- round(fpc$p.value, digits = 4)

# Path B
model_outcome <- rlm(PC1 ~ NEUROINFLAM_PRS + CRP,
    data = df_filt, psi = psi.huber
)
pc_dd <- data.frame(summary(model_outcome)$coefficients)
moutcome_b <- round(pc_dd$Value[2], digits = 2)
moutcome_se <- round(pc_dd$Std..Error[2], digits = 5)
moutcome_conflow <- round(confint.default(object = model_outcome, parm = "NEUROINFLAM_PRS", level = 0.95)[1], digits = 3)
moutcome_confhigh <- round(confint.default(object = model_outcome, parm = "NEUROINFLAM_PRS", level = 0.95)[2], digits = 3)
fpc <- f.robftest(model_outcome, var = "NEUROINFLAM_PRS")
fpcFval <- round(fpc$statistic, digits = 2)
moutcome_F <- fpcFval[[1]]
moutcome_p <- round(fpc$p.value, digits = 4)

model_rev <- rlm(NEUROINFLAM_PRS ~ PC1 + CRP,
    data = df_filt, psi = psi.huber
)

# Path C'
a_coefficient <- coef(model_mediator)["NEUROINFLAM_PRS"]
b_coefficient <- coef(model_outcome)["CRP"]
indirect_effect <- round(a_coefficient * b_coefficient, digits = 3)
indirect_effect <- indirect_effect[[1]]
set.seed(123) # For reproducibility
boot_indirect <- function(data, indices) {
    d <- data[indices, ] # Resample the data
    med_model <- rlm(CRP ~ NEUROINFLAM_PRS, data = d, psi = psi.huber)
    out_model <- rlm(PC1 ~ NEUROINFLAM_PRS + CRP, data = d, psi = psi.huber)
    a <- coef(med_model)["NEUROINFLAM_PRS"]
    b <- coef(out_model)["CRP"]
    return(a * b)
}
results <- boot(data = df_filt, statistic = boot_indirect, R = 1000)
mediate_ci <- boot.ci(results, type = "perc") # Calculate percentile confidence intervals
medCI_low <- round(mediate_ci$percent[4], digits = 3)
medCI_high <- round(mediate_ci$percent[5], digits = 4)
mediate_p_value <- boot.pval(boot_res = results, type = "perc", theta_null = 0)


## Associations between individuals subfield's thickness and PRSSZ while controlling for confounds
df_ret <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), c("Macular subfield", "hm", "Coef", "CI(low)", "CI(high)", "SE", "BPstat", "BPp", "F", "p"))
# ret_list <- c("GC_IPL_left","GC_IPL_right","Macula_left","Macula_right", "RNFL_left","RNFL_right")
ret_list <- c(
    "Inner_Inferior_left", "Inner_Inferior_right",
    "Outer_Inferior_left", "Outer_Inferior_right",
    "Inner_Nasal_left", "Inner_Nasal_right",
    "Outer_Nasal_left", "Outer_Nasal_right",
    "Inner_Superior_left", "Inner_Superior_right",
    "Outer_Superior_left", "Outer_Superior_right",
    "Inner_Temporal_left", "Inner_Temporal_right",
    "Outer_Temporal_left", "Outer_Temporal_right",
    "Central_left", "Central_right"
)
ret_list_mod <- gsub("_", " ", ret_list)
cov <- "Age + Age_squared + Sex + Smoking_status + BMI + Eye_disorders + Genotype_array + Townsend_index + Genetic.PC1+Genetic.PC2+ Genetic.PC3+ Genetic.PC4+ Genetic.PC5+Genetic.PC6+ Genetic.PC7+ Genetic.PC8+ Genetic.PC9+ Genetic.PC10"
# Loop over the new variable names
i <- 0
prsvar <- "PRSSZ"
pvals_pc <- c()
for (var_name in ret_list) {
    # Get eye label
    lat_lst <- as.list(strsplit(var_name, "_")[[1]])[-1]
    lat <- lat_lst[length(lat_lst)]
    octqc <- paste0("OCT_quality_", lat)
    maccen <- paste0("Macula_centered_", lat)

    # Create the formula string
    pc_str <- paste(var_name, "~", prsvar, "+", cov, "+", octqc, "+", maccen)
    # Convert the string to a formula
    pc_obj <- as.formula(pc_str)
    # Run the regression
    m_pc <- rlm(pc_obj, data = df_filt, psi = psi.huber)

    # pc_intercept<- round(map(m_pc$coefficients["(Intercept)"], 1)[[1]],digits = 2)
    pc_intercept <- round(m_pc$coefficients["(Intercept)"][[1]], digits = 2)
    pc1bp <- bptest(m_pc)
    tab1_pcbppval <- format_p(map(pc1bp$p.value, 1)[[1]], stars = FALSE)
    # tab1_pcbpstat <- round(map(pc1bp$statistic, 1)[[1]],digits = 2)
    tab1_pcbpstat <- round(pc1bp$statistic[[1]], digits = 2)

    pc_dd <- data.frame(summary(m_pc)$coefficients)
    # pc_b <- round(pc_dd$Value[2],digits = 2)
    pc_b <- pc_dd$Value[2]
    pc_se <- round(pc_dd$Std..Error[2], digits = 5)
    lgc_conflow <- round(confint.default(object = m_pc, parm = prsvar, level = 0.95)[1], digits = 3)
    lgc_confhigh <- round(confint.default(object = m_pc, parm = prsvar, level = 0.95)[2], digits = 3)

    fpc <- f.robftest(m_pc, var = prsvar)
    fpcFval <- round(fpc$statistic, digits = 2)
    fpcFval <- map(fpcFval, 1)[[1]]
    fpcpval <- round(fpc$p.value, digits = 5)

    pvals_pc <- append(pvals_pc, fpcpval)
    data <- unlist(list(var_name, lat, pc_b, lgc_conflow, lgc_confhigh, pc_se, tab1_pcbpstat, tab1_pcbppval, fpcFval, fpcpval), recursive = FALSE)
    df_it <- as.data.frame(t(data))
    new <- rep(i, ncol(df_it))
    df_ret[nrow(df_it) + i, ] <- df_it
    i <- i + 1
}
pc_fdr <- p.adjust(df_ret$p, method = "holm")
# pfdr_cols_ast <- paste(pc_fdr,p_ast,sep="")
df_ret$Coef <- as.numeric(df_ret$Coef)
df_ret$F <- as.numeric(df_ret$F)
df_ret$p <- as.numeric(df_ret$p)
df_ret["Macular subfield"] <- ret_list_mod
df_ret["pFWER"] <- pc_fdr # pfdr_cols_ast
df_ret["Groups"] <- c(
    "Inner Inferior", "Inner Inferior",
    "Outer Inferior", "Outer Inferior",
    "Inner Nasal", "Inner Nasal",
    "Outer Nasal", "Outer Nasal",
    "Inner Superior", "Inner Superior",
    "Outer Superior", "Outer Superior",
    "Inner Temporal", "Inner Temporal",
    "Outer Temporal", "Outer Temporal",
    "Central", "Central"
)

# test if F-stat is significantly different between left and right eye
lh_coef <- df_ret[df_ret$hm == "left", "Coef"]
rh_coef <- df_ret[df_ret$hm == "right", "Coef"]
normt_lh <- shapiro.test(lh_coef)
normt_rh <- shapiro.test(rh_coef)
paired_test_result <- t.test(lh_coef, rh_coef, paired = TRUE)
ttest_subf_stat <- round(paired_test_result$statistic[[1]], digits = 2)
ttest_subf_pval <- paired_test_result$p.value
ttest_subf_df <- paired_test_result$df

allfvals <- df_ret$Coef
color_range <- colorRampPalette(c("#065535", "#d3ffce"))(length(allfvals))[rank(allfvals)]
df_ret_col <- df_ret
df_ret_col["Fcolor"] <- color_range

subf_lbl <- c("Inner Superior", "Outer Superior", "Inner Temporal", "Outer Temporal", "Inner Inferior", "Outer Inferior", "Inner Nasal", "Outer Nasal")
# subf_num <- c('2','6','5','9','4','8','3','7')
subf_num <- c("IS", "OS", "IT", "OT", "II", "OI", "IN", "ON")
pie_data <- data.frame(
    subfield = subf_lbl,
    position = c(2020, 2021, 2020, 2021, 2020, 2021, 2020, 2021),
    thickness = c(10, 10, 10, 10, 10, 10, 10, 10)
)

# Run for each eye seperately
Fleft <- df_ret_col[df_ret_col$hm == "left", ]
Fright <- df_ret_col[df_ret_col$hm == "right", ]

df_subf_right <- Fright %>% arrange(factor(Groups, levels = subf_lbl))
df_subf_right$Coef <- round(df_subf_right$Coef, digits = 3)
# subfvals <- df_subf_sort$mean_F
color_code <- df_subf_right$Fcolor
# color_code = colorRampPalette(c('#d3ffce', '#065535'))(length(subfvals))[rank(subfvals)]
# color_code = c('black','black','grey','grey','white', 'white','white','white','white')


## Subfield associations results fundus plots (right eye)
subfheatmap_right <- pie_data %>%
    ggplot(aes(x = position, y = thickness, fill = subfield)) +
    geom_col(
        position = "fill", width = 1,
        color = "white", show.legend = FALSE
    ) +
    coord_polar(theta = "y", start = 225 * pi / 180) +
    lims(x = c(2019, 2022)) +
    scale_fill_manual(
        values = color_code[-length(color_code)],
        breaks = subf_lbl
    ) +
    geom_point(aes(x = 2019, y = 0),
        size = 49,
        shape = 21, fill = color_code[9], colour = "white"
    ) +
    geom_text(
        label = "CS", x = 2019, y = 0, size = 12, fontface = 2,
        show.legend = FALSE
    ) +
    geom_textpath(
        position = position_fill(vjust = .5),
        angle = 90, alpha = 1,
        aes(color = "black", label = subf_num), color = "black",
        size = 12, fontface = 2, show.legend = FALSE
    ) +
    theme_void()

ggsave("../output/figures/ETDRS_heatmap_right.png",
    subfheatmap_right,
    bg = "transparent",
    width = 3042, height = 3042, dpi = 300, units = "px"
)

img_inset <- image_read("../output/figures/ETDRS_heatmap_right.png")
img_inset <- image_scale(img_inset, "70%x")
img <- image_read("../data/retina_bg/macula_2D_right.png")

img_with_inset <- img %>% image_composite(
    img_inset,
    operator = "Atop",
    offset = "-300-200",
    gravity = "Center"
)
image_write(img_with_inset, "../output/figures/msubf_heatmap_right.png")


## Subfield associations results fundus plots (left eye)
df_subf_left <- Fleft %>% arrange(factor(Groups, levels = subf_lbl))
df_subf_left$Coef <- round(df_subf_left$Coef, digits = 3)
color_code <- df_subf_left$Fcolor
# color_code = colorRampPalette(c('#d3ffce', '#065535'))(length(subfvals))[rank(subfvals)]

subfheatmap_left <- pie_data %>% ggplot(aes(x = position, y = thickness, fill = subfield)) +
    geom_col(position = "fill", width = 1, color = "white", show.legend = FALSE) +
    coord_polar(theta = "y", start = 225 * pi / 180) +
    lims(x = c(2019, 2022)) +
    scale_fill_manual(values = color_code[-length(color_code)], breaks = subf_lbl) +
    geom_point(aes(x = 2019, y = 0), size = 49, shape = 21, fill = color_code[9], colour = "white") +
    geom_text(label = "CS", x = 2019, y = 0, size = 12, fontface = 2, show.legend = FALSE) +
    geom_textpath(
        position = position_fill(vjust = .5), angle = 90, alpha = 1,
        aes(color = "black", label = subf_num), color = "black", size = 12, fontface = 2, show.legend = FALSE
    ) +
    theme_void()

ggsave("../output/figures/ETDRS_heatmap_left.png", subfheatmap_left, bg = "transparent", width = 3042, height = 3042, dpi = 300, units = "px")

img_inset <- image_read("../output/figures/ETDRS_heatmap_left.png")
img_inset <- image_scale(img_inset, "70%x")
img <- image_read("../data/retina_bg/macula_2D_left.png")

img_with_inset <- img %>% image_composite(
    img_inset,
    operator = "Atop",
    offset = "+300-200",
    gravity = "Center"
)
image_write(img_with_inset, "../output/figures/msubf_heatmap_left.png")



## Draw graph of mediation analysis results
mediation_graph <- grViz("digraph {
graph [layout = dot, rankdir = LR]
 node [shape = rectangle, fontname = Helvetica]
rec1 [label = '@@1']
rec2 [label = '@@2']
rec3 [label = '@@3']
 # edge definitions with the node IDs
rec1 -> rec2 [label= '@@4', fontname='Helvetica']
rec1 -> rec3 [label= '@@6\n @@7', fontname='Helvetica']
rec2 -> rec3 [label= '@@5', fontname='Helvetica']
}
[1]:  paste0('Neuroinflammatory-specific PRSSZ')
[2]:  paste0('CRP')
[3]:  paste0('PC1')
[4]:  paste0('coef (A) = ', as.character(round(a_coefficient[[1]],digits=3)), '**')
[5]:  paste0('coef (B) = ', as.character(round(b_coefficient[[1]],digits=3)), '***')
[6]:  paste0('coef (C) = ', as.character(mdirect_b), '***')
[7]:  paste0('coef (C`) = ', as.character(indirect_effect), '**')
")
mediation_graph %>%
    export_svg() %>%
    charToRaw() %>%
    rsvg_pdf("../output/figures/Mediation_Diagramme.pdf")


## Diagram of exclusion/inclusion numbers
results <- read.csv("../output/data/retinflam_results_m.csv")
# create participant diagramme
graph <- grViz("digraph {
graph [layout = dot, rankdir = TB]
 node [shape = rectangle, fontname = Helvetica]
rec1 [label = '@@1\n @@2']
rec2 [label = '@@3\n-@@4\n-@@5']
rec3 [label = '@@6\n @@7']
rec4 [label = '@@8\n @@9']
rec5 [label = '@@10\n @@11']
rec6 [label = '@@12\n-@@13\n-@@14']
rec7 [label = '@@15\n @@16']
 # edge definitions with the node IDs
rec1 -> rec2 -> rec3 -> rec4 -> rec6 -> rec7
rec3 -> rec5
}
[1]:  paste0('Individuals')
[2]:  paste0('(n = ', results$n_total[1], ')')
[3]:  paste0('Excluded (n = ', results$n_ethn[1]+results$n_qen_qc[1], ')')
[4]:  paste0('SNP QC/Sample QC (n = ', results$n_qen_qc[1], ')')
[5]:  paste0('Non British/Irish ancestry (n = ', results$n_ethn[1], ')')
[6]:  paste0('Individuals')
[7]:  paste0('(n = ', results$n_total[1]-results$n_ethn[1]-results$n_qen_qc[1], ')')
[8]:  paste0('Not Diagnosed')
[9]:  paste0('(n = ', results$n_total[1]-results$n_ethn[1]-results$n_qen_qc[1]-results$n_diag[1], ')')
[10]: paste0('Diagnosed ICD-10, F20-29')
[11]: paste0('(n = ', results$n_diag[1], ')')
[12]: paste0('Excluded (n = ', results$n_nondiag_nan+results$n_antipsy[1], ')')
[13]: paste0('Using antipsychotics (n = ', results$n_antipsy[1], ')')
[14]: paste0('Incomplete data (n = ', results$n_nondiag_nan, ')')
[15]: paste0('Individuals')
[16]: paste0('(n = ', n_model_left, ')')
")
graph %>%
    export_svg() %>%
    charToRaw() %>%
    rsvg_pdf("../output/figures/Diagramme1.pdf")


## Define all variables for markdown
diag1 <- "../output/figures/Diagramme1.pdf"
mediat_graph <- "../output/figures/Mediation_Diagramme.pdf"

fig1 <- "../data/retina_bg/Fig1_eye_anat.png"
fig2 <- "../output/figures/Fig1_overall_retinap_partialreg.png"
fig3a <- "../output/figures/msubf_heatmap_left.png"
fig3b <- "../output/figures/msubf_heatmap_right.png"
fig5 <- "../output/figures/Fig3Aa_subfields_CRP_partialreg.png"

tab1 <- "../output/figures/Table1_PCRLM.png"
tabs1 <- "../output/figures/Appx_Table1_subfields_sex_ratio.png"
tabs2 <- "../output/figures/Appx_Table2_subfields_vif.png"

figs1 <- "../output/figures/Fig1_subfields_corrmatrix.png"
figs2 <- "../output/figures/Appx_Fig1_subfields_regmatrix.png"
figs3 <- "../output/figures/Appx_Fig2A_subfields_pca.png"
figs4 <- "../output/figures/Appx_Fig2B_subfields_pca_loadingscores.png"

figs6 <- "../output/figures/Appx_FigureS1_CRP_logtransform.png"

myImage <- png::readPNG(fig1)
myImage <- grid::rasterGrob(myImage, interpolate = TRUE)

# convert vif sumary to table
df_vifl <- read.csv("../output/data/Appx_Table2_subfields_vif.csv")
df_vif <- dplyr::select(df_vifl, Phenotype, VIF)
df_vif$Phenotype <- gsub("_", " ", df_vif$Phenotype)

# collect mutiple plots into figure 1
p1a <- readPNG(fig1)
p1ap <- rasterGrob(p1a, interpolate = TRUE)
p1b <- readPNG(fig2)
p1bp <- rasterGrob(p1b, interpolate = TRUE)
p1c <- readPNG(fig3a)
p1d <- readPNG(fig3b)
p1cp <- rasterGrob(p1c, interpolate = TRUE)
p1dp <- rasterGrob(p1d, interpolate = TRUE)
p1cnew <- cowplot::plot_grid(p1cp, p1dp, nrow = 1)
p1bcp <- cowplot::plot_grid(p1bp, p1cnew,
    nrow = 1, labels = c("b", "c"),
    label_size = 25
)

fig1new <- cowplot::plot_grid(p1ap, p1bcp,
    nrow = 2, labels = c("a"),
    label_size = 25
)
fig1newp <- ggsave(filename = "../output/figures/retinflam_fig1.png", fig1new)
# plot(fig1new)

# collect multi pics for pca figure s2
p3c <- readPNG(figs3)
p3d <- readPNG(figs4)
p3cp <- rasterGrob(p3c, interpolate = TRUE)
p3dp <- rasterGrob(p3d, interpolate = TRUE)
# p3cnew <- cowplot::plot_grid(p3cp, p3dp, nrow=1)
p3bcp <- cowplot::plot_grid(p3cp, p3dp,
    nrow = 1, labels = c("a", "b"), label_y = 0.7,
    label_size = 25
)
figs3newp <- ggsave(filename = "../output/figures/retinflam_pca_figs3.png", p3bcp)
figs3_comb <- figs3newp


# collect sample stats and transform into table 1
columns <- c("Characteristic", "N", "Mean (SD)")
df_nsex <- data.frame(matrix(ncol = length(columns)))
colnames(df_nsex) <- columns
df_nsex[nrow(df_nsex) + 1, ] <- list("Female", results$female, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Male", results$male, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Eye disorder", results$n_eyedis, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Current smoker", results$curr_smoker, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Previous smoker", results$prev_smoker, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Non smoker", results$non_smoker, "-")
df_nsex[nrow(df_nsex) + 1, ] <- list("Age", n_model_left, results$mean_age)
df_nsex[nrow(df_nsex) + 1, ] <- list("BMI", n_model_left, results$mean_bmi)
df_nsex[nrow(df_nsex) + 1, ] <- list("Townsend Index", n_model_left, results$mean_townsend_index)
# append retinal phenotypes measures
df_samp_retphen <- read_csv("../output/data/retphen_samp.csv")
for (i in 1:nrow(df_samp_retphen)) {
    row <- df_samp_retphen[i, ]
    df_nsex[nrow(df_nsex) + 1, ] <- list(row$retinal_phenotype, n_model_left, row$thickness)
}
df_nsex <- df_nsex[-(1), ]
