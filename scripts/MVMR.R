library(tidyverse)   
library(cowplot)
library(TwoSampleMR) 
library(LDlinkR)    
library(RadialMR)     
library(phenoscanner) 
library(gt)
library(remotes)
library(MVMR)         
library(RMVMR)


### EXPOSURE DATA 1: SLEEP APNEA ----

## Import exposure GWAS SumStats
sa_path = "data/SA_reformatted.tsv.gz"
sa_exp_dat <- read_tsv(sa_path)
head(sa_exp_dat)

## Format data to TwoSampleMR format
exposure1 <- sa_exp_dat %>%
  mutate(TRAIT= "sa") %>%
  format_data(.,
              type = "exposure",
              snps = NULL,
              header = TRUE,
              phenotype_col = "TRAIT",
              snp_col = "SNP",
              beta_col = "BETA",
              se_col = "SE",
              eaf_col = "FRQ",
              effect_allele_col = "A2",
              other_allele_col = "A1",
              pval_col = "P",
              samplesize_col = "N",
              chr_col = "CHR",
              pos_col = "BP",
              z_col = "Z",
              log_pval = FALSE
  ) %>%
  as_tibble()


## Perform LD clumping on SNP data, filter SNPs to make it run faster
exposure1_clump <- exposure1 %>% 
  rename(rsid = SNP, pval = pval.exposure) %>%
  ieugwasr::ld_clump(
    dat = .,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p = 0.99,
    bfile = 'resources/EUR_All_Chr',
    plink_bin = 'plink_win64_20230116/plink'
  )

exposure1_dat <- exposure1_clump %>%
  rename(SNP = rsid, pval.exposure= pval)


### EXPOSURE DATA 2 - BMI ----

## Reformatted BMI dataset 

# Define column types for summary statistics
coltypes = cols(
  ID = col_character(),
  CHROM = col_double(),
  POS = col_double(),
  REF = col_character(),
  ALT = col_character(),
  AF = col_double(),
  TRAIT = col_character(),
  BETA = col_double(),
  SE = col_double(),
  Z = col_double(),
  P = col_double(),
  N = col_double(),
  OR = col_double(),
  OR_L95 = col_double(),
  OR_U95 = col_double(),
  DIR = col_character(),
  G1000_ID = col_character(),
  G1000_VARIANT = col_character(),
  DBSNP_ID = col_character(),
  DBSNP_VARIANT = col_character(),
  OLD_ID = col_character(),
  OLD_VARIANT = col_character()
)

## Import exposure GWAS SumStats
bmi_path = "data/Yengo2018bmi.chrall.CPRA_b37.tsv.gz"
bmi_exp_dat <- read_tsv(bmi_path, comment = "##",  col_types = coltypes, 
                        col_select = c(DBSNP_ID, CHROM, POS, REF, ALT, AF, BETA, SE, Z, P, N, TRAIT))
head(bmi_exp_dat)

## Format data to TwoSampleMR format
exposure2 <- bmi_exp_dat %>%
  mutate(TRAIT= "BMI") %>%
  format_data(.,
              type = "exposure",
              snps = NULL,
              header = TRUE,
              phenotype_col = "TRAIT",
              snp_col = "DBSNP_ID",
              beta_col = "BETA",
              se_col = "SE",
              eaf_col = "AF",
              effect_allele_col = "ALT",
              other_allele_col = "REF",
              pval_col = "P",
              samplesize_col = "N",
              z_col = "Z",
              chr_col = "CHROM",
              pos_col = "POS",
              log_pval = FALSE
  ) %>%
  as_tibble()

# Perform clumping to obtain independent genome-wide significant variants for 
# each exposure. This step involves identifying the SNPs that are independently
# associated with each exposure and are significant at the genome-wide level

# Clump BMI
exposure2_clump <- exposure2 %>% 
  rename(rsid = SNP, pval = pval.exposure) %>%
  ieugwasr::ld_clump(
    dat = .,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p = 0.99,
    bfile = 'resources/EUR_All_Chr',
    plink_bin = 'plink_win64_20230116/plink'
  )

exposure2_dat <- exposure2_clump %>%
  rename(SNP = rsid, pval.exposure= pval)


# Combine exposures

# Make list of distinct SNPs from the exposures independent GWS snps 
mvmr_snps <- bind_rows(
  select(exposure1_dat, SNP), 
  select(exposure2_dat, SNP), 
) %>%
  distinct(SNP) %>%
  pull(SNP)

# extract combined SNP list from exposure datasets 
comb_exp <- bind_rows(
  exposure1 %>% filter(SNP %in% mvmr_snps),
  exposure2 %>% filter(SNP %in% mvmr_snps),
)


### OUTCOME DATA: CHD ----

## Import outcome GWAS SumStats - clean ones
chd_path = "data/CHD_reformatted.tsv.gz"
chd_out_dat <- read_tsv(chd_path)
head(chd_out_dat)

## Format data to TwoSampleMR format
outcome <- chd_out_dat %>%
  format_data(.,
              type = "outcome",
              snps = NULL,
              header = TRUE,
              snp_col = "SNP",
              beta_col = "BETA",
              se_col = "SE",
              eaf_col = "FRQ",
              effect_allele_col = "A2",
              other_allele_col = "A1",
              pval_col = "P",
              samplesize_col = "N",
              chr_col = "CHR",
              z_col = "Z",
              pos_col = "BP",
              log_pval = FALSE
  ) %>%
  as_tibble()


# Combine the exposure SNP lists and extract all the SNPs from each exposure. 
# This step combines the lists of SNPs for each exposure and extracts all the 
# SNPs that are present in each exposure.
## Unique SNP ids across exposures
exp_snp_list <- full_join(
  select(exposure1_dat, SNP, exposure), 
  select(exposure2_dat, SNP, exposure), 
  by = 'SNP'
)

exp1_unq_snps <- filter(exp_snp_list, !is.na(exposure.x) & is.na(exposure.y)) %>% nrow()
exp2_unq_snps <- filter(exp_snp_list, !is.na(exposure.y) & is.na(exposure.x)) %>% nrow()
both_unq_snps <- filter(exp_snp_list, !is.na(exposure.y) & !is.na(exposure.x)) %>% nrow()


# Combined exposures SNPs

## exposure and outcome snps lists 
comb_exp_out <- bind_rows(
  exposure1 %>% filter(SNP %in% mvmr_snps) %>% 
    rename(pval = pval.exposure, trait = exposure, chr = chr.exposure, pos = pos.exposure), 
  exposure2 %>% filter(SNP %in% mvmr_snps) %>% 
    rename(pval = pval.exposure, trait = exposure, chr = chr.exposure, pos = pos.exposure), 
  outcome %>% filter(SNP %in% mvmr_snps) %>% 
    rename(pval = pval.outcome, trait = outcome, chr = chr.outcome, pos = pos.outcome)
)

## How many SNPs are GWS across the exposures 
comb_exp_wide <- comb_exp_out %>%
  arrange(chr, pos) %>%
  select(SNP, trait, pval) %>%
  pivot_wider(names_from = trait, values_from = pval) %>%
  mutate(pval.exposure = pmin(sa, BMI, na.rm = T))

comb_exp_wide %>% count(sa < 5e-8, BMI < 5e-8)  %>%
  gt() 


## Identify Proxy Variants - sleep apnea
exposure1_clump2 <- semi_join(exposure1, exposure2_dat, by = "SNP")
exp4_snps_wo <- anti_join(exposure2_dat, exposure1, by = "SNP")

LDproxy_batch(exp4_snps_wo$SNP, 
              pop = "CEU",             # Match population ancestries
              r2d = "r2", 
              token = 'a6deee62cc4a', 
              append = TRUE,           # We appended the results of each LDlink query to a single file
              genome_build = "grch37") # Select genome build based on summary stats

system("mv combined_query_snp_list_grch37.txt data/exposure1_proxy_snps.txt")

## MVMR Munge proxies function
mvmr_munge_proxies <- function(LDLink_file, exposure){
  LDLink_file_path <- LDLink_file
  proxy_snps <- read_tsv(LDLink_file_path, skip = 1, col_names = F) %>%
    dplyr::rename(id = X1, func = X2, proxy_snp = X3, coord = X4, alleles = X5, maf = X6, 
                  distance = X7, dprime = X8, rsq = X9, correlated_alleles = X10, FORGEdb = X11, RegulomeDB = X12) %>%
    separate(coord, c('chr', 'pos'), sep = ":") %>%
    mutate(snp = ifelse(id == 1, proxy_snp, NA), 
           chr = str_replace(chr, 'chr', ""), 
           chr = as.numeric(chr), 
           pos = as.numeric(pos)) %>%
    fill(snp, .direction = 'down') %>%
    relocate(snp, .before = proxy_snp) %>%
    dplyr::select(-id, -func, -FORGEdb, -RegulomeDB) %>%
    filter(rsq >= 0.8)
  
  # Rename columns 
  outcome <- exposure %>% 
    magrittr::set_colnames(str_replace_all(colnames(.), "exposure", "outcome"))
  
  # Munge proxy snp and outcome data
  proxy_outcome <- left_join(
    proxy_snps, outcome, by = c("proxy_snp" = "SNP")
  ) %>%
    separate(correlated_alleles, c("target_a1.outcome", "proxy_a1.outcome", 
                                   "target_a2.outcome", "proxy_a2.outcome"), sep = ",|=") %>%
    filter(!is.na(chr.outcome)) %>%
    arrange(snp, -rsq, abs(distance)) %>%
    group_by(snp) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    mutate(
      proxy.outcome = TRUE,
      target_snp.outcome = snp,
      proxy_snp.outcome = proxy_snp, 
    ) %>% 
    mutate(
      new_effect_allele.outcome = case_when(
        proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a1.outcome,
        proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a2.outcome,
        TRUE ~ NA_character_
      ), 
      new_other_allele.outcome = case_when(
        proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a2.outcome,
        proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a1.outcome,
        TRUE ~ NA_character_
      ), 
      effect_allele.outcome = new_effect_allele.outcome, 
      other_allele.outcome = new_other_allele.outcome
    ) %>%
    dplyr::select(-proxy_snp, -chr, -pos, -alleles, -maf, -distance, -rsq, -dprime,  
                  -new_effect_allele.outcome, -new_other_allele.outcome) %>%
    relocate(target_a1.outcome, proxy_a1.outcome, target_a2.outcome, proxy_a2.outcome, .after = proxy_snp.outcome) %>%
    dplyr::rename(SNP = snp) %>%
    relocate(SNP, .after = samplesize.outcome)
  
  # Merge outcome and proxy outcomes
  out <- proxy_outcome %>%
    arrange(chr.outcome, pos.outcome) %>%
    magrittr::set_colnames(str_replace_all(colnames(.), "outcome", "exposure")) 
  
  out
  
}

# Munge proxy snp file
munge_proxy_exposure1_dat <- mvmr_munge_proxies("data/exposure1_proxy_snps.txt", exposure1)


## Identify Proxy Variants - BMI
outcome_clump1 <- semi_join(exposure2, exposure1_dat, by = "SNP")
exp3_snps_wo <- anti_join(exposure1_dat, exposure2, by = "SNP")

LDproxy_batch(exp3_snps_wo$SNP, 
              pop = "CEU",             # Match population ancestries
              r2d = "r2", 
              token = 'a6deee62cc4a', 
              append = TRUE,           # We appended the results of each LDlink query to a single file
              genome_build = "grch37") # Select genome build based on summary stats
system("mv combined_query_snp_list_grch37.txt data/exposure2_proxy_snps.txt") 

# Munge proxy snp file
munge_proxy_exposure2_dat <- mvmr_munge_proxies("data/exposure2_proxy_snps.txt", exposure2)


## Combine proxy SNPs and clump
comb_exp_proxy <- bind_rows(
  exposure1 %>% filter(SNP %in% mvmr_snps), 
  munge_proxy_exposure1_dat %>%
    select(-c(X13, proxy.exposure, target_snp.exposure, proxy_snp.exposure, 
              target_a1.exposure, proxy_a1.exposure, target_a2.exposure, proxy_a2.exposure)),
  exposure2 %>% filter(SNP %in% mvmr_snps),
  munge_proxy_exposure2_dat %>%
    select(-c(X13, proxy.exposure, target_snp.exposure, proxy_snp.exposure, 
              target_a1.exposure, proxy_a1.exposure, target_a2.exposure, proxy_a2.exposure)),
) 

exposure1_chd_bmi_proxy <- comb_exp_proxy %>%
  select(SNP, pval.exposure, exposure) %>%
  pivot_wider(names_from = exposure, values_from = pval.exposure) %>%
  mutate(pval.exposure = pmin(sa, BMI, na.rm = T)
  )

## Perform LD clumping on the combined SNP list to retain indepdent SNPs. 
## This step is performed to reduce the risk of spurious results arising due to multi-collinearity by including correlated SNPs. 
## LD clumping identifies lead SNPs in each locus and discards the remaining SNPs that are in linkage disequilibrium with them.

exposure1_chd_bmi_proxy_clumb <- exposure1_chd_bmi_proxy %>% 
  clump_data(.,
             clump_kb = 10000,
             clump_r2 = 0.001,
             clump_p1 = 1,
             clump_p2 = 1,
             pop = "EUR"
  )

exposure1_chd_bmi_proxy_clumb_independent <- comb_exp_proxy %>% filter(SNP %in% 
                                                                         exposure1_chd_bmi_proxy_clumb$SNP)


write_tsv(exposure1_chd_bmi_proxy_clumb_independent, 
          "data/combined_exposures.tsv.gz")

## Can start here
exposure1_chd_bmi_proxy_clumb_independent <- read_tsv("data/combined_exposures.tsv.gz")


## Extract Exposure SNPs from Outcome ----

# extract exposure SNPs present in outcome
outcome_clump <- semi_join(
  outcome, exposure1_chd_bmi_proxy_clumb_independent, by = "SNP"
)

# Exposure SNPs not present in outcome
outcome_wo_snps <- anti_join(
  exposure1_chd_bmi_proxy_clumb_independent, outcome, by = "SNP"
) %>%
  distinct(SNP)


## Proxy
# Use LDLinkR to identify proxy snps
LDproxy_batch(outcome_wo_snps, 
              pop = "CEU", 
              r2d = "r2", 
              token = 'd897616e6fa6', 
              append = TRUE,
              genome_build = "grch37")

system("mv combined_query_snp_list_grch37.txt data/mvmr_proxy_snps.txt")

# Munge proxy snp file

# Function for munging LDlink output
munge_proxies <- function(LDLink_file, outcome, outcome_clump){
  LDLink_file_path <- LDLink_file
  proxy_snps <- read_tsv(LDLink_file_path, skip = 1, col_names = F) %>%
    dplyr::rename(id = X1, func = X2, proxy_snp = X3, coord = X4, alleles = X5, maf = X6, 
                  distance = X7, dprime = X8, rsq = X9, correlated_alleles = X10, FORGEdb = X11, RegulomeDB = X12) %>%
    separate(coord, c('chr', 'pos'), sep = ":") %>%
    mutate(snp = ifelse(id == 1, proxy_snp, NA), 
           chr = str_replace(chr, 'chr', ""), 
           chr = as.numeric(chr), 
           pos = as.numeric(pos)) %>%
    fill(snp, .direction = 'down') %>%
    relocate(snp, .before = proxy_snp) %>%
    dplyr::select(-id, -func, -FORGEdb, -RegulomeDB) %>%
    filter(rsq >= 0.8)
  
  # Munge proxy snp and outcome data
  proxy_outcome <- left_join(
    proxy_snps, outcome, by = c("proxy_snp" = "SNP")
  ) %>%
    separate(correlated_alleles, c("target_a1.outcome", "proxy_a1.outcome", 
                                   "target_a2.outcome", "proxy_a2.outcome"), sep = ",|=") %>%
    filter(!is.na(chr.outcome)) %>%
    arrange(snp, -rsq, abs(distance)) %>%
    group_by(snp) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    mutate(
      proxy.outcome = TRUE,
      target_snp.outcome = snp,
      proxy_snp.outcome = proxy_snp, 
    ) %>% 
    mutate(
      new_effect_allele.outcome = case_when(
        proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a1.outcome,
        proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a2.outcome,
        TRUE ~ NA_character_
      ), 
      new_other_allele.outcome = case_when(
        proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a2.outcome,
        proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a1.outcome,
        TRUE ~ NA_character_
      ), 
      effect_allele.outcome = new_effect_allele.outcome, 
      other_allele.outcome = new_other_allele.outcome
    ) %>%
    dplyr::select(-proxy_snp, -chr, -pos, -alleles, -maf, -distance, -rsq, -dprime,  
                  -new_effect_allele.outcome, -new_other_allele.outcome) %>%
    relocate(target_a1.outcome, proxy_a1.outcome, target_a2.outcome, proxy_a2.outcome, .after = proxy_snp.outcome) %>%
    dplyr::rename(SNP = snp) %>%
    relocate(SNP, .after = samplesize.outcome)
  
  # Merge outcome and proxy outcomes
  outcome_dat <- bind_rows(
    outcome_clump, proxy_outcome
  ) %>% 
    arrange(chr.outcome, pos.outcome)
  
  outcome_dat
}

outcome_dat <- munge_proxies("data/mvmr_proxy_snps.txt", outcome, outcome_clump) 


## HARMONIZE ----

## Remove APOE region
comb_exposure_dat <- exposure1_chd_bmi_proxy_clumb_independent 
# %>%
#   filter(!(chr.exposure == 19 & between(pos.exposure, 44912079, 45912079)))

## Univariate MR 
mrdat <- harmonise_data(comb_exposure_dat, outcome_dat, action = 2) %>%
  as_tibble() %>%
  filter(pval.exposure < 5e-8)


#Let's remove the SNP that was identified as outlier in the univariate analysis
comb_exposure_dat2 <- comb_exposure_dat[comb_exposure_dat$SNP != "rs10878269",]
outcome_dat2 <- outcome_dat[outcome_dat$SNP != "rs10878269",]

mvdat <- mv_harmonise_data(comb_exposure_dat2, outcome_dat2)
mvdat

## Export datasets
write_csv(mrdat, 'data/harmonized_mvmr_uni.csv')
write_rds(mvdat, 'data/harmonized_mvmr.rds')


#Read in harmonized data (can start here)
mrdat <- read_csv("data/harmonized_mvmr_uni.csv")
mvdat <- read_rds('data/harmonized_mvmr.rds')

mvmr_res <- mv_multiple(mvdat, plots = F)
mvmr_res


###############################################################
###############################################################

## MVMR PACKAGE

# Format for MVMR
mvmr_dat <- bind_cols(
  rownames(mvdat$exposure_beta),
  mvdat$exposure_beta, 
  mvdat$exposure_se, 
  mvdat$outcome_beta, 
  mvdat$outcome_se
) %>%
  as.matrix()

F.data <- format_mvmr(BXGs = mvmr_dat[,c(2,3)],
                      BYG = mvmr_dat[,6],
                      seBXGs = mvmr_dat[,c(4,5)],
                      seBYG = mvmr_dat[,7],
                      RSID = mvmr_dat[,1])

# Results for MVMR
res_mvmr <- ivw_mvmr(r_input = F.data)

mvmr_out <- res_mvmr %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>% 
  janitor::clean_names() %>%
  rename(exposure = rowname, b = estimate, se = std_error, p = pr_t) %>%
  mutate(
    exposure = ifelse(exposure == 'exposure1', 'sa', 'obesity'), 
    outcome = 'CHD',
    method = 'IVW-MVMR'
  ) %>%
  relocate(outcome, method, .after = exposure)


## MVMR diagnostics
sres <- strength_mvmr(r_input = F.data, gencov = 0)

pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)

## Radial MVMR
F.data_rmvmr <- format_rmvmr(BXGs = mvmr_dat[,c(2,3)],
                             BYG = mvmr_dat[,6],
                             seBXGs = mvmr_dat[,c(4,5)],
                             seBYG = mvmr_dat[,7],
                             RSID = mvmr_dat[,1])


res_rmvrm <- ivw_rmvmr(F.data_rmvmr)

sres_rmvrm <- strength_rmvmr(F.data_rmvmr, gencov = 0)

rmvrm_p <- plot_rmvmr(F.data_rmvmr, res_rmvrm) %>% 
  labs(title = "Radial MVMR")
rmvrm_p