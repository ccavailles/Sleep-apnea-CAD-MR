library(readr)
library(tidyverse)
library(TwoSampleMR)
library(LDlinkR)
library(gt)


### EXPOSURE DATA ----

## Import exposure GWAS SumStats
chd_path = "data/CHD_reformatted.tsv.gz"
chd_exp_dat <- read_tsv(chd_path)
head(chd_exp_dat)

## Format data to TwoSampleMR format
exposure <- chd_exp_dat %>%
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
              z_col = "Z",
              chr_col = "CHR",
              pos_col = "BP",
              log_pval = FALSE
  ) %>%
  as_tibble()

## Perform LD clumping on SNP data, filter SNPs to make it run faster
exposure_clump <- exposure %>% 
  rename(rsid = SNP, pval = pval.exposure) %>%
  ieugwasr::ld_clump(
    dat = .,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p = 0.99,
    bfile = 'resources/EUR_All_Chr',
    plink_bin = 'plink_win64_20230116/plink'
  )

exposure_dat <- exposure_clump %>%
  rename(SNP = rsid, pval.exposure= pval)



### OUTCOME DATA ----

## Import outcome GWAS SumStats
sa_path = "data/SA_reformatted.tsv.gz"
sa_exp_dat <- read_tsv(sa_path)
head(sa_exp_dat)

outcome <- sa_exp_dat %>%
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
              pos_col = "BP",
              log_pval = FALSE
  ) %>%
  as_tibble()


## Identify Proxy Variants

# extract exposure SNPs present in outcome
outcome_clump <- semi_join(outcome, exposure_dat, by = "SNP")

# Exposure SNPs not present in outcome
exp_snps_wo <- anti_join(exposure_dat, outcome, by = "SNP")


## Munge proxy snp file

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


LDproxy_batch(exp_snps_wo$SNP, 
              pop = "CEU",             # Match population ancestries
              r2d = "r2", 
              token = 'a6deee62cc4a', 
              append = TRUE,           # We appended the results of each LDlink query to a single file
              genome_build = "grch37") # Select genome build based on summary stats

outcome_dat <- munge_proxies("combined_query_snp_list_grch37.txt", outcome, outcome_clump)



### HARMONIZE EXPOSURE - OUTCOME DATASETS ----

mr_dat <- harmonise_data(exposure_dat, outcome_dat, action = 2) %>% as_tibble() %>%
  mutate(
    gws.outcome = ifelse(pval.outcome < 5e-8, TRUE, FALSE),
    mr_keep = ifelse(mr_keep == FALSE | gws.outcome == TRUE, FALSE, TRUE)
  ) %>%
  filter(pval.exposure < 5e-8)



### TWO SAMPLE MR ----

## MR

mr_res <- mr(mr_dat, method_list = c(
  "mr_egger_regression", "mr_weighted_median", "mr_ivw_fe","mr_weighted_mode"))
generate_odds_ratios(mr_res)

res_single <- mr_singlesnp(mr_dat, all_method = c("mr_ivw_fe", "mr_egger_regression", 
                                                  "mr_weighted_median", "mr_weighted_mode")) %>% as_tibble()

## Visualization 

# MR scatter plot

scatter_p <- mr_scatter_plot(mr_res, mr_dat) 
scatter_p

scatter_out_p <- scatter_p[[1]] + theme_bw() + 
  guides(color=guide_legend(ncol =1)) + 
  theme(
    text = element_text(size = 8), 
  )
scatter_out_p

# MR forrest plot 

forrest_p <- mr_forest_plot(res_single)
forrest_p[[1]]



### DIAGNOSTICS AND SENSITIVITY ANALYSES ----

## Pleiotropy 

res_pleio <- mr_pleiotropy_test(mr_dat)

res_pleio %>%
  select(-id.exposure, -id.outcome, -outcome, -exposure) %>%
  gt() %>%
  fmt_number(
    columns = c('egger_intercept', 'se')
  ) %>%
  fmt_number(
    columns = pval,
    rows = pval > 0.001,
    decimals = 3
  ) %>% 
  fmt_scientific(
    columns = pval,
    rows = pval <= 0.001,
    decimals = 1
  )


# Funnel plots

funnel_p <- mr_funnel_plot(res_single)
funnel_out_p <- funnel_p[[1]] + theme_bw() + 
  guides(color=guide_legend(ncol =1)) + 
  theme(
    text = element_text(size = 8), 
  )
funnel_out_p


## Heterogeneity 

res_het <- mr_heterogeneity(mr_dat, method_list = c("mr_egger_regression", "mr_ivw"))

res_het %>%
  select(-id.exposure, -id.outcome, -outcome, -exposure) %>%
  gt() %>%
  fmt_number(
    columns = Q
  ) %>%
  fmt_number(
    columns = Q_pval,
    rows = Q_pval > 0.001,
    decimals = 3
  ) %>% 
  fmt_scientific(
    columns = Q_pval,
    rows = Q_pval <= 0.001,
    decimals = 1
  )



## Outliers

# Leave-one-out
res_loo <- mr_leaveoneout(mr_dat, method = mr_ivw_fe) %>% as_tibble()

loo_p <- mr_leaveoneout_plot(res_loo)
loo_p[[1]]



### RADIAL-MR ----

library(RadialMR)

## Radial-MR

# Format data 
radial_dat <- mr_dat %>% filter(mr_keep == T) %>% TwoSampleMR::dat_to_RadialMR()
rmr_dat <- radial_dat$exposure.outcome

# Run Radial MR
bonff = 0.05/nrow(rmr_dat) # Bonferonni Correction

radial_ivw_res <- ivw_radial(rmr_dat, alpha = bonff) 

radial_egger_res <- egger_radial(rmr_dat, alpha = bonff) 


## Radial plots

ivw_radial_p <- plot_radial(radial_ivw_res, radial_scale = F, show_outliers = F)
egger_radial_p <- plot_radial(radial_egger_res, radial_scale = F, show_outliers = F)

cowplot::plot_grid(
  ivw_radial_p + coord_fixed(ratio=0.25) + theme(legend.position = 'bottom'), 
  egger_radial_p + theme(legend.position = 'bottom'), 
  align = 'h'
)


## Remove outliers

# Phenoscanner 

library(devtools)
library(phenoscanner)

phewas_dat <- phenoscanner(
  snpquery=radial_ivw_res$outliers$SNP, 
  proxies='EUR', 
  r2 = 0.8
) 

radial_ivw_res$outliers %>% 
  left_join(filter(phewas_dat$snps, snp == rsid), by = c("SNP" = "snp")) %>%
  select(SNP, Q_statistic, p.value, hg19_coordinates, a1, a2, hgnc) %>% 
  DT::datatable()


## Modify the mrkeep variable to flag variants as outliers from radial MR for removal
mr_dat_outlier <- mr_dat %>%
  left_join(radial_ivw_res$dat) %>%
  mutate(
    mr_keep = case_when(
      Outliers == "Outlier" ~ FALSE, 
      TRUE ~ mr_keep
    )
  )


## MR analysis w/o outliers

# Standard MR
mr_res <- mr(mr_dat, method_list = c(
  c("mr_ivw_fe", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
))
mr_res

## MR analysis 
mr_res_outlier <- mr(mr_dat_outlier, method_list = c("mr_ivw_fe", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
mr_res_outlier
generate_odds_ratios(mr_res_outlier)

# Heterogeneity statistics 
mr_het_outlier <- mr_heterogeneity(mr_dat_outlier, method_list = c("mr_egger_regression", "mr_ivw"))
mr_het_outlier

res_single <- mr_singlesnp(mr_dat_outlier, all_method = c("mr_ivw_fe", "mr_egger_regression", 
                                                          "mr_weighted_median", "mr_weighted_mode")) %>% as_tibble()
funnel_p <- mr_funnel_plot(res_single)
funnel_out_p <- funnel_p[[1]] + theme_bw() + 
  guides(color=guide_legend(ncol =1)) + 
  theme(
    text = element_text(size = 8), 
  )
funnel_out_p


# Horizontal pleitropy
mr_plt_outlier <- mr_pleiotropy_test(mr_dat_outlier)
mr_plt_outlier

# MR-PRESSO
library(MRPRESSO)
data_test <- as.data.frame(mr_dat_outlier)
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", 
          SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
          data = data_test, NbDistribution = 1000,  SignifThreshold = 0.05)

## F stat
f_stat = function(N, K, R){
  f = ((N-K-1) / K) * (R/(1-R))
  f
}

snp.pve <- function(eaf, beta, se, n){
  (2*eaf*(1 - eaf)*beta^2) / (2 * beta * eaf * (1-eaf) + se^2 * 2 * n * eaf * (1-eaf))
}

f_res <- mr_dat_outlier %>%
  group_by(exposure, outcome) %>%
  filter(mr_keep == TRUE) %>%
  select(SNP, exposure, outcome, effect_allele.exposure, eaf.exposure, beta.exposure, se.exposure) %>%
  mutate(
    samplesize.exposure = 1165690, 
    pve = snp.pve(eaf.exposure, beta.exposure, se.exposure, samplesize.exposure), 
    f = f_stat(samplesize.exposure, 1, pve),
    # f = abs(beta.exposure)^2 / se.exposure^2
  ) %>% 
  summarise(
    pve = sum(pve, na.rm = T), 
    k = n(), 
    samplesize = max(samplesize.exposure), 
    f = mean(f, na.rm = T),
  ) 
f_res


# Radial MR 
radial_dat_outlier <- mr_dat_outlier %>% filter(mr_keep == T) %>% dat_to_RadialMR()
radial_res_outlier <- ivw_radial(radial_dat_outlier$exposure.outcome, alpha = 0.05/nrow(radial_dat_outlier$exposure.outcome)) 

radial_res_outlier


bind_rows(
  mr_res %>% mutate(model = "w/ Outliers"), 
  mr_res_outlier %>% mutate(model = "w/o Outliers")) %>%
  select(model, method, nsnp, b, se, pval) %>%
  gt(
    groupname_col = "model"
  ) %>%
  fmt_number(
    columns = c("b", "se")
  ) %>%
  fmt_number(
    columns = pval,
    rows = pval > 0.001,
    decimals = 3
  ) %>% 
  fmt_scientific(
    columns = pval,
    rows = pval <= 0.001,
    decimals = 1
  ) 


# MR plots w/o outliers

# Plots
scatter_p <- mr_scatter_plot(mr_res, mr_dat) 
scatter_out_p <- scatter_p[[1]] + theme_bw() + 
  guides(color=guide_legend(nrow =1)) + 
  labs(title = "w/ outliers") +
  theme(
    text = element_text(size = 8), 
    legend.position = 'bottom'
  )

scatter_outlier_p <- mr_scatter_plot(mr_res_outlier, mr_dat_outlier)

scatter_outlier_out_p <- scatter_outlier_p[[1]] + theme_bw() + 
  guides(color=guide_legend(nrow =1)) + 
  labs(title = "w/o outliers") +
  theme(
    text = element_text(size = 8), 
    legend.position = 'bottom'
  )

radial_outlier_p <- plot_radial(radial_res_outlier, radial_scale = F, show_outliers = F)

cowplot::plot_grid(
  scatter_out_p + theme(legend.position = 'none'), 
  scatter_outlier_out_p + theme(legend.position = 'none'), 
  ivw_radial_p + coord_fixed(ratio=0.25) + theme(legend.position = 'none') + labs(title = "w/ outliers"),
  radial_outlier_p + coord_fixed(ratio=0.2) + theme(legend.position = 'none') + labs(title = "w/o outliers") 
)