
library(tidyverse)
library(survival)
library(survminer)
library(finalfit)

multiomics_clusters <- readRDS(file = "others/umap_optimized_clinical.Rds") %>% 
  dplyr::rename(id = sample)

## Clinical data
clinical <- readxl::read_xlsx("/Users/pol/Dropbox/compare_clusters/data/target/TARGET_AML_ClinicalData_Discovery_20221108.xlsx") %>% 
  mutate(id = gsub("TARGET-20-", "", `TARGET USI`)) %>% 
  right_join(multiomics_clusters, by = "id") %>% 
  filter(`First Event` %in% c("Censored", "Relapse")) %>% 
  dplyr::mutate(status_relapse = ifelse(`First Event` == "Censored", 0, 1),
                status_overall = ifelse(`Vital Status` == "Alive", 0, 1),
                time_relapse = `Event Free Survival Time in Days`/365,
                time_overall = `Overall Survival Time in Days`/365) %>% 
  dplyr::rename(Risk = risk_group,
                Age = `Age at Diagnosis in Days`,
                mutation = `Primary Cytogenetic Code`) %>% 
  # We carefully choose our reference groups for "Risk", "clust", and "mutation" (see below)
  dplyr::mutate(Risk = factor(Risk, levels = c("LR1", "LR2", "HR")),
                # clust = factor(clust, levels = c("others", "low_cluster", "high_cluster")),
                mutation = factor(mutation, levels = c("Normal", "inv(16)", "MLL", "t(8;21)", "Other", "Unknown"))) %>% 
  dplyr::select(id, clust = cluster, status_relapse, status_overall, time_relapse, time_overall, Risk, Age, mutation)

table(clinical$status_relapse, clinical$clust)

### Reference group for "Risk"
# table(clinical$Risk)
# 
# clinical %>%
#   dplyr::group_by(Risk) %>%
#   dplyr::summarise(median_relapse = median(time_relapse, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::arrange(median_relapse)

#### "Standard Risk" group is the normative category and contains a reasonable number of samples (46)

### Reference group for "clust"
table(clinical$clust)

clinical %>%
  dplyr::group_by(clust) %>%
  dplyr::summarise(median_relapse = median(time_relapse, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(median_relapse)

#### "Cluster_9" is the largest group (29 samples) and is also the group with a median relapse time 
#### closer to the dataset median relapse time 

### Reference group for "mutation"
# table(clinical$mutation)
# 
# clinical %>%
#   dplyr::group_by(mutation) %>%
#   dplyr::summarise(median_relapse = median(time_relapse, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::arrange(median_relapse)

#### "Normal" group is the normative category and contains a reasonable number of samples (21)

# Cox Proportional Hazard models 
### Univariate
fit_relapse_clu <- coxph(Surv(time_relapse, status_relapse) ~ clust, data = clinical) %>%
  finalfit::fit2df() %>% 
  mutate(group = "RF")

fit_relapse_ris <- coxph(Surv(time_relapse, status_relapse) ~ Risk, data = clinical) %>%
  finalfit::fit2df() %>% 
  mutate(group = "RF")

fit_overall_clu <- coxph(Surv(time_overall, status_overall) ~ clust, data = clinical) %>%
  finalfit::fit2df() %>% 
  mutate(group = "OS")

fit_overall_ris <- coxph(Surv(time_overall, status_overall) ~ Risk, data = clinical) %>%
  finalfit::fit2df() %>% 
  mutate(group = "OS")

fit_univ_results <- bind_rows(fit_relapse_clu, fit_relapse_ris, fit_overall_clu, fit_overall_ris) %>% 
  dplyr::mutate(explanatory = gsub("clust", "", explanatory),
                explanatory = gsub("Risk", "", explanatory),
                explanatory = gsub(" ", " Risk", explanatory),
                explanatory = gsub("Unknown", "Unknown Risk", explanatory),
                lwr = sub(".*\\(", "", HR),
                lwr = as.numeric(sub("\\-.*", "", lwr)),
                upr = sub(".*\\-", "", HR),
                upr = as.numeric(sub("\\,.*", "", upr)),
                pval = sub(".*, ", "", HR),
                pval = gsub("\\)", "", pval),
                pval = gsub("p=", "p = ", pval),
                pval = gsub("p<", "p < ", pval),
                HR = as.numeric(sub(" .*", "", HR)),
                log_hr = log(HR),
                log_lwr = log(lwr),
                log_upr = log(upr))

fit_univ_results_plot <- fit_univ_results #%>% 
# filter(explanatory %in% c("High Risk", "Low Risk", "Cluster_8", "Cluster_1"))

ggplot(fit_univ_results_plot, aes(log_hr, explanatory, color = group, label = pval)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") + # log(1) 
  geom_point(size = 2, position = position_dodge(width = 1)) +
  geom_errorbarh(aes(xmin = log_lwr, xmax = log_upr), height = 0.2, position = position_dodge(width = 1)) +
  geom_text(aes(log_upr - 0.5, group = group), vjust = -1, size = 5, color = "black", position = position_dodge(width = 1)) +
  labs(x = "log(Hazard Ratio)",
       y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = c(0.15, 0.08),
        axis.text = element_text(color = "black", size = 15),
        legend.title = element_blank()) +
  scale_color_manual(values = ggsci::pal_npg(alpha = 1)(2))

# ggsave(filename = "clone/results_target/survival_hr.png", width = 5)

### Multivariate
# We don't include mutations in the multivariate model because the "Risk Group" variable is based on them. 

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2214754/

# Risk stratification typically classifies as favorable those patients with Down syndrome or
# with AML characterized by t(8;21), t(15;17), or inv(16) cytogenetic abnormalities and rapid 
# early response to induction therapy. Unfavorable features include high white blood cell count, 
# −7/7q−, −5/5q−, or complex cytogenetics, and slow or no early response. The emerging consensus 
# is that patients with favorable AML do not benefit from MRD BMT in first remission.

fit_relapse <- coxph(Surv(time_relapse, status_relapse) ~ Age + Risk + clust, data = clinical) %>%
  finalfit::fit2df() %>% 
  mutate(group = "RF")

fit_overall <- coxph(Surv(time_overall, status_overall) ~ Age + Risk + clust, data = clinical) %>%
  finalfit::fit2df() %>% 
  mutate(group = "OS")

fit_multi_results <- bind_rows(fit_relapse, fit_overall) %>% 
  dplyr::mutate(explanatory = gsub("clust", "", explanatory),
                explanatory = gsub("Risk", "", explanatory),
                explanatory = gsub(" ", " Risk", explanatory),
                explanatory = gsub("GenderMale", "Gender: Male", explanatory),
                explanatory = gsub("mutation", "Mutation: ", explanatory),
                lwr = sub(".*\\(", "", HR),
                lwr = as.numeric(sub("\\-.*", "", lwr)),
                upr = sub(".*\\-", "", HR),
                upr = as.numeric(sub("\\,.*", "", upr)),
                pval = sub(".*, ", "", HR),
                pval = gsub("\\)", "", pval),
                pval = gsub("p=", "p = ", pval),
                pval = gsub("p<", "p < ", pval),
                HR = as.numeric(sub(" .*", "", HR)),
                log_hr = log(HR),
                log_lwr = log(lwr),
                log_upr = log(upr)) %>% 
  dplyr::filter(!(grepl("Unknown", explanatory) | grepl("Other", explanatory))) #%>% 
  # dplyr::filter(explanatory %in% c("low_cluster", "high_cluster", "High Risk", "Low Risk"))

ggplot(fit_multi_results, aes(log_hr, explanatory, label = pval)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") + # log(1)
  geom_point(size = 2, position = position_dodge(width = 1)) +
  geom_errorbarh(aes(xmin = log_lwr, xmax = log_upr), height = 0.2, position = position_dodge(width = 1)) +
  geom_text(aes(log_upr - 0.6), vjust = -1.5, color = "black",  size = 5, position = position_dodge(width = 1)) +
  labs(x = "log(Hazard Ratio)",
       y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.position = c(0.1, 0.08),
        axis.text = element_text(color = "black", size = 15),
        legend.title = element_blank()) +
  facet_wrap(~ group)

# ggsave(filename = "clone/results_target/survival_hr_multi.png", width = 8.5)

# Kaplan–Meier plots
clinical_sub <- clinical %>% 
  # dplyr::filter(clust %in% c("low_cluster", "high_cluster")) #%>%
dplyr::filter(clust %in% c("Cluster 1", "Cluster 2", "Cluster 8", "Cluster 5", "Cluster 3"))
# dplyr::filter(clust %in% c("Cluster_2", "Cluster_5"))

## Relapse
model_relapse <- survfit(Surv(time_relapse, status_relapse) ~ clust, data = clinical_sub)
names(model_relapse$strata) <- gsub("clust=", "", names(model_relapse$strata))

p <- ggsurvplot(
  model_relapse,
  # palette = ggsci::pal_npg(alpha = 1)(2),
  legend = "none",
  legend.title = element_blank(),
  surv.median.line = "none",
  data = clinical_sub, 
  pval = TRUE,
  pval.method = FALSE,
  risk.table = TRUE,
  tables.theme = theme_cleantable(),
  xlim = c(0, max(clinical_sub$time_relapse))
) +
  labs(x = "Time (years)",
       y = "Relapse Free Probability")

p$plot <- p$plot + theme(legend.title = element_text(size = 15, color = "black"),
                         legend.text = element_text(size = 15, color = "black"),
                         axis.text.x = element_text(size = 15, color = "black"),
                         axis.text.y = element_text(size = 15, color = "black"),
                         axis.title.x = element_text(size = 15, color = "black"),
                         axis.title.y = element_text(size = 15, color = "black"))

p$table$theme$text$size <- 15
p

# png("clone/results_target/survival_relapse_all.png", width = 6, height = 5, units = "in", res = 300)
# print(p, newpage = FALSE)
# dev.off()

## Overall survival
model_overall <- survfit(Surv(time_overall, status_overall) ~ clust, data = clinical_sub)
names(model_overall$strata) <- gsub("clust=", "", names(model_overall$strata))

p <- ggsurvplot(
  model_overall,
  # palette = ggsci::pal_npg(alpha = 1)(2),
  legend = "none",
  legend.title = element_blank(),
  surv.median.line = "none",
  data = clinical_sub, 
  pval = TRUE,
  pval.method = FALSE,
  risk.table = TRUE,
  tables.theme = theme_cleantable(),
  xlim = c(0, max(clinical_sub$time_overall))
) +
  labs(x = "Time (years)",
       y = "Overall Survival Probability")

p$plot <- p$plot + theme(legend.title = element_text(size = 15, color = "black"),
                         legend.text = element_text(size = 15, color = "black"),
                         axis.text.x = element_text(size = 15, color = "black"),
                         axis.text.y = element_text(size = 15, color = "black"),
                         axis.title.x = element_text(size = 15, color = "black"),
                         axis.title.y = element_text(size = 15, color = "black"))

p$table$theme$text$size <- 15
p

# png("clone/results_target/survival_overall_all.png", width = 6, height = 5, units = "in", res = 300)
# print(p, newpage = FALSE)
# dev.off()

# Compare Risk and Cluster
# clinical_sub <- clinical_sub %>% 
#   dplyr::mutate(Risk = factor(Risk, levels = c("Low Risk", "Standard Risk", "High Risk")))

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

ggplot(clinical_sub %>% drop_na(), aes(Risk, after_stat(count), fill = clust)) +
  geom_bar(position = ggplot2::position_dodge2()) +
  labs(x = NULL,
       y = "# Samples") +
  # scale_fill_manual(values = ggsci::pal_npg(alpha = 1)(2)) +
  scale_y_continuous(breaks = integer_breaks()) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.9),
        legend.title = element_blank(),
        text = element_text(size = 15),
        axis.text = element_text(size = 15, color = "black"),
        panel.grid = element_blank()) 

####
poma_obj <- POMA::PomaCreateObject(metadata = data.frame(id = multiomics_clusters$id, 
                                                         cluster = multiomics_clusters$clust), 
                                   features = do.call(cbind, omics))

limma_res <- POMA::PomaLimma(poma_obj, contrast = "high_cluster-low_cluster")

POMA::PomaBoxplots(poma_obj, x = "features", feature_name = limma_res$feature[1:10])

aaa <- ddh::make_cca_genes_table(input = list(content = c("MAK16", "NOP16", "KCNMB1")))

