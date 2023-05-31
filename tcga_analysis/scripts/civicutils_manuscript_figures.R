#!/usr/bin/env Rscript
################################################################################
## Figures for CIViCutils manuscript
## Author: Lourdes Rosano
## Date created: April 2022
## R Version: 4.2.2
################################################################################

library(optparse)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(stringr)


##################
## FUNCTIONS
##################

# Check that a given column name can be found in the provided header
check_column_name = function(column_name, header){
  if (sum(header %in% column_name) != 1){
    stop(paste0("Error! Column '", column_name, "' could not be found in the header of the provided input files!"))
    quit(save="no", status=1)
  }
}


# Check that a list of column names are all found in the provided header
check_column_names = function(column_name_list, header){
  for (column_name in column_name_list){
    check_column_name(column_name, header)
  }
}


# Check that two tables have identical headers
check_identical_headers = function(dataset_snv, dataset_cnv){
  if(!identical(colnames(dataset_snv), colnames(dataset_cnv))){
    stop("Error! Column names of '--infile_snv' and '--infile_cnv' are not identical!")
    quit(save="no", status=1)
  }  
}


# General helper function used for retrieving a subset of the input column names
prepare_data_for_plotting = function(dataset, column_names, header, sort_variable){
  check_column_names(column_names, header)
  column_names = c("sample_name", column_names)
  subdata = subset(dataset, select=column_names)
  df = melt(subdata, id="sample_name")
  df_sorted = arrange(df, sample_name)
  df_sorted = ddply(df_sorted, "sample_name", transform, label_ypos = value)
  # Sort samples by provided variable (decreasing order)
  sample_order = df_sorted[df_sorted$variable == sort_variable,,][order(df_sorted$value[df_sorted$variable == sort_variable], decreasing=T), "sample_name"]
  df_sorted$sample_order = factor(df_sorted$sample_name, levels=sample_order)
  return(df_sorted)
}


# General helper function for combining SNV and CNV datasets with identical internal structures
combine_snv_and_cnv_data_for_plotting = function(dataset_snv, dataset_cnv, column_names){
  check_column_names(column_names, colnames(dataset_snv))
  check_column_names(column_names, colnames(dataset_cnv))
  check_identical_headers(dataset_snv, dataset_cnv)
  
  # Keep track of which dataframe each entry originated from
  dataset_snv$dataset = "SNVs"
  dataset_cnv$dataset = "CNVs"
  
  # Re-sort samples in the CNV dataset by the order of the SNV dataset
  snv_sample_order = levels(dataset_snv$sample_order)
  cnv_sample_order = levels(dataset_cnv$sample_order)
  dataset_snv$sample_order = as.character(dataset_snv$sample_order)
  dataset_cnv$sample_order = as.character(dataset_cnv$sample_order)
  dataset_snv$sample_order = factor(dataset_snv$sample_order, levels=snv_sample_order)
  dataset_cnv$sample_order = factor(dataset_cnv$sample_order, levels=snv_sample_order)
  
  # Combine both dataframes
  df_combined = rbind(dataset_snv, dataset_cnv)
  df_combined$dataset = factor(df_combined$dataset, levels=c("SNVs", "CNVs"))
  return(df_combined)  
}


# General helper function to manually calculate outlier points per group in box plots
get_outliers_per_group = function(df, categories){
  df$index = 1:nrow(df)
  df$is_outlier = FALSE
  all_outliers = sapply(categories, function(category){
    subdf = subset(df, variable == category)
    lower_condition = subdf$value < quantile(subdf$value, 0.25) - 1.5 * IQR(subdf$value)
    upper_condition = subdf$value > quantile(subdf$value, 0.75) + 1.5 * IQR(subdf$value)
    is_out = lower_condition | upper_condition
    subdf$index[is_out]
  })
  all_indxs = unlist(all_outliers)
  df$is_outlier[all_indxs] = TRUE
  df = subset(df, select = -index)
  return(df)
}


# General helper function to manually calculate outlier points per group and dataset in box plots
get_outliers_per_group_and_dataset = function(df, categories){
  df_snv = subset(df, dataset == "SNVs")
  df_cnv = subset(df, dataset == "CNVs")
  df_snv_outliers = get_outliers_per_group(df_snv, categories)
  df_cnv_outliers = get_outliers_per_group(df_cnv, categories)
  df = rbind(df_snv_outliers, df_cnv_outliers)
  return(df)
}


# Helper function for generating Figure 1A
# Computes fraction of variants with CIViC information available per sample, in preparation for plotting
prepare_civic_fraction_data_for_plotting = function(dataset, header){
  n_civic_columns = c("sample_name", "all_variants", "all_civic_variants")
  check_column_names(n_civic_columns, header)
  subdata = subset(dataset, select=n_civic_columns)
  # Avoid having NaN in the results due to patients having all_variants=0
  skip_rows = subdata$all_variants!=0
  subdata$fraction = 0
  subdata[skip_rows,]$fraction = (subdata[skip_rows,]$all_civic_variants / subdata[skip_rows,]$all_variants)*100
  df = subset(subdata, select=c("sample_name", "fraction"))
  # Sort samples by computed fraction (decreasing order)
  sample_order = df[order(df$fraction, decreasing=T), "sample_name"]
  df$sample_order = factor(df$sample_name, levels=sample_order)
  return(df)
}


# Helper function for generating Figure 1B
# Aggregates numbers of CIViC variants across all samples in the cohort
aggregate_civic_data_across_patients = function(df_civic){
  # Sum all variants and CIViC variants across patients
  all_variants = sum(subset(df_civic, subset= variable == "all_variants", select = value))
  all_civic_variants = sum(subset(df_civic, subset= variable == "all_civic_variants", select = value))
  # Compute total fraction of variants with and without CIViC info across entire patient cohort
  fraction_civic_variants = (all_civic_variants / all_variants) * 100
  return(fraction_civic_variants)
}


# Helper function for generating Figure 2A
# Computes fractions of CIViC variants per tier class for each sample, in preparation for plotting
prepare_tier_fraction_data_for_plotting = function(dataset, header){
  n_tier_columns = c("sample_name", "all_civic_variants", "n_tier_1", "n_tier_1b", "n_tier_2", "n_tier_3")
  check_column_names(n_tier_columns, header)
  subdata = subset(dataset, select=n_tier_columns)
  # Avoid having NaN in the results due to patients having all_civic_variants=0
  skip_rows = subdata$all_civic_variants!=0
  subdata$fraction_tier_1 = 0
  subdata$fraction_tier_1b = 0
  subdata$fraction_tier_2 = 0
  subdata$fraction_tier_3 = 0
  subdata[skip_rows,]$fraction_tier_1 = (subdata[skip_rows,]$n_tier_1 / subdata[skip_rows,]$all_civic_variants)*100
  subdata[skip_rows,]$fraction_tier_1b = (subdata[skip_rows,]$n_tier_1b / subdata[skip_rows,]$all_civic_variants)*100
  subdata[skip_rows,]$fraction_tier_2 = (subdata[skip_rows,]$n_tier_2 / subdata[skip_rows,]$all_civic_variants)*100
  subdata[skip_rows,]$fraction_tier_3 = (subdata[skip_rows,]$n_tier_3 / subdata[skip_rows,]$all_civic_variants)*100
  new_subdata = subset(subdata, select=c("sample_name", "fraction_tier_1", "fraction_tier_1b", "fraction_tier_2", "fraction_tier_3"))
  new_subdata_sorted = new_subdata[with(new_subdata, order(fraction_tier_1, fraction_tier_1b, fraction_tier_2, fraction_tier_3)),]
  sample_order = new_subdata_sorted$sample_name
  
  # Reformat dataset in preparation for plotting
  df = melt(new_subdata, id="sample_name")
  df_sorted = arrange(df, sample_name) #, increasing = TRUE)
  df_sorted = ddply(df_sorted, "sample_name", transform, label_ypos = value)
  
  # Sort samples by provided variable  (decreasing order)
  df_sorted$sample_order = factor(df_sorted$sample_name, levels=sample_order)
  # Sort tiers by increasing order of priority
  tier_order = c("fraction_tier_3", "fraction_tier_2", "fraction_tier_1b", "fraction_tier_1")
  df_sorted$tier_order = factor(df_sorted$variable, levels=tier_order)
  return(df_sorted)
}


# Helper function for generating Figure 2B
# Aggregates absolute numbers of CIViC variants per tier class across all samples in the cohort
aggregate_tier_data_across_patients = function(df_tier){
  # Per tier, sum all variants across patients
  all_tier1_variants = sum(subset(df_tier, subset= variable == "n_tier_1", select = value))
  all_tier1b_variants = sum(subset(df_tier, subset= variable == "n_tier_1b", select = value))
  all_tier2_variants = sum(subset(df_tier, subset= variable == "n_tier_2", select = value))
  all_tier3_variants = sum(subset(df_tier, subset= variable == "n_tier_3", select = value))
  all_civic_variants = all_tier1_variants + all_tier1b_variants + all_tier2_variants + all_tier3_variants
  # Per tier, compute total fraction of variants across entire patient cohort
  fraction_tier1_variants = (all_tier1_variants / all_civic_variants) * 100
  fraction_tier1b_variants = (all_tier1b_variants / all_civic_variants) * 100
  fraction_tier2_variants = (all_tier2_variants / all_civic_variants) * 100
  fraction_tier3_variants = (all_tier3_variants / all_civic_variants) * 100
  return(list(tier1=fraction_tier1_variants, tier1b=fraction_tier1b_variants, tier2=fraction_tier2_variants, tier3=fraction_tier3_variants))
}


# Helper function for generating Figure 2B
# Aggregates fractions of CIViC variants per tier class across all samples in the cohort
aggregate_tier_fraction_data_across_patients = function(df_tier){
  # Per tier, sum all variants across patients
  tier1_fraction_sum = sum(subset(df_tier, subset= variable == "fraction_tier_1", select = value))
  tier1b_fraction_sum = sum(subset(df_tier, subset= variable == "fraction_tier_1b", select = value))
  tier2_fraction_sum = sum(subset(df_tier, subset= variable == "fraction_tier_2", select = value))
  tier3_fraction_sum = sum(subset(df_tier, subset= variable == "fraction_tier_3", select = value))
  n_samples = nrow(subset(df_tier, subset= variable == "fraction_tier_1", select = value))
  # Per tier, compute total fraction of variants across entire patient cohort
  mean_fraction_tier1 = (tier1_fraction_sum / n_samples)
  mean_fraction_tier1b = (tier1b_fraction_sum / n_samples)
  mean_fraction_tier2 = (tier2_fraction_sum / n_samples)
  mean_fraction_tier3 = (tier3_fraction_sum / n_samples)
  return(list(tier1=mean_fraction_tier1, tier1b=mean_fraction_tier1b, tier2=mean_fraction_tier2, tier3=mean_fraction_tier3))
}


# Helper function for generating Figure 4A
# Computes fractions of available disease names per ct class for each sample, in preparation for plotting
# Two versions of stats are provided: including vs. excluding tier3 variants
prepare_disease_fraction_data_for_plotting = function(dataset, header){
  n_disease_columns = c("sample_name", "n_diseases", "n_diseases_ct", "n_diseases_gt", "n_diseases_nct", 
                        "n_diseases_no_tier3", "n_diseases_ct_no_tier3", "n_diseases_gt_no_tier3", "n_diseases_nct_no_tier3")
  check_column_names(n_disease_columns, header)
  subdata = subset(dataset, select=n_disease_columns)
  
  # Avoid having NaN in the results due to patients having n_diseases=0
  skip_rows = subdata$n_diseases!=0
  subdata$fraction_ct = 0
  subdata$fraction_gt = 0
  subdata$fraction_nct = 0
  subdata[skip_rows,]$fraction_ct = (subdata[skip_rows,]$n_diseases_ct / subdata[skip_rows,]$n_diseases)*100
  subdata[skip_rows,]$fraction_gt = (subdata[skip_rows,]$n_diseases_gt / subdata[skip_rows,]$n_diseases)*100
  subdata[skip_rows,]$fraction_nct = (subdata[skip_rows,]$n_diseases_nct / subdata[skip_rows,]$n_diseases)*100
  
  # Also apply computation to version of stats excluding tier3 variants
  skip_rows_no_tier3 = subdata$n_diseases_no_tier3!=0
  subdata$fraction_ct_no_tier3 = 0
  subdata$fraction_gt_no_tier3 = 0
  subdata$fraction_nct_no_tier3 = 0
  subdata[skip_rows_no_tier3,]$fraction_ct_no_tier3 = (subdata[skip_rows_no_tier3,]$n_diseases_ct_no_tier3 / subdata[skip_rows_no_tier3,]$n_diseases_no_tier3)*100
  subdata[skip_rows_no_tier3,]$fraction_gt_no_tier3 = (subdata[skip_rows_no_tier3,]$n_diseases_gt_no_tier3 / subdata[skip_rows_no_tier3,]$n_diseases_no_tier3)*100
  subdata[skip_rows_no_tier3,]$fraction_nct_no_tier3 = (subdata[skip_rows_no_tier3,]$n_diseases_nct_no_tier3 / subdata[skip_rows_no_tier3,]$n_diseases_no_tier3)*100
  
  new_subdata = subset(subdata, select=c("sample_name", "fraction_ct", "fraction_gt", "fraction_nct", "fraction_ct_no_tier3", "fraction_gt_no_tier3", "fraction_nct_no_tier3"))
  new_subdata_sorted = new_subdata[with(new_subdata, order(fraction_ct, fraction_gt, fraction_nct, fraction_ct_no_tier3, fraction_gt_no_tier3, fraction_nct_no_tier3)),]
  sample_order = new_subdata_sorted$sample_name
  
  # Reformat dataset in preparation for plotting
  df = melt(new_subdata, id="sample_name")
  df_sorted = arrange(df, sample_name) #, increasing = TRUE)
  df_sorted = ddply(df_sorted, "sample_name", transform, label_ypos = value)
  
  # Sort samples by provided variable (decreasing order)
  df_sorted$sample_order = factor(df_sorted$sample_name, levels=sample_order)
  return(df_sorted)
}


# Helper function for generating Figure 4B
# Aggregates absolute numbers of disease names per ct class across all samples in the cohort
# Two versions of stats are provided: including vs. excluding tier3 variants
aggregate_disease_data_across_patients = function(df_disease){
  # Per ct class, sum all disease names across patients
  all_ct_diseases = sum(subset(df_disease, subset= variable == "n_diseases_ct", select = value))
  all_gt_diseases = sum(subset(df_disease, subset= variable == "n_diseases_gt", select = value))
  all_nct_diseases = sum(subset(df_disease, subset= variable == "n_diseases_nct", select = value))
  all_diseases = sum(subset(df_disease, subset= variable == "n_diseases", select = value))
  
  # Version of ct stats excluding tier3 variant matches
  all_ct_diseases_no_tier3 = sum(subset(df_disease, subset= variable == "n_diseases_ct_no_tier3", select = value))
  all_gt_diseases_no_tier3 = sum(subset(df_disease, subset= variable == "n_diseases_gt_no_tier3", select = value))
  all_nct_diseases_no_tier3 = sum(subset(df_disease, subset= variable == "n_diseases_nct_no_tier3", select = value))
  all_diseases_no_tier3 = sum(subset(df_disease, subset= variable == "n_diseases_no_tier3", select = value))
  
  # Per ct class, compute total fractions of ct classes across entire patient cohort
  fraction_ct_diseases = (all_ct_diseases / all_diseases) * 100
  fraction_gt_diseases = (all_gt_diseases / all_diseases) * 100
  fraction_nct_diseases = (all_nct_diseases / all_diseases) * 100
  fraction_ct_diseases_no_tier3 = (all_ct_diseases_no_tier3 / all_diseases_no_tier3) * 100
  fraction_gt_diseases_no_tier3 = (all_gt_diseases_no_tier3 / all_diseases_no_tier3) * 100
  fraction_nct_diseases_no_tier3 = (all_nct_diseases_no_tier3 / all_diseases_no_tier3) * 100
  return(list(ct=fraction_ct_diseases, gt=fraction_gt_diseases, nct=fraction_nct_diseases,
              ct_no_tier3=fraction_ct_diseases_no_tier3, gt_no_tier3=fraction_gt_diseases_no_tier3, nct_no_tier3=fraction_nct_diseases_no_tier3))
}


# Helper function for generating Figure 4B
# Aggregates fractions of disease names per ct class across all samples in the cohort
# Two versions of stats are provided: including vs. excluding tier3 variants
aggregate_disease_fraction_data_across_patients = function(df_disease){
  # Sum all disease names per ct class across patients
  ct_fraction_sum = sum(subset(df_disease, subset= variable == "fraction_ct", select = value))
  gt_fraction_sum = sum(subset(df_disease, subset= variable == "fraction_gt", select = value))
  nct_fraction_sum = sum(subset(df_disease, subset= variable == "fraction_nct", select = value))
  n_samples = nrow(subset(df_disease, subset= variable == "fraction_ct", select = value))
  # Compute total fraction of disease names across entire patient cohort
  mean_fraction_ct = (ct_fraction_sum / n_samples)
  mean_fraction_gt = (gt_fraction_sum / n_samples)
  mean_fraction_nct = (nct_fraction_sum / n_samples)
  # Version of drug stats excluding tier3 variant matches
  ct_fraction_sum_no_tier3 = sum(subset(df_disease, subset= variable == "fraction_ct_no_tier3", select = value))
  gt_fraction_sum_no_tier3 = sum(subset(df_disease, subset= variable == "fraction_gt_no_tier3", select = value))
  nct_fraction_sum_no_tier3 = sum(subset(df_disease, subset= variable == "fraction_nct_no_tier3", select = value))
  mean_fraction_ct_no_tier3 = (ct_fraction_sum_no_tier3 / n_samples)
  mean_fraction_gt_no_tier3 = (gt_fraction_sum_no_tier3 / n_samples)
  mean_fraction_nct_no_tier3 = (nct_fraction_sum_no_tier3 / n_samples)
  return(list(ct=mean_fraction_ct, gt=mean_fraction_gt, nct=mean_fraction_nct,
              ct_no_tier3=mean_fraction_ct_no_tier3, gt_no_tier3=mean_fraction_gt_no_tier3, nct_no_tier3=mean_fraction_nct_no_tier3))
}


# Helper function for generating Figure 5A
# Computes fraction of CIViC variants with drug info available for each sample, in preparation for plotting
# Two versions of stats are provided: including vs. excluding tier3 variants
prepare_drug_fraction_data_for_plotting = function(dataset, header){
  n_drug_columns = c("sample_name", "n_vars_drug_avail", "all_civic_variants", 
                     "n_vars_drug_avail_no_tier3", "all_civic_variants_no_tier3")
  check_column_names(n_drug_columns, header)
  subdata = subset(dataset, select=n_drug_columns)
  
  # Avoid having NaN in the results due to patients having all_civic_variants=0
  skip_rows = subdata$all_civic_variants!=0
  subdata$fraction = 0
  subdata[skip_rows,]$fraction = (subdata[skip_rows,]$n_vars_drug_avail / subdata[skip_rows,]$all_civic_variants)*100
  
  # Also apply computation to version of stats excluding tier3 variants
  skip_rows_no_tier3 = subdata$all_civic_variants_no_tier3!=0
  subdata$fraction_no_tier3 = 0
  subdata[skip_rows_no_tier3,]$fraction_no_tier3 = (subdata[skip_rows_no_tier3,]$n_vars_drug_avail_no_tier3 / subdata[skip_rows_no_tier3,]$all_civic_variants_no_tier3)*100
  
  new_subdata = subset(subdata, select=c("sample_name", "fraction", "fraction_no_tier3"))
  new_subdata_sorted = new_subdata[with(new_subdata, order(fraction, fraction_no_tier3)),]
  sample_order = new_subdata_sorted$sample_name
  
  # Reformat dataset in preparation for plotting
  df = melt(new_subdata, id="sample_name")
  df_sorted = arrange(df, sample_name) #, increasing = TRUE)
  df_sorted = ddply(df_sorted, "sample_name", transform, label_ypos = value)
  
  # Sort samples by provided variable (decreasing order)
  df_sorted$sample_order = factor(df_sorted$sample_name, levels=sample_order)
  return(df_sorted)
}


# Helper function for generating Figure 5B
# Aggregates absolute numbers of CIViC variants with drug info across all samples in the cohort
# Two versions of stats are provided: including vs. excluding tier3 variants
aggregate_drug_data_across_patients = function(df_drug){
  # Sum all CIViC variants and variants with consensus drug support across patients
  all_civic_variants = sum(subset(df_drug, subset= variable == "all_civic_variants", select = value))
  all_drug_variants = sum(subset(df_drug, subset= variable == "n_vars_drug_avail", select = value))
  # Version of drug stats excluding tier3 variant matches
  all_civic_variants_no_tier3 = sum(subset(df_drug, subset= variable == "all_civic_variants_no_tier3", select = value))
  all_drug_variants_no_tier3 = sum(subset(df_drug, subset= variable == "n_vars_drug_avail_no_tier3", select = value))
  # Compute total fraction of variants with and without CIViC info across entire patient cohort
  fraction_drug_variants = (all_drug_variants / all_civic_variants) * 100
  fraction_drug_variants_no_tier3 = (all_drug_variants_no_tier3 / all_civic_variants_no_tier3) * 100
  return(list(drug=fraction_drug_variants, drug_no_tier3=fraction_drug_variants_no_tier3))
}


# Helper function for generating Figure 5B
# Aggregates fractions of CIViC variants with drug info across all samples in the cohort
# Two versions of stats are provided: including vs. excluding tier3 variants
aggregate_drug_fraction_data_across_patients = function(df_drug){
  # Sum all variants across patients
  fraction_sum = sum(subset(df_drug, subset= variable == "fraction", select = value))
  n_samples = nrow(subset(df_drug, subset= variable == "fraction", select = value))
  # Compute total fraction of variants across entire patient cohort
  mean_fraction = (fraction_sum / n_samples)
  # Version of drug stats excluding tier3 variant matches
  fraction_sum_no_tier3 = sum(subset(df_drug, subset= variable == "fraction_no_tier3", select = value))
  mean_fraction_no_tier3 = (fraction_sum_no_tier3 / n_samples)
  return(list(drug=mean_fraction, drug_no_tier3=mean_fraction_no_tier3))
}


# Helper function for generating Figure 6A
# Computes mean number of consensus drug predictions per variant for each sample, in preparation for plotting
# Two versions of stats are provided: including vs. excluding tier3 variants
prepare_mean_n_consensus_data_for_plotting = function(dataset, header){
  n_consensus_columns = c("sample_name", "n_total_consensus", "n_vars_drug_avail", 
                          "n_total_consensus_no_tier3", "n_vars_drug_avail_no_tier3")
  check_column_names(n_consensus_columns, header)
  subdata = subset(dataset, select=n_consensus_columns)
  
  # Avoid having NaN in the results due to patients having n_vars_drug_avail=0
  skip_rows = subdata$n_vars_drug_avail!=0
  subdata$mean = 0
  subdata[skip_rows,]$mean = (subdata[skip_rows,]$n_total_consensus / subdata[skip_rows,]$n_vars_drug_avail)
  
  # Also apply computation to version of stats excluding tier3 variants
  skip_rows_no_tier3 = subdata$n_vars_drug_avail_no_tier3!=0
  subdata$mean_no_tier3 = 0
  subdata[skip_rows_no_tier3,]$mean_no_tier3 = (subdata[skip_rows_no_tier3,]$n_total_consensus_no_tier3 / subdata[skip_rows_no_tier3,]$n_vars_drug_avail_no_tier3)
  
  new_subdata = subset(subdata, select=c("sample_name", "mean", "mean_no_tier3"))
  new_subdata_sorted = new_subdata[with(new_subdata, order(mean, mean_no_tier3)),]
  sample_order = new_subdata_sorted$sample_name
  
  # Reformat dataset in preparation for plotting
  df = melt(new_subdata, id="sample_name")
  df_sorted = arrange(df, sample_name) #, increasing = TRUE)
  df_sorted = ddply(df_sorted, "sample_name", transform, label_ypos = value)
  
  # Sort samples by provided variable (decreasing order)
  df_sorted$sample_order = factor(df_sorted$sample_name, levels=sample_order)
  return(df_sorted)
}


# Helper function for generating Figure 6B
# Computes fractions of available consensus drug predictions per variant for each sample, in preparation for plotting
# Two versions of stats are provided: including vs. excluding tier3 variants
prepare_consensus_fraction_data_for_plotting = function(dataset, header){
  n_consensus_columns = c("sample_name", "n_total_consensus", "n_total_support", "n_total_resistance", "n_total_conflict", "n_total_unknown",
                          "n_total_consensus_no_tier3", "n_total_support_no_tier3", "n_total_resistance_no_tier3", "n_total_conflict_no_tier3", "n_total_unknown_no_tier3")
  check_column_names(n_consensus_columns, header)
  subdata = subset(dataset, select=n_consensus_columns)
  # Avoid having NaN in the results due to patients having n_total_consensus=0
  skip_rows = subdata$n_total_consensus!=0
  subdata$fraction_support = 0
  subdata$fraction_resistance = 0
  subdata$fraction_conflict = 0
  subdata$fraction_unknown = 0
  subdata[skip_rows,]$fraction_support = (subdata[skip_rows,]$n_total_support / subdata[skip_rows,]$n_total_consensus)*100
  subdata[skip_rows,]$fraction_resistance = (subdata[skip_rows,]$n_total_resistance / subdata[skip_rows,]$n_total_consensus)*100
  subdata[skip_rows,]$fraction_conflict = (subdata[skip_rows,]$n_total_conflict / subdata[skip_rows,]$n_total_consensus)*100
  subdata[skip_rows,]$fraction_unknown = (subdata[skip_rows,]$n_total_unknown / subdata[skip_rows,]$n_total_consensus)*100
  # Version of stats excluding tier3 variants
  skip_rows_no_tier3 = subdata$n_total_consensus_no_tier3!=0
  subdata$fraction_support_no_tier3 = 0
  subdata$fraction_resistance_no_tier3 = 0
  subdata$fraction_conflict_no_tier3 = 0
  subdata$fraction_unknown_no_tier3 = 0
  subdata[skip_rows_no_tier3,]$fraction_support_no_tier3 = (subdata[skip_rows_no_tier3,]$n_total_support_no_tier3 / subdata[skip_rows_no_tier3,]$n_total_consensus_no_tier3)*100
  subdata[skip_rows_no_tier3,]$fraction_resistance_no_tier3 = (subdata[skip_rows_no_tier3,]$n_total_resistance_no_tier3 / subdata[skip_rows_no_tier3,]$n_total_consensus_no_tier3)*100
  subdata[skip_rows_no_tier3,]$fraction_conflict_no_tier3 = (subdata[skip_rows_no_tier3,]$n_total_conflict_no_tier3 / subdata[skip_rows_no_tier3,]$n_total_consensus_no_tier3)*100
  subdata[skip_rows_no_tier3,]$fraction_unknown_no_tier3 = (subdata[skip_rows_no_tier3,]$n_total_unknown_no_tier3 / subdata[skip_rows_no_tier3,]$n_total_consensus_no_tier3)*100
  
  new_subdata = subset(subdata, select=c("sample_name", "fraction_support", "fraction_resistance", "fraction_conflict", "fraction_unknown",
                                         "fraction_support_no_tier3", "fraction_resistance_no_tier3", "fraction_conflict_no_tier3", "fraction_unknown_no_tier3"))
  new_subdata_sorted = new_subdata[with(new_subdata, order(fraction_support, fraction_resistance, fraction_conflict, fraction_unknown)),]
  sample_order = new_subdata_sorted$sample_name
  
  # Reformat dataset in preparation for plotting
  df = melt(new_subdata, id="sample_name")
  df_sorted = arrange(df, sample_name) #, increasing = TRUE)
  df_sorted = ddply(df_sorted, "sample_name", transform, label_ypos = value)
  
  # Sort samples by provided variable (decreasing order)
  df_sorted$sample_order = factor(df_sorted$sample_name, levels=sample_order)
  return(df_sorted)
}



##################
## MAIN SCRIPT
##################


# Parse command line arguments
option_list = list(
  make_option("--infile_snv", type = "character", help = "Input table listing diverse statistics from processing CIViCutils SNV annotations reported for a set of samples."),
  make_option("--infile_cnv", type = "character", help = "Input table listing diverse statistics from processing CIViCutils CNV annotations reported for a set of samples."),
  make_option("--outfile_tag", type = "character", help = "Name prefix of the output figures.")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


data_snv = read.table(opt$infile_snv, header=TRUE, sep="\t", check.names = FALSE)
data_cnv = read.table(opt$infile_cnv, header=TRUE, sep="\t", check.names = FALSE)

# Sanity check that the two tables have identical headers
check_identical_headers(data_snv, data_cnv)

# Sanity check that required input column can be found
input_colnames = colnames(data_snv)
check_column_name("sample_name", input_colnames)

# Sanity check for identical order of rows (sample names) in both tables
samples_names_snv = data_snv$sample_name
samples_names_cnv = data_cnv$sample_name
if (!identical(samples_names_snv, samples_names_cnv)){
  stop("Error! Order of samples in column 'sample_name' of the provided input files are not identical!")
  quit(save="no", status=1)
}



## S1) BOX PLOT NUMBER OF VARIANTS MATCHED IN CIVIC

# Use fraction of total variants having CIViCutils results
# Sanity check that required input columns can be found
df_civic_snv_fraction = prepare_civic_fraction_data_for_plotting(data_snv, input_colnames)
df_civic_cnv_fraction = prepare_civic_fraction_data_for_plotting(data_cnv, input_colnames)

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
plotting_columns_ps1 = c("sample_name", "fraction", "sample_order")
df_civic_fraction =  combine_snv_and_cnv_data_for_plotting(df_civic_snv_fraction, df_civic_cnv_fraction, plotting_columns_ps1)


# Combined version of box plot
# Showing SNVs and CNVs as facets
outfile_ps1 = paste0(opt$outfile_tag, ".s1.boxplot_fraction_variants_in_civic.pdf")
# outfile_ps1 = paste0(opt$outfile_tag, ".s1.boxplot_fraction_variants_in_civic.png")
ps1 = ggplot(df_civic_fraction, aes(x=dataset, y=fraction, fill=dataset)) +
  geom_boxplot(outlier.shape=NA, outlier.size=0, notch=FALSE, width=0.4) +
  # geom_jitter(aes(color=dataset), size=1, alpha=0.9, position=position_jitter(seed=123)) +
  # position=position_jitterdodge(jitter.width=0.65, jitter.height=0.65, seed=123)
  ylab("Percent of variants matched in CIViC") + 
  xlab("") +
  theme(axis.text = element_text(size=22, color="black"),
        axis.title = element_text(size=24, color="black")) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("grey", "grey")) +
  scale_color_manual(values=c("black", "black"))
# plot(ps1)
ggsave(outfile_ps1, ps1, width=10, height=8, units="in", dpi=1200)



## 2A) PIE CHART TOTAL NUMBER OF VARIANTS MATCHED IN CIVIC ACROSS PATIENTS

# Use absolute variant numbers
# Sanity check that required input columns can be found
n_civic_columns = c("all_variants", "all_civic_variants")
df_civic_snv_abs = prepare_data_for_plotting(data_snv, n_civic_columns, input_colnames, "all_variants")
df_civic_cnv_abs = prepare_data_for_plotting(data_cnv, n_civic_columns, input_colnames, "all_variants")

# Aggregate absolute numbers across patients in the cohort
fraction_civic_snv = aggregate_civic_data_across_patients(df_civic_snv_abs)
# Compute mean fraction for SNVs
mean_fraction_civic_snv = (sum(df_civic_snv_fraction$fraction) / nrow(df_civic_snv_fraction))

# Aggregate absolute numbers across patients in the cohort
fraction_civic_cnv = aggregate_civic_data_across_patients(df_civic_cnv_abs)
# Compute mean fraction for CNVs
mean_fraction_civic_cnv = (sum(df_civic_cnv_fraction$fraction) / nrow(df_civic_cnv_fraction))

# Manually combine data for pie chart
df_civic_pie = data.frame(variant_type = rep(c("No match", "Matched in CIViC"), 2),
                          dataset = c(rep("SNVs", 2), rep("CNVs", 2)),
                          value = c((100 - fraction_civic_snv), fraction_civic_snv, (100 - fraction_civic_cnv), fraction_civic_cnv))

# Assign percentages as labels for the pie chart
df_civic_pie$variant_type = factor(df_civic_pie$variant_type, levels=c("No match", "Matched in CIViC"))
df_civic_pie$dataset = factor(df_civic_pie$dataset, levels=c("SNVs", "CNVs"))
df_civic_pie$label = percent(df_civic_pie$value/100, accuracy=0.1)


# Pie chart
# Text labels introduced post-generation
outfile_p2a = paste0(opt$outfile_tag, ".2a.piechart_variants_in_civic_across_patients.pdf")
# outfile_p2a = paste0(opt$outfile_tag, ".2a.piechart_variants_in_civic_across_patients.png")
p2a = ggplot(df_civic_pie, aes(x="", y=value, fill=variant_type)) +
  geom_bar(width=1, stat="identity", color="black") +
  facet_grid(. ~ dataset) +
  coord_polar("y", start=0) +
  ylab(paste0("Total patients = ", nrow(df_civic_snv_fraction))) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size=40, hjust=0.5),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=35),
        strip.text = element_text(size=40)) +
  scale_fill_manual(values=c("grey", "#0072B2"), name="")
# plot(p2a)
ggsave(outfile_p2a, p2a, width=10, height=8, units="in", dpi=1200)



## S2) BOX PLOT DISTRIBUTION OF TIER CLASSIFICATION (only for variants matched in CIVIC)

# Use fraction of total variants having CIViCutils results
df_tier_snv_fraction = prepare_tier_fraction_data_for_plotting(data_snv, input_colnames)
df_tier_cnv_fraction = prepare_tier_fraction_data_for_plotting(data_cnv, input_colnames)

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
plotting_columns_ps2 = c("sample_name", "variable", "value", "label_ypos", "sample_order", "tier_order")
df_tier_fraction =  combine_snv_and_cnv_data_for_plotting(df_tier_snv_fraction, df_tier_cnv_fraction, plotting_columns_ps2)

# Choose color scheme for plot
tier_cols = c("fraction_tier_1" = "black", "fraction_tier_1b" = "black", "fraction_tier_2" = "black", "fraction_tier_3" = "black")
# tier_cols = c("fraction_tier_1" = "red", "fraction_tier_1b" = "#0072B2", "fraction_tier_2" = "gold", "fraction_tier_3" = "skyblue")

# Combined version of box plot
# Showing SNVs and CNVs as facets
outfile_ps2 = paste0(opt$outfile_tag, ".s2.boxplot_fraction_variants_per_tier.pdf")
# outfile_ps2 = paste0(opt$outfile_tag, ".s2.boxplot_fraction_variants_per_tier.png")
ps2 = ggplot(df_tier_fraction, aes(x=variable, y=value)) +
  geom_boxplot(fill="grey", outlier.shape=NA, outlier.size=0, notch=FALSE, width=0.3) +
  geom_jitter(aes(color=variable), size=1, alpha=0.9, position=position_jitter(seed=123)) +
  facet_grid(dataset ~ .) +
  ylab("Percent of variants matched in CIViC") + 
  xlab("") +
  theme(axis.text = element_text(size=22, color="black"),
        axis.title = element_text(size=24, color="black")) +
  theme(legend.justification = c(1, 1),
        legend.position = "none",
        strip.text = element_text(size=20)) +
  scale_x_discrete(breaks=c("fraction_tier_1", "fraction_tier_1b", "fraction_tier_2", "fraction_tier_3"), 
                   labels=c("Tier 1", "Tier 1b", "Tier 2", "Tier 3")) +
  scale_color_manual(name="Tier class", breaks=c("fraction_tier_1", "fraction_tier_1b", "fraction_tier_2", "fraction_tier_3"), 
                     labels=c("Tier 1", "Tier 1b", "Tier 2", "Tier 3"), values=tier_cols)
# plot(ps2)
ggsave(outfile_ps2, ps2, width=10, height=8, units="in", dpi=1200)



## 2B) PIE CHART TOTAL DISTRIBUTION OF TIER CLASSIFICATION ACROSS PATIENTS

# Use absolute variant numbers
# Sanity check that required input columns can be found
n_tier_columns = c("n_tier_1", "n_tier_1b", "n_tier_2", "n_tier_3", "all_civic_variants")
df_tier_snv_abs = prepare_data_for_plotting(data_snv, n_tier_columns, input_colnames, "all_civic_variants")
df_tier_cnv_abs = prepare_data_for_plotting(data_cnv, n_tier_columns, input_colnames, "all_civic_variants")

# Sort tiers by increasing order of priority
tier_order = c("n_tier_3", "n_tier_2", "n_tier_1b", "n_tier_1")
df_tier_snv_abs$tier_order = factor(df_tier_snv_abs$variable, levels=tier_order)
df_tier_cnv_abs$tier_order = factor(df_tier_cnv_abs$variable, levels=tier_order)

# Remove number of variants matched in CIViC (variable used only for sorting samples)
df_tier_snv_abs = df_tier_snv_abs[df_tier_snv_abs$variable != "all_civic_variants",,drop=F]
df_tier_cnv_abs = df_tier_cnv_abs[df_tier_cnv_abs$variable != "all_civic_variants",,drop=F]

# Aggregate absolute numbers across patients in the cohort
fractions_per_tier_snv = aggregate_tier_data_across_patients(df_tier_snv_abs)
# Compute mean fraction for SNVs
mean_fractions_per_tier_snv = aggregate_tier_fraction_data_across_patients(df_tier_snv_fraction)

# Aggregate absolute numbers across patients in the cohort
fractions_per_tier_cnv = aggregate_tier_data_across_patients(df_tier_cnv_abs)
# Compute mean fraction for CNVs
mean_fractions_per_tier_cnv = aggregate_tier_fraction_data_across_patients(df_tier_cnv_fraction)

# Manually combine data for pie chart
df_tier_pie = data.frame(tier = rep(c("Tier 1", "Tier 1b", "Tier 2", "Tier 3"), 2),
                         dataset = c(rep("SNVs", 4), rep("CNVs", 4)),
                         value = c(unlist(fractions_per_tier_snv), unlist(fractions_per_tier_cnv)))

# Assign percentages as labels for the pie chart
df_tier_pie$dataset = factor(df_tier_pie$dataset, levels=c("SNVs", "CNVs"))
df_tier_pie$label = percent(df_tier_pie$value/100, accuracy=0.1)
# Show all available classes (even if 0% frequency)
df_tier_pie$tier = factor(df_tier_pie$tier, levels=c("Tier 1", "Tier 1b", "Tier 2", "Tier 3"))

# Choose color scheme for plot
tier_cols = c("Tier 1"="red", "Tier 1b"="#0072B2", "Tier 2"="gold", "Tier 3"="skyblue")


# Pie chart
# Text labels introduced post-generation
outfile_p2b = paste0(opt$outfile_tag, ".2b.piechart_variants_per_tier_across_patients.pdf")
# outfile_p2b = paste0(opt$outfile_tag, ".2b.piechart_variants_per_tier_across_patients.png")
p2b = ggplot(df_tier_pie, aes(x="", y=value, fill=tier)) +
  geom_bar(width=1, stat="identity", color="black") +
  facet_grid(. ~ dataset) +
  coord_polar("y", start=0) +
  ylab(paste0("Total patients = ", nrow(df_civic_snv_fraction))) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size=40, hjust=0.5),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=35),
        strip.text = element_text(size=40)) +
  scale_fill_manual(values=tier_cols, name="")
# plot(p2b)
ggsave(outfile_p2b, p2b, width=10, height=8, units="in", dpi=1200)



## S3) BOX PLOT DISTRIBUTION OF EVIDENCE TYPES (only for variants matched in CIVIC)

# Use absolute variant numbers
# Note: critical in this case to avoid interpretation bias (as the same variant can be accounted for across several evidence types)

# Sanity check that required input columns can be found
n_evidence_columns = c("n_predictive_vars", "n_diagnostic_vars", "n_prognostic_vars", "n_predisposing_vars",
                       "n_predictive_vars_no_tier3", "n_diagnostic_vars_no_tier3", "n_prognostic_vars_no_tier3", "n_predisposing_vars_no_tier3")
df_evidence_snv_abs = prepare_data_for_plotting(data_snv, n_evidence_columns, input_colnames, "n_predictive_vars_no_tier3")
df_evidence_cnv_abs = prepare_data_for_plotting(data_cnv, n_evidence_columns, input_colnames, "n_predictive_vars_no_tier3")

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
plotting_columns = c("sample_name", "variable", "value", "label_ypos", "sample_order")
df_evidence =  combine_snv_and_cnv_data_for_plotting(df_evidence_snv_abs, df_evidence_cnv_abs, plotting_columns)

# Extract relevant features from the variable names in preparation for plotting
df_evidence$tier3_variants = ifelse(grepl("_no_tier3", df_evidence$variable), "Excluding", "Including")
df_evidence$tier3_variants = factor(df_evidence$tier3_variants, levels=c("Including", "Excluding"))
df_evidence$evidence_type = gsub("_no_tier3", "", gsub("_vars", "", gsub("n_", "", df_evidence$variable)))
# Use camel case for displaying the Evidence Types in the facets of the plot
df_evidence$evidence_type = str_to_title(df_evidence$evidence_type)
df_evidence$evidence_type = as.character(df_evidence$evidence_type)
df_evidence$evidence_type = factor(df_evidence$evidence_type, levels=c("Predictive", "Prognostic", "Diagnostic", "Predisposing"))

# Choose color scheme for plot
tier3_colors = c("Including" = "grey", "Excluding" = "#0072B2")


# Combined version of box plot
# Showing Evidence types in X-axis (variant numbers grouped by tier3 exclusion/inclusion)
# Showing SNVs and CNVs as facets in Y-axis
outfile_ps3 = paste0(opt$outfile_tag, ".s3.boxplot_n_variants_per_evidence_type.pdf")
# outfile_ps3 = paste0(opt$outfile_tag, ".s3.boxplot_n_variants_per_evidence_type.png")
ps3 = ggplot(df_evidence, aes(x=evidence_type, y=value, fill=tier3_variants)) +
  geom_boxplot(outlier.shape=NA, outlier.size=0, notch=FALSE) +
  facet_grid(dataset ~ .) +
  # coord_cartesian(ylim = c(0, 40)) +    # version with correct scale for CNVs
  # coord_cartesian(ylim = c(0, 15)) +    # version with correct scale for SNVs
  ylab("Number of variants") +
  xlab("") +
  theme(axis.text = element_text(size=18, color="black"),
        axis.title = element_text(size=24, color="black")) +
  theme(legend.justification = c(1, 1),
        legend.text = element_text(size=20, color="black"),
        legend.title = element_text(size=22, color="black"),
        strip.text = element_text(size=20)) +
  scale_fill_manual(name="Tier3 variants", values=tier3_colors) +
  scale_color_manual(name="Tier3 variants", values=tier3_colors)
# plot(ps3)
ggsave(outfile_ps3, ps3, width=10, height=8, units="in", dpi=1200)



## 3A) BOX PLOT DISEASE DISTRIBUTION OF CT CLASSIFICATION

# Use fraction of total diseases across CIViCutils results
df_disease_snv_fraction = prepare_disease_fraction_data_for_plotting(data_snv, input_colnames)
df_disease_cnv_fraction = prepare_disease_fraction_data_for_plotting(data_cnv, input_colnames)

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
df_disease_fraction =  combine_snv_and_cnv_data_for_plotting(df_disease_snv_fraction, df_disease_cnv_fraction, plotting_columns)

# Extract relevant features from the variable names in preparation for plotting
df_disease_fraction$ct_class = gsub("_no_tier3", "", gsub("fraction_", "", df_disease_fraction$variable))
df_disease_fraction$ct_class = factor(df_disease_fraction$ct_class, levels=c("ct", "gt", "nct"))
df_disease_fraction$tier3_variants = ifelse(grepl("_no_tier3", df_disease_fraction$variable), "Excluding", "Including")
df_disease_fraction$tier3_variants = factor(df_disease_fraction$tier3_variants, levels=c("Including", "Excluding"))

# Retrieve subset of outlier points for each group being plotted
disease_fraction_categories = names(table(df_disease_fraction$variable))
df_disease_fraction_outliers = get_outliers_per_group_and_dataset(df_disease_fraction, disease_fraction_categories)

# Include dummy outlier point from missing category to ensure the position and jittering of points work correctly
df_disease_fraction_outliers$is_outlier[5] = TRUE


# Combined version of box plot
# Showing ct classes in X-axis (disease numbers grouped by tier3 exclusion/inclusion)
# Showing SNVs and CNVs as facets in Y-axis
outfile_p3a = paste0(opt$outfile_tag, ".3a.boxplot_fraction_diseases_per_ct.pdf")
# outfile_p3a = paste0(opt$outfile_tag, ".3a.boxplot_fraction_diseases_per_ct.png")
p3a = ggplot(df_disease_fraction_outliers, aes(x=ct_class, y=value, fill=tier3_variants)) +
  geom_boxplot(outlier.shape=NA, outlier.size=0, notch=FALSE) +
  # geom_point(position=position_jitterdodge(jitter.width=0.65, jitter.height=0.65, seed=123),
  #            aes(color=tier3_variants), size=1, alpha=0.9) +
  geom_point(data=subset(df_disease_fraction_outliers, is_outlier==TRUE), 
             position=position_jitterdodge(jitter.width=0.50, jitter.height=0.50, seed=123),
             aes(color=tier3_variants), size=1, alpha=0.9) +
  facet_grid(dataset ~ .) +
  ylab("Percent of disease names") + 
  xlab("") +
  theme(axis.text = element_text(size=28, color="black"),
        axis.title = element_text(size=32, color="black"),
        strip.text = element_text(size=24)) +
  theme(legend.justification = c(1, 1),
        legend.text = element_text(size=26, color="black"),
        legend.title = element_text(size=30, color="black")) +
  scale_fill_manual(name="Tier3 variants", values=tier3_colors) +
  scale_color_manual(name="Tier3 variants", values=tier3_colors)
# plot(p3a)
ggsave(outfile_p3a, p3a, width=10, height=8, units="in", dpi=1200)



## 3B) PIE CHART TOTAL DISTRIBUTION OF CT CLASSES ACROSS PATIENTS

# Use absolute disease numbers
# Sanity check that required input columns can be found
n_disease_columns = c("n_diseases", "n_diseases_ct", "n_diseases_gt", "n_diseases_nct", 
                      "n_diseases_no_tier3", "n_diseases_ct_no_tier3", "n_diseases_gt_no_tier3", "n_diseases_nct_no_tier3")
df_disease_snv_abs = prepare_data_for_plotting(data_snv, n_disease_columns, input_colnames, "n_diseases")
df_disease_cnv_abs = prepare_data_for_plotting(data_cnv, n_disease_columns, input_colnames, "n_diseases")

# Aggregate absolute numbers across patients in the cohort
fraction_disease_snv = aggregate_disease_data_across_patients(df_disease_snv_abs)
# Compute mean fraction for SNVs
mean_fraction_disease_snv = aggregate_disease_fraction_data_across_patients(df_disease_snv_fraction)

# Aggregate absolute numbers across patients in the cohort
fraction_disease_cnv = aggregate_disease_data_across_patients(df_disease_cnv_abs)
# Compute mean fraction for CNVs
mean_fraction_disease_cnv = aggregate_disease_fraction_data_across_patients(df_disease_cnv_fraction)

# Manually combine data for pie chart
df_disease_pie = data.frame(ct_class = rep(c("ct", "gt", "nct"), 4),
                            dataset = c(rep("SNVs", 3), rep("CNVs", 3), rep("SNVs", 3), rep("CNVs", 3)),
                            tier3_variants = c(rep("Including\ntier3", 6), rep("Excluding\ntier3", 6)),
                            value = c(fraction_disease_snv$ct, fraction_disease_snv$gt, fraction_disease_snv$nct,
                                      fraction_disease_cnv$ct, fraction_disease_cnv$gt, fraction_disease_cnv$nct,
                                      fraction_disease_snv$ct_no_tier3, fraction_disease_snv$gt_no_tier3, fraction_disease_snv$nct_no_tier3,
                                      fraction_disease_cnv$ct_no_tier3, fraction_disease_cnv$gt_no_tier3, fraction_disease_cnv$nct_no_tier3))

# Assign percentages as labels for the pie chart
df_disease_pie$ct_class = factor(df_disease_pie$ct_class, levels=c("ct", "gt", "nct"))
df_disease_pie$dataset = factor(df_disease_pie$dataset, levels=c("SNVs", "CNVs"))
df_disease_pie$tier3_variants = factor(df_disease_pie$tier3_variants, levels=c("Including\ntier3", "Excluding\ntier3"))
df_disease_pie$label = percent(df_disease_pie$value/100, accuracy=0.1)

# Choose color scheme for plot
disease_cols = c("ct" = "red", "gt" = "gold", "nct" = "skyblue")


# Pie chart
# Text labels introduced post-generation
outfile_p3b = paste0(opt$outfile_tag, ".3b.piechart_diseases_per_ct_across_patients.pdf")
# outfile_p3b = paste0(opt$outfile_tag, ".3b.piechart_diseases_per_ct_across_patients.png")
p3b = ggplot(df_disease_pie, aes(x="", y=value, fill=ct_class)) +
  geom_bar(width=1, stat="identity", color="black") +
  facet_grid(tier3_variants ~ dataset, switch = "y") +
  coord_polar("y", start=0) +
  ylab(paste0("Total patients = ", nrow(df_civic_snv_fraction))) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size=44, hjust=0.5),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=40),
        legend.title = element_text(size=30),
        strip.text.x = element_text(size=40),
        strip.text.y.left = element_text(size=44, angle=0),
        panel.spacing = unit(1, "lines")) +
  scale_fill_manual(name="Cancer specificity", values=disease_cols)
# plot(p3b)
ggsave(outfile_p3b, p3b, width=10, height=8, units="in", dpi=1200)



## 4A) BOX PLOT DISTRIBUTION OF CONSENSUS DRUG SUPPORT

# Use fraction of CIViC variants having consensus drug support info available
df_drug_snv_fraction = prepare_drug_fraction_data_for_plotting(data_snv, input_colnames)
df_drug_cnv_fraction = prepare_drug_fraction_data_for_plotting(data_cnv, input_colnames)

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
df_drug_fraction =  combine_snv_and_cnv_data_for_plotting(df_drug_snv_fraction, df_drug_cnv_fraction, plotting_columns)

# Extract relevant features from the variable names in preparation for plotting
df_drug_fraction$tier3_variants = ifelse(grepl("_no_tier3", df_drug_fraction$variable), "Excluding", "Including")
df_drug_fraction$tier3_variants = factor(df_drug_fraction$tier3_variants, levels=c("Including", "Excluding"))

# Retrieve subset of outlier points for each group being plotted
drug_fraction_categories = names(table(df_drug_fraction$variable))
df_drug_fraction_outliers = get_outliers_per_group_and_dataset(df_drug_fraction, drug_fraction_categories)


# Combined version of box plot
# Showing SNVs and CNVs in X-axis (variant numbers grouped by tier3 exclusion/inclusion)
outfile_p4a = paste0(opt$outfile_tag, ".4a.boxplot_fraction_variants_with_drug.pdf")
# outfile_p4a = paste0(opt$outfile_tag, ".4a.boxplot_fraction_variants_with_drug.png")
p4a = ggplot(df_drug_fraction_outliers, aes(x=dataset, y=value, fill=tier3_variants)) +
  geom_boxplot(outlier.shape=NA, outlier.size=0, notch=FALSE) +
  # geom_point(position=position_jitterdodge(jitter.width=0.65, jitter.height=0.65, seed=123), 
  #            aes(color=tier3_variants), size=1, alpha=0.9) +
  geom_point(data=subset(df_drug_fraction_outliers, is_outlier==TRUE), 
             position=position_jitterdodge(jitter.width=0.50, jitter.height=0.50, seed=123),
             aes(color=tier3_variants), size=1, alpha=0.9) +
  ylab("Percent of variants with consensus drug predictions") + 
  xlab("") +
  theme(axis.text = element_text(size=28, color="black"),
        axis.title = element_text(size=21, color="black")) +
  theme(legend.justification = c(1, 1),
        legend.text = element_text(size=24, color="black"),
        legend.title = element_text(size=30, color="black")) +
  scale_fill_manual(name="Tier3 variants", values=tier3_colors) +
  scale_color_manual(name="Tier3 variants", values=tier3_colors)
# plot(p4a)
ggsave(outfile_p4a, p4a, width=10, height=8, units="in", dpi=1200)



## S4) PIE CHART TOTAL NUMBER OF VARIANTS WITH CONSENSUS DRUG SUPPORT ACROSS PATIENTS

# Use absolute variant numbers
# Sanity check that required input columns can be found
n_drug_columns = c("n_vars_drug_avail", "n_vars_drug_avail_no_tier3", "all_civic_variants", "all_civic_variants_no_tier3")
df_drug_snv_abs = prepare_data_for_plotting(data_snv, n_drug_columns, input_colnames, "all_civic_variants")
df_drug_cnv_abs = prepare_data_for_plotting(data_cnv, n_drug_columns, input_colnames, "all_civic_variants")

# Aggregate absolute numbers across patients in the cohort
fraction_drug_snv = aggregate_drug_data_across_patients(df_drug_snv_abs)
# Compute mean fraction for SNVs
mean_fraction_drug_snv = aggregate_drug_fraction_data_across_patients(df_drug_snv_fraction)

# Aggregate absolute numbers across patients in the cohort
fraction_drug_cnv = aggregate_drug_data_across_patients(df_drug_cnv_abs)
# Compute mean fraction for CNVs
mean_fraction_drug_cnv = aggregate_drug_fraction_data_across_patients(df_drug_cnv_fraction)

# Manually combine data for pie chart
df_drug_pie = data.frame(variant_type = rep(c("No drug info available", "With consensus drug support"), 4),
                         dataset = c(rep("SNVs", 2), rep("CNVs", 2), rep("SNVs", 2), rep("CNVs", 2)),
                         tier3_variants = c(rep("Including\ntier3", 4), rep("Excluding\ntier3", 4)),
                         value = c((100 - fraction_drug_snv$drug), fraction_drug_snv$drug, (100 - fraction_drug_cnv$drug), fraction_drug_cnv$drug,
                                   (100 - fraction_drug_snv$drug_no_tier3), fraction_drug_snv$drug_no_tier3, (100 - fraction_drug_cnv$drug_no_tier3), fraction_drug_cnv$drug_no_tier3))

# Assign percentages as labels for the pie chart
df_drug_pie$variant_type = factor(df_drug_pie$variant_type, levels=c("No drug info available", "With consensus drug support"))
df_drug_pie$dataset = factor(df_drug_pie$dataset, levels=c("SNVs", "CNVs"))
df_drug_pie$tier3_variants = factor(df_drug_pie$tier3_variants, levels=c("Including\ntier3", "Excluding\ntier3"))
df_drug_pie$label = percent(df_drug_pie$value/100, accuracy=0.1)

# Pie chart
# Text labels introduced post-generation
outfile_ps4 = paste0(opt$outfile_tag, ".s4.piechart_variants_with_drug_across_patients.pdf")
# outfile_ps4 = paste0(opt$outfile_tag, ".s4.piechart_variants_with_drug_across_patients.png")
ps4 = ggplot(df_drug_pie, aes(x="", y=value, fill=variant_type)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  facet_grid(tier3_variants ~ dataset, switch = "y") +
  coord_polar("y", start=0) +
  ylab(paste0("Total patients = ", nrow(df_civic_snv_fraction))) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size=22, hjust=0.5),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y.left = element_text(size=22, angle=0),
        panel.spacing = unit(1, "lines")) +
  scale_fill_manual(values=c("grey", "#0072B2"), name="")
# plot(ps4)
ggsave(outfile_ps4, ps4, width=10, height=8, units="in", dpi=1200)



## S5) BOX PLOT DISTRIBUTION OF MEAN CONSENSUS DRUG SUPPORT PREDICTIONS AVAILABLE PER VARIANT

# Use mean number of consensus predictions available per CIViC variant
df_n_consensus_snv_mean = prepare_mean_n_consensus_data_for_plotting(data_snv, input_colnames)
df_n_consensus_cnv_mean = prepare_mean_n_consensus_data_for_plotting(data_cnv, input_colnames)

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
df_n_consensus_mean =  combine_snv_and_cnv_data_for_plotting(df_n_consensus_snv_mean, df_n_consensus_cnv_mean, plotting_columns)

# Extract relevant features from the variable names in preparation for plotting
df_n_consensus_mean$tier3_variants = ifelse(grepl("_no_tier3", df_n_consensus_mean$variable), "Excluding", "Including")
df_n_consensus_mean$tier3_variants = factor(df_n_consensus_mean$tier3_variants, levels=c("Including", "Excluding"))

# Retrieve subset of outlier points for each group being plotted
n_consensus_mean_categories = names(table(df_n_consensus_mean$variable))
df_n_consensus_mean_outliers = get_outliers_per_group_and_dataset(df_n_consensus_mean, n_consensus_mean_categories)


# Combined version of box plot
# Showing SNVs and CNVs in X-axis (prediction numbers grouped by tier3 exclusion/inclusion)
outfile_ps5 = paste0(opt$outfile_tag, ".s5.boxplot_mean_n_consensus_predictions_per_variant.pdf")
# outfile_ps5 = paste0(opt$outfile_tag, ".s5.boxplot_mean_n_consensus_predictions_per_variant.png")
ps5 = ggplot(df_n_consensus_mean_outliers, aes(x=dataset, y=value, fill=tier3_variants)) +
  geom_boxplot(outlier.shape=NA, outlier.size=0, notch=FALSE) +
  geom_point(data=subset(df_n_consensus_mean_outliers, is_outlier==TRUE), 
             position=position_jitterdodge(jitter.width=0.50, jitter.height=0.50, seed=123),
             aes(color=tier3_variants), size=1, alpha=0.9) +  
  ylab("Mean number of consensus predictions per variant") + 
  xlab("") +
  theme(axis.text = element_text(size=22, color="black"),
        axis.title = element_text(size=22, color="black")) +
  theme(legend.justification = c(1, 1),
        legend.text = element_text(size=20, color="black"),
        legend.title = element_text(size=22, color="black")) +
  scale_fill_manual(name="Tier3 variants", values=tier3_colors) +
  scale_color_manual(name="Tier3 variants", values=tier3_colors)
# plot(ps5)
ggsave(outfile_ps5, ps5, width=10, height=8, units="in", dpi=1200)



## S6) BOX PLOT DISTRIBUTION OF CONSENSUS PREDICTIONS AVAILABLE PER VARIANT

# Use mean fraction of total consensus drug predictions per CIViC variant
# Sanity check that required input columns can be found
mean_fraction_consensus_columns = c("mean_percent_support", "mean_percent_resistance", "mean_percent_conflict", "mean_percent_unknown",
                                    "mean_percent_support_no_tier3", "mean_percent_resistance_no_tier3", "mean_percent_conflict_no_tier3", "mean_percent_unknown_no_tier3")
df_mean_fraction_consensus_snv = prepare_data_for_plotting(data_snv, mean_fraction_consensus_columns, input_colnames, "mean_percent_support")
df_mean_fraction_consensus_cnv = prepare_data_for_plotting(data_cnv, mean_fraction_consensus_columns, input_colnames, "mean_percent_support")

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
df_mean_fraction_consensus =  combine_snv_and_cnv_data_for_plotting(df_mean_fraction_consensus_snv, df_mean_fraction_consensus_cnv, plotting_columns)

# Extract relevant features from the variable names in preparation for plotting
df_mean_fraction_consensus$tier3_variants = ifelse(grepl("_no_tier3", df_mean_fraction_consensus$variable), "Excluding", "Including")
df_mean_fraction_consensus$tier3_variants = factor(df_mean_fraction_consensus$tier3_variants, levels=c("Including", "Excluding"))
df_mean_fraction_consensus$support_type = unlist(lapply(strsplit(as.character(df_mean_fraction_consensus$variable), "_"), `[`, 3))
# Sort types of consensus drug predictions available
consensus_prediction_order = c("support", "resistance", "conflict", "unknown")
df_mean_fraction_consensus$support_type = factor(df_mean_fraction_consensus$support_type, levels=consensus_prediction_order)

# Retrieve subset of outlier points for each group being plotted
mean_fraction_consensus_categories = names(table(df_mean_fraction_consensus$variable))
df_mean_fraction_consensus_outliers = get_outliers_per_group_and_dataset(df_mean_fraction_consensus, mean_fraction_consensus_categories)

# Include dummy outlier point from missing category to ensure the position and jittering of points work correctly
df_mean_fraction_consensus_outliers$is_outlier[5] = TRUE


# Combined version of box plot
# Showing SNVs and CNVs as facets (prediction numbers grouped by tier3 exclusion/inclusion)
outfile_ps6 = paste0(opt$outfile_tag, ".s6.boxplot_mean_fractions_consensus_predictions_per_variant.pdf")
# outfile_ps6 = paste0(opt$outfile_tag, ".s6.boxplot_mean_fractions_consensus_predictions_per_variant.png")
ps6 = ggplot(df_mean_fraction_consensus_outliers, aes(x=support_type, y=value, fill=tier3_variants)) +
  geom_boxplot(outlier.shape=NA, outlier.size=0, notch=FALSE) +
  # geom_point(position=position_jitterdodge(jitter.width=0.65, jitter.height=0.65, seed=123),
  #            aes(color=tier3_variants), size=1, alpha=0.9) +
  geom_point(data=subset(df_mean_fraction_consensus_outliers, is_outlier==TRUE), 
             position=position_jitterdodge(jitter.width=0.50, jitter.height=0.50, seed=123),
             aes(color=tier3_variants), size=1, alpha=0.9) +
  facet_grid(dataset ~ .) +
  ylab("Mean percent of consensus predictions per variant") +
  xlab("") +
  theme(axis.text = element_text(size=21, color="black"),
        axis.title = element_text(size=23, color="black")) +
  theme(legend.justification = c(1, 1),
        legend.text = element_text(size=20, color="black"),
        legend.title = element_text(size=22, color="black"),
        strip.text = element_text(size=20)) +
  scale_fill_manual(name="Tier3 variants", values=tier3_colors) +
  scale_color_manual(name="Tier3 variants", values=tier3_colors) +
  scale_x_discrete(breaks=c("support", "resistance", "conflict", "unknown"), labels=c("Support", "Resistance", "Conflict", "Unknown"))
# plot(ps6)
ggsave(outfile_ps6, ps6, width=10, height=8, units="in", dpi=1200)



## 4B) BOX PLOT DISTRIBUTION PERCENT OF PREDICTED CONSENSUS DRUGS OVERALL

# Use fractions of predicted consensus drugs overall
# Sanity check that required input columns can be found
overall_consensus_drugs_columns = c("percent_all_support_drugs", "percent_all_resistance_drugs", "percent_all_conflict_drugs", 
                                    "percent_all_unknown_drugs", "percent_mixed_drugs", "percent_all_support_drugs_no_tier3",
                                    "percent_all_resistance_drugs_no_tier3", "percent_all_conflict_drugs_no_tier3", 
                                    "percent_all_unknown_drugs_no_tier3", "percent_mixed_drugs_no_tier3")
df_overall_consensus_drugs_snv = prepare_data_for_plotting(data_snv, overall_consensus_drugs_columns, input_colnames, "percent_all_support_drugs")
df_overall_consensus_drugs_cnv = prepare_data_for_plotting(data_cnv, overall_consensus_drugs_columns, input_colnames, "percent_all_support_drugs")

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
df_overall_consensus_drugs =  combine_snv_and_cnv_data_for_plotting(df_overall_consensus_drugs_snv, df_overall_consensus_drugs_cnv, plotting_columns)

# Extract relevant features from the variable names in preparation for plotting
df_overall_consensus_drugs$tier3_variants = ifelse(grepl("_no_tier3", df_overall_consensus_drugs$variable), "Excluding", "Including")
df_overall_consensus_drugs$tier3_variants = factor(df_overall_consensus_drugs$tier3_variants, levels=c("Including", "Excluding"))
df_overall_consensus_drugs$support_type = unlist(lapply(strsplit(as.character(df_overall_consensus_drugs$variable), "_"), `[`, 3))
df_overall_consensus_drugs$support_type[df_overall_consensus_drugs$support_type == "drugs"] = "mixed"
# Sort types of consensus drugs according to the underlying predictions
consensus_drug_order = c("support", "resistance", "conflict", "unknown", "mixed")
df_overall_consensus_drugs$support_type = factor(df_overall_consensus_drugs$support_type, levels=consensus_drug_order)

# Retrieve subset of outlier points for each group being plotted
overall_consensus_drugs_categories = names(table(df_overall_consensus_drugs$variable))
df_overall_consensus_drugs_outliers = get_outliers_per_group_and_dataset(df_overall_consensus_drugs, overall_consensus_drugs_categories)

# Include dummy outlier point from missing category to ensure the position and jittering of points work correctly
df_overall_consensus_drugs_outliers$is_outlier[6] = TRUE


# Combined version of box plot
# Showing SNVs and CNVs as facets (drug numbers grouped by tier3 exclusion/inclusion)
outfile_p4b = paste0(opt$outfile_tag, ".4b.boxplot_fraction_overall_consensus_drugs.pdf")
# outfile_p4b = paste0(opt$outfile_tag, ".4b.boxplot_fraction_overall_consensus_drugs.png")
p4b = ggplot(df_overall_consensus_drugs_outliers, aes(x=support_type, y=value, fill=tier3_variants)) +
  geom_boxplot(outlier.shape=NA, outlier.size=0, notch=FALSE) +
  # geom_point(position=position_jitterdodge(jitter.width=0.65, jitter.height=0.65, seed=123),
  #            aes(color=tier3_variants), size=1, alpha=0.9) +
  geom_point(data=subset(df_overall_consensus_drugs_outliers, is_outlier==TRUE), 
             position=position_jitterdodge(jitter.width=0.50, jitter.height=0.50, seed=123),
             aes(color=tier3_variants), size=1, alpha=0.9) +
  facet_grid(dataset ~ .) +
  ylab("Percent of overall consensus drugs") + 
  xlab("") +
  theme(axis.text.x = element_text(size=28, color="black", angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=28, color="black"),
        axis.title = element_text(size=25, color="black")) +
  theme(legend.justification = c(1, 1),
        legend.text = element_text(size=24, color="black"),
        legend.title = element_text(size=30, color="black"),
        strip.text = element_text(size=24)) +
  scale_fill_manual(name="Tier3 variants", values=tier3_colors) +
  scale_color_manual(name="Tier3 variants", values=tier3_colors) +
  scale_x_discrete(breaks=c("support", "resistance", "conflict", "unknown", "mixed"),
                   labels=c("All-support", "All-resistance", "All-conflict", "All-unknown", "Mixed"))
# plot(p4b)
ggsave(outfile_p4b, p4b, width=10, height=8, units="in", dpi=1200)



########################################################################################################################


# NOT PART OF THE MANUSCRIPT


## S7A) DISTRIBUTION OF OVERALL CONSENSUS DRUG PREDICTIONS FOR CLASS "CT"

# Use fractions of predicted consensus drugs for class "ct" overall
# Sanity check that required input columns can be found
overall_consensus_drugs_ct_columns = c("percent_all_support_drugs_ct", "percent_all_resistance_drugs_ct", "percent_all_conflict_drugs_ct", "percent_all_unknown_drugs_ct", "percent_mixed_drugs_ct")
df_overall_consensus_drugs_ct_snv_sorted = prepare_data_for_plotting(data_snv, overall_consensus_drugs_ct_columns, input_colnames, "percent_all_support_drugs_ct")
df_overall_consensus_drugs_ct_cnv_sorted = prepare_data_for_plotting(data_cnv, overall_consensus_drugs_ct_columns, input_colnames, "percent_all_support_drugs_ct")
# overall_consensus_drugs_ct_columns = c("percent_all_support_drugs_ct_no_tier3", "percent_all_resistance_drugs_ct_no_tier3", "percent_all_conflict_drugs_ct_no_tier3", "percent_all_unknown_drugs_ct_no_tier3", "percent_mixed_drugs_ct_no_tier3")
# df_overall_consensus_drugs_ct_snv_sorted = prepare_data_for_plotting(data_snv, overall_consensus_drugs_ct_columns, input_colnames, "percent_all_support_drugs_ct_no_tier3")
# df_overall_consensus_drugs_ct_cnv_sorted = prepare_data_for_plotting(data_cnv, overall_consensus_drugs_ct_columns, input_colnames, "percent_all_support_drugs_ct_no_tier3")

# Re-sort all samples in the CNV dataset according to all types of consensus drugs predicted (simultaneously)
tmp_df_overall_consensus_drugs_ct_cnv = subset(data_cnv, select=c("sample_name", overall_consensus_drugs_ct_columns))
tmp_df_overall_consensus_drugs_ct_cnv_sorted = tmp_df_overall_consensus_drugs_ct_cnv[with(tmp_df_overall_consensus_drugs_ct_cnv, order(-percent_all_support_drugs_ct, -percent_all_resistance_drugs_ct, -percent_all_conflict_drugs_ct, -percent_all_unknown_drugs_ct, -percent_mixed_drugs_ct)),]
# tmp_df_overall_consensus_drugs_ct_cnv_sorted = tmp_df_overall_consensus_drugs_ct_cnv[with(tmp_df_overall_consensus_drugs_ct_cnv, order(-percent_all_support_drugs_ct_no_tier3, -percent_all_resistance_drugs_ct_no_tier3, -percent_all_conflict_drugs_ct_no_tier3, -percent_all_unknown_drugs_ct_no_tier3, -percent_mixed_drugs_ct_no_tier3)),]

# Keep track of the order of samples for later sorting the combined dataset
cnv_sample_order_ct = tmp_df_overall_consensus_drugs_ct_cnv_sorted$sample_name

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
df_overall_consensus_drugs_ct =  combine_snv_and_cnv_data_for_plotting(df_overall_consensus_drugs_ct_snv_sorted, df_overall_consensus_drugs_ct_cnv_sorted, plotting_columns)

# Re-sort all samples by the order of the CNV dataset
df_overall_consensus_drugs_ct$sample_order = as.character(df_overall_consensus_drugs_ct$sample_order)
df_overall_consensus_drugs_ct$sample_order = factor(df_overall_consensus_drugs_ct$sample_order, levels=cnv_sample_order_ct)

# Extract relevant features from the variable names in preparation for plotting
df_overall_consensus_drugs_ct$ct_class = unlist(lapply(strsplit(as.character(df_overall_consensus_drugs_ct$variable), "_"), tail, n = 1L))
df_overall_consensus_drugs_ct$support_type = unlist(lapply(strsplit(as.character(df_overall_consensus_drugs_ct$variable), "_"), `[`, 3))
df_overall_consensus_drugs_ct$support_type[df_overall_consensus_drugs_ct$support_type == "drugs"] = "mixed"
# Sort "ct" classes
ct_order = c("ct", "gt", "nct")
df_overall_consensus_drugs_ct$ct_class = factor(df_overall_consensus_drugs_ct$ct_class, levels=ct_order)
# Sort types of consensus drugs according to the underlying predictions
consensus_drug_order = c("support", "resistance", "conflict", "unknown", "mixed")
df_overall_consensus_drugs_ct$support_type = factor(df_overall_consensus_drugs_ct$support_type, levels=consensus_drug_order)

# Choose color scheme for plot
consensus_cols = c("support" = "limegreen", "resistance" = "red", "conflict" = "gold1", "unknown" = "navyblue", "mixed" = "lightblue2")
# consensus_cols = c("support" = "limegreen", "resistance" = "red", "conflict" = "orange", "unknown" = "navyblue", "mixed" = "lightblue2")
# consensus_cols = c("support" = "limegreen", "resistance" = "red", "conflict" = "orange", "unknown" = "gray55", "mixed" = "deepskyblue3")


# Combined version of barplot for class "ct" (stacked and not overlaying)
# Showing SNVs and CNVs as facets (sorting of sample based on CNV dataset)
outfile_overall_consensus_drugs_fraction_ct_bar = paste0(opt$outfile_tag, ".barplot_fraction_overall_consensus_drugs_ct.pdf")
# outfile_overall_consensus_drugs_fraction_ct_bar = paste0(opt$outfile_tag, ".barplot_fraction_overall_consensus_drugs_ct_no_tier3.pdf")
ps7a = ggplot(df_overall_consensus_drugs_ct, aes(x=sample_order, y=value, fill=support_type)) +
  geom_bar(stat='identity', position = "stack", width=1) +
  facet_grid(dataset ~ .) +
  ylab("Percent of overall consensus drugs") + 
  xlab("Patient") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=15, color="black")) +
  theme(legend.justification = c(1, 1),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=15, color="black"),
        strip.text=element_text(size=12)) +
  scale_fill_manual(name="Type of prediction", breaks=c("support", "resistance", "conflict", "unknown", "mixed"),
                    labels=c("All-support", "All-resistance", "All-conflict", "All-unknown", "Mixed"), values=consensus_cols)
# plot(ps7a)
ggsave(outfile_overall_consensus_drugs_fraction_ct_bar, ps7a, width=10, height=8, units="in", dpi=1200)



## S7B) DISTRIBUTION OF OVERALL CONSENSUS DRUG PREDICTIONS FOR CLASS "GT"

# Use fractions of predicted consensus drugs for class "gt" overall
# Sanity check that required input columns can be found
overall_consensus_drugs_gt_columns = c("percent_all_support_drugs_gt", "percent_all_resistance_drugs_gt", "percent_all_conflict_drugs_gt", "percent_all_unknown_drugs_gt", "percent_mixed_drugs_gt")
df_overall_consensus_drugs_gt_snv_sorted = prepare_data_for_plotting(data_snv, overall_consensus_drugs_gt_columns, input_colnames, "percent_all_support_drugs_gt")
df_overall_consensus_drugs_gt_cnv_sorted = prepare_data_for_plotting(data_cnv, overall_consensus_drugs_gt_columns, input_colnames, "percent_all_support_drugs_gt")
# overall_consensus_drugs_gt_columns = c("percent_all_support_drugs_gt_no_tier3", "percent_all_resistance_drugs_gt_no_tier3", "percent_all_conflict_drugs_gt_no_tier3", "percent_all_unknown_drugs_gt_no_tier3", "percent_mixed_drugs_gt_no_tier3")
# df_overall_consensus_drugs_gt_snv_sorted = prepare_data_for_plotting(data_snv, overall_consensus_drugs_gt_columns, input_colnames, "percent_all_support_drugs_gt_no_tier3")
# df_overall_consensus_drugs_gt_cnv_sorted = prepare_data_for_plotting(data_cnv, overall_consensus_drugs_gt_columns, input_colnames, "percent_all_support_drugs_gt_no_tier3")

# Re-sort all samples in the CNV dataset according to all types of consensus drugs predicted (simultaneously)
tmp_df_overall_consensus_drugs_gt_cnv = subset(data_cnv, select=c("sample_name", overall_consensus_drugs_gt_columns))
tmp_df_overall_consensus_drugs_gt_cnv_sorted = tmp_df_overall_consensus_drugs_gt_cnv[with(tmp_df_overall_consensus_drugs_gt_cnv, order(-percent_all_support_drugs_gt, -percent_all_resistance_drugs_gt, -percent_all_conflict_drugs_gt, -percent_all_unknown_drugs_gt, -percent_mixed_drugs_gt)),]
# tmp_df_overall_consensus_drugs_gt_cnv_sorted = tmp_df_overall_consensus_drugs_gt_cnv[with(tmp_df_overall_consensus_drugs_gt_cnv, order(-percent_all_support_drugs_gt_no_tier3, -percent_all_resistance_drugs_gt_no_tier3, -percent_all_conflict_drugs_gt_no_tier3, -percent_all_unknown_drugs_gt_no_tier3, -percent_mixed_drugs_gt_no_tier3)),]

# Keep track of the order of samples for later sorting the combined dataset
cnv_sample_order_gt = tmp_df_overall_consensus_drugs_gt_cnv_sorted$sample_name

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
df_overall_consensus_drugs_gt =  combine_snv_and_cnv_data_for_plotting(df_overall_consensus_drugs_gt_snv_sorted, df_overall_consensus_drugs_gt_cnv_sorted, plotting_columns)

# Re-sort all samples by the order of the CNV dataset
df_overall_consensus_drugs_gt$sample_order = as.character(df_overall_consensus_drugs_gt$sample_order)
df_overall_consensus_drugs_gt$sample_order = factor(df_overall_consensus_drugs_gt$sample_order, levels=cnv_sample_order_gt)

# Extract relevant features from the variable names in preparation for plotting
df_overall_consensus_drugs_gt$ct_class = unlist(lapply(strsplit(as.character(df_overall_consensus_drugs_gt$variable), "_"), tail, n = 1L))
df_overall_consensus_drugs_gt$support_type = unlist(lapply(strsplit(as.character(df_overall_consensus_drugs_gt$variable), "_"), `[`, 3))
df_overall_consensus_drugs_gt$support_type[df_overall_consensus_drugs_gt$support_type == "drugs"] = "mixed"
# Sort "ct" classes
df_overall_consensus_drugs_gt$ct_class = factor(df_overall_consensus_drugs_gt$ct_class, levels=ct_order)
# Sort types of consensus drugs according to the underlying predictions
df_overall_consensus_drugs_gt$support_type = factor(df_overall_consensus_drugs_gt$support_type, levels=consensus_drug_order)


# Combined version of barplot for class "gt" (stacked and not overlaying)
# Showing SNVs and CNVs as facets (sorting of sample based on CNV dataset)
outfile_overall_consensus_drugs_fraction_gt_bar = paste0(opt$outfile_tag, ".barplot_fraction_overall_consensus_drugs_gt.pdf")
# outfile_overall_consensus_drugs_fraction_gt_bar = paste0(opt$outfile_tag, ".barplot_fraction_overall_consensus_drugs_gt_no_tier3.pdf")
ps7b = ggplot(df_overall_consensus_drugs_gt, aes(x=sample_order, y=value, fill=support_type)) +
  geom_bar(stat='identity', position = "stack", width=1) +
  facet_grid(dataset ~ .) +
  ylab("Percent of overall consensus drugs") + 
  xlab("Patient") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=15, color="black")) +
  theme(legend.justification = c(1, 1),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=15, color="black"),
        strip.text=element_text(size=12)) +
  scale_fill_manual(name="Type of prediction", breaks=c("support", "resistance", "conflict", "unknown", "mixed"),
                    labels=c("All-support", "All-resistance", "All-conflict", "All-unknown", "Mixed"), values=consensus_cols)
# plot(ps7b)
ggsave(outfile_overall_consensus_drugs_fraction_gt_bar, ps7b, width=10, height=8, units="in", dpi=1200)



## S7C) DISTRIBUTION OF OVERALL CONSENSUS DRUG PREDICTIONS FOR CLASS "NCT"

# Use fractions of predicted consensus drugs for class "nct" overall
# Sanity check that required input columns can be found
overall_consensus_drugs_nct_columns = c("percent_all_support_drugs_nct", "percent_all_resistance_drugs_nct", "percent_all_conflict_drugs_nct", "percent_all_unknown_drugs_nct", "percent_mixed_drugs_nct")
df_overall_consensus_drugs_nct_snv_sorted = prepare_data_for_plotting(data_snv, overall_consensus_drugs_nct_columns, input_colnames, "percent_all_support_drugs_nct")
df_overall_consensus_drugs_nct_cnv_sorted = prepare_data_for_plotting(data_cnv, overall_consensus_drugs_nct_columns, input_colnames, "percent_all_support_drugs_nct")
# overall_consensus_drugs_nct_columns = c("percent_all_support_drugs_nct_no_tier3", "percent_all_resistance_drugs_nct_no_tier3", "percent_all_conflict_drugs_nct_no_tier3", "percent_all_unknown_drugs_nct_no_tier3", "percent_mixed_drugs_nct_no_tier3")
# df_overall_consensus_drugs_nct_snv_sorted = prepare_data_for_plotting(data_snv, overall_consensus_drugs_nct_columns, input_colnames, "percent_all_support_drugs_nct_no_tier3")
# df_overall_consensus_drugs_nct_cnv_sorted = prepare_data_for_plotting(data_cnv, overall_consensus_drugs_nct_columns, input_colnames, "percent_all_support_drugs_nct_no_tier3")

# Re-sort all samples in the SNV dataset according to all types of consensus drugs predicted (simultaneously)
tmp_df_overall_consensus_drugs_nct_snv = subset(data_snv, select=c("sample_name", overall_consensus_drugs_nct_columns))
tmp_df_overall_consensus_drugs_nct_snv_sorted = tmp_df_overall_consensus_drugs_nct_snv[with(tmp_df_overall_consensus_drugs_nct_snv, order(-percent_all_support_drugs_nct, -percent_all_resistance_drugs_nct, -percent_all_conflict_drugs_nct, -percent_all_unknown_drugs_nct, -percent_mixed_drugs_nct)),]
# tmp_df_overall_consensus_drugs_nct_snv_sorted = tmp_df_overall_consensus_drugs_nct_snv[with(tmp_df_overall_consensus_drugs_nct_snv, order(-percent_all_support_drugs_nct_no_tier3, -percent_all_resistance_drugs_nct_no_tier3, -percent_all_conflict_drugs_nct_no_tier3, -percent_all_unknown_drugs_nct_no_tier3, -percent_mixed_drugs_nct_no_tier3)),]

# Keep track of the order of samples for later sorting the combined dataset
snv_sample_order_nct = tmp_df_overall_consensus_drugs_nct_snv_sorted$sample_name

# Combine SNV and CNV datasets
# Sanity check that required input columns can be found
df_overall_consensus_drugs_nct =  combine_snv_and_cnv_data_for_plotting(df_overall_consensus_drugs_nct_snv_sorted, df_overall_consensus_drugs_nct_cnv_sorted, plotting_columns)

# Re-sort all samples by the order of the SNV dataset
df_overall_consensus_drugs_nct$sample_order = as.character(df_overall_consensus_drugs_nct$sample_order)
df_overall_consensus_drugs_nct$sample_order = factor(df_overall_consensus_drugs_nct$sample_order, levels=snv_sample_order_nct)

# Extract relevant features from the variable names in preparation for plotting
df_overall_consensus_drugs_nct$ct_class = unlist(lapply(strsplit(as.character(df_overall_consensus_drugs_nct$variable), "_"), tail, n = 1L))
df_overall_consensus_drugs_nct$support_type = unlist(lapply(strsplit(as.character(df_overall_consensus_drugs_nct$variable), "_"), `[`, 3))
df_overall_consensus_drugs_nct$support_type[df_overall_consensus_drugs_nct$support_type == "drugs"] = "mixed"
# Sort "ct" classes
df_overall_consensus_drugs_nct$ct_class = factor(df_overall_consensus_drugs_nct$ct_class, levels=ct_order)
# Sort types of consensus drugs according to the underlying predictions
df_overall_consensus_drugs_nct$support_type = factor(df_overall_consensus_drugs_nct$support_type, levels=consensus_drug_order)


# Combined version of barplot for class "gt" (stacked and not overlaying)
# Showing SNVs and CNVs as facets (sorting of sample based on SNV dataset)
outfile_overall_consensus_drugs_fraction_nct_bar = paste0(opt$outfile_tag, ".barplot_fraction_overall_consensus_drugs_nct.pdf")
# outfile_overall_consensus_drugs_fraction_nct_bar = paste0(opt$outfile_tag, ".barplot_fraction_overall_consensus_drugs_nct_no_tier3.pdf")
ps7c = ggplot(df_overall_consensus_drugs_nct, aes(x=sample_order, y=value, fill=support_type)) +
  geom_bar(stat='identity', position = "stack", width=1) +
  facet_grid(dataset ~ .) +
  ylab("Percent of overall consensus drugs") + 
  xlab("Patient") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=15, color="black")) +
  theme(legend.justification = c(1, 1),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=15, color="black"),
        strip.text=element_text(size=12)) +
  scale_fill_manual(name="Type of prediction", breaks=c("support", "resistance", "conflict", "unknown", "mixed"),
                    labels=c("All-support", "All-resistance", "All-conflict", "All-unknown", "Mixed"), values=consensus_cols)
# plot(ps7c)
ggsave(outfile_overall_consensus_drugs_fraction_nct_bar, ps7c, width=10, height=8, units="in", dpi=1200)
