#!/usr/bin/env Rscript
################################################################################
## Generate figures for CIViCutils manuscript
## Author: Lourdes Rosano
## Date created: April 2022
## R Version: 4.1.2
################################################################################


library(optparse)
library(scales)
library(reshape2)
library(plyr)
library(ggplot2)


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

# 
prepare_data_for_plotting = function(dataset, column_names, header, sort_variable){
  check_column_names(column_names, header)
  column_names = c("sample_name", column_names)
  subdata = subset(dataset, select=column_names)
  df = melt(subdata, id="sample_name")
  df_sorted = arrange(df, sample_name) #, increasing = TRUE)
  df_sorted = ddply(df_sorted, "sample_name", transform, label_ypos = value)
  # Sort samples by provided variable (decreasing order)
  sample_order = df_sorted[df_sorted$variable == sort_variable,,][order(df_sorted$value[df_sorted$variable == sort_variable], decreasing=T), "sample_name"]
  df_sorted$sample_order = factor(df_sorted$sample_name, levels=sample_order)
  return(df_sorted)
}


#
prepare_fraction_civic_data_for_plotting = function(dataset, header){
  n_civic_columns = c("sample_name", "all_variants", "all_civic_variants")
  check_column_names(n_civic_columns, header)
  subdata = subset(dataset, select=n_civic_columns)
  # Avoid having NaN in the results due to patients having all_variants=0
  # subdata$fraction = subdata$all_civic_variants / subdata$all_variants
  skip_rows = subdata$all_variants!=0
  subdata$fraction = 0
  # subdata[skip_rows,]$fraction = subdata[skip_rows,]$all_civic_variants / subdata[skip_rows,]$all_variants
  subdata[skip_rows,]$fraction = (subdata[skip_rows,]$all_civic_variants / subdata[skip_rows,]$all_variants)*100
  df = subset(subdata, select=c("sample_name", "fraction"))
  # Sort samples by computed fraction (decreasing order)
  sample_order = df[order(df$fraction, decreasing=T), "sample_name"]
  df$sample_order = factor(df$sample_name, levels=sample_order)
  return(df)
}


# 
aggregate_civic_across_patients = function(df_civic){
  # Sum all variants and CIViC variants across patients
  all_variants = sum(subset(df_civic, subset= variable == "all_variants", select = value))
  all_civic_variants = sum(subset(df_civic, subset= variable == "all_civic_variants", select = value))
  # Compute total fraction of variants with and without CIViC info across entire patient cohort
  fraction_civic_variants = (all_civic_variants / all_variants) * 100
  return(fraction_civic_variants)
}


#
aggregate_tier_across_patients = function(df_tier){
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


# 
aggregate_tier_fraction_across_patients = function(df_tier){
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

# 
prepare_fraction_tier_data_for_plotting = function(dataset, header){
  n_tier_columns = c("sample_name", "all_civic_variants", "n_tier_1", "n_tier_1b", "n_tier_2", "n_tier_3")
  # n_tier_columns = c("sample_name", "all_civic_variants", "n_tier_1", "n_tier_1b", "n_tier_2", "n_tier_3")
  check_column_names(n_tier_columns, header)
  subdata = subset(dataset, select=n_tier_columns)
  # Avoid having NaN in the results due to patients having all_variants=0
  skip_rows = subdata$all_civic_variants!=0
  subdata$fraction_tier_1 = 0
  subdata$fraction_tier_1b = 0
  subdata$fraction_tier_2 = 0
  subdata$fraction_tier_3 = 0
  subdata[skip_rows,]$fraction_tier_1 = (subdata[skip_rows,]$n_tier_1 / subdata[skip_rows,]$all_civic_variants)*100
  subdata[skip_rows,]$fraction_tier_1b = (subdata[skip_rows,]$n_tier_1b / subdata[skip_rows,]$all_civic_variants)*100
  subdata[skip_rows,]$fraction_tier_2 = (subdata[skip_rows,]$n_tier_2 / subdata[skip_rows,]$all_civic_variants)*100
  subdata[skip_rows,]$fraction_tier_3 = (subdata[skip_rows,]$n_tier_3 / subdata[skip_rows,]$all_civic_variants)*100
  # new_subdata = subset(subdata, select=c("sample_name", "all_civic_variants", "fraction_tier_1", "fraction_tier_1b", "fraction_tier_2", "fraction_tier_3"))
  new_subdata = subset(subdata, select=c("sample_name", "fraction_tier_1", "fraction_tier_1b", "fraction_tier_2", "fraction_tier_3"))
  new_subdata_sorted = new_subdata[with(new_subdata, order(fraction_tier_1, fraction_tier_1b, fraction_tier_2, fraction_tier_3)),]
  # new_subdata_sorted = new_subdata[with(new_subdata, order(fraction_tier_3, fraction_tier_2, fraction_tier_1b, fraction_tier_1)),]
  new_subdata_sorted = new_subdata[with(new_subdata, order(fraction_tier_1, fraction_tier_1b, fraction_tier_2, fraction_tier_3)),]
  sample_order = new_subdata_sorted$sample_name
  
  # Reformat dataset in preparation for plotting
  df = melt(new_subdata, id="sample_name")
  df_sorted = arrange(df, sample_name) #, increasing = TRUE)
  df_sorted = ddply(df_sorted, "sample_name", transform, label_ypos = value)
  
  # Sort samples by provided variable (decreasing order)
  # sample_order = df_sorted[df_sorted$variable == "fraction_tier_3",,][order(df_sorted$value[df_sorted$variable == "fraction_tier_3"], decreasing=T), "sample_name"]
  # sample_order = df_sorted[df_sorted$variable == "all_civic_variants",,][order(df_sorted$value[df_sorted$variable == "all_civic_variants"], decreasing=T), "sample_name"]
  df_sorted$sample_order = factor(df_sorted$sample_name, levels=sample_order)
  # Remove number of variants matched in CIViC (used only for sorting samples)
  # df_sorted = df_sorted[df_sorted$variable != "all_civic_variants",,drop=F]
  # df_sorted$variable = as.character(df_sorted$variable)
  # Sort tiers by increasing order of priority
  tier_order = c("fraction_tier_3", "fraction_tier_2", "fraction_tier_1b", "fraction_tier_1")
  df_sorted$tier_order = factor(df_sorted$variable, levels=tier_order)
  return(df_sorted)
}


# 
combine_snv_and_cnv_for_plotting = function(dataset_snv, dataset_cnv, column_names){
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


opt = list(
  infile_snv = "/Users/rlourdes/Documents/PROJECTS/CIVIC/tests/test_complete_tcga_blca.snvs.tsv",
  infile_cnv = "/Users/rlourdes/Documents/PROJECTS/CIVIC/tests/test_complete_tcga_blca.cnvs.tsv",
  outfile_tag = "/Users/rlourdes/Documents/PROJECTS/CIVIC/tests/test_figures"
)

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


## A1) NUMBER OF VARIANTS MATCHED IN CIVIC

## Plot version showing absolute variant numbers

# Sanity check that required input columns can be found
n_civic_columns = c("all_variants", "all_civic_variants")
df_civic_snv_sorted = prepare_data_for_plotting(data_snv, n_civic_columns, input_colnames, "all_variants")
df_civic_cnv_sorted = prepare_data_for_plotting(data_cnv, n_civic_columns, input_colnames, "all_variants")

# Choose color scheme for plot
# civic_colors = c("all_variants" = "grey", "all_civic_variants" = "#0072B2")

# NOTE: barplots shown below are overlaying and not stacked!

# outfile_civic_snv = paste0(opt$outfile_tag, ".snvs_in_civic.pdf")
# p1a = ggplot(df_civic_snv_sorted, aes(x=sample_order, y=value, fill=variable)) +
#   geom_bar(stat='identity', position = "identity", width=.7) +
#   ylab("Number of SNVs") + 
#   xlab("Patient") +
#   # ggtitle("") +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#   theme(legend.justification = c(1, 1),legend.position = c(1,1)) +
#   scale_y_continuous(breaks=seq(0, 600, by = 100)) +
#   scale_fill_manual(name="SNVs", breaks=c("all_variants", "all_civic_variants"), labels=c("all", "matched in CIViC"), values=civic_colors)
# plot(p1a)
# ggsave(outfile_civic_snv, p1a, width=10, height=8, units="in", dpi=1200)

# outfile_civic_cnv = paste0(opt$outfile_tag, ".cnvs_in_civic.pdf")
# p1b = ggplot(df_civic_cnv_sorted, aes(x=sample_order, y=value, fill=variable)) +
#   geom_bar(stat='identity', position = "identity", width=.7) +
#   ylab("Number of CNVs") + 
#   xlab("Patient") +
#   # ggtitle("") +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#   theme(legend.justification = c(1, 1),legend.position = c(1,1)) +
#   # coord_cartesian(ylim=c(0, 200)) + 
#   scale_fill_manual(name="CNVs", breaks=c("all_variants", "all_civic_variants"), labels=c("all", "matched in CIViC"), values=civic_colors)
# plot(p1b)
# ggsave(outfile_civic_cnv, p1b, width=10, height=8, units="in", dpi=1200)


## Plot version showing fraction of total variants having CIViCutils results

# Sanity check that required input columns can be found
df_civic_snv_fraction = prepare_fraction_civic_data_for_plotting(data_snv, input_colnames)
df_civic_cnv_fraction = prepare_fraction_civic_data_for_plotting(data_cnv, input_colnames)

# outfile_civic_snv = paste0(opt$outfile_tag, ".fraction_snvs_in_civic.pdf")
# p1a = ggplot(df_civic_snv_fraction, aes(x=sample_order, y=fraction, fill="#0072B2")) +
#   geom_bar(stat='identity', position = "identity", width=.7) +
#   ylab("Percent of SNVs") + 
#   xlab("Patient") +
#   # ggtitle("") +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#   theme(legend.justification = c(1, 1),legend.position = c(1,1)) +
#   scale_fill_manual(name="SNVs", labels=c("Matched in CIViC"), values=c("#0072B2"))
# plot(p1a)
# ggsave(outfile_civic_snv, p1a, width=10, height=8, units="in", dpi=1200)

# outfile_civic_cnv = paste0(opt$outfile_tag, ".fraction_cnvs_in_civic.pdf")
# p1b = ggplot(df_civic_cnv_fraction, aes(x=sample_order, y=fraction, fill="#0072B2")) +
#   geom_bar(stat='identity', position = "identity", width=.7) +
#   ylab("Percent of CNVs") + 
#   xlab("Patient") +
#   # ggtitle("") +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#   theme(legend.justification = c(1, 1),legend.position = c(1,1)) +
#   scale_fill_manual(name="CNVs", labels=c("Matched in CIViC"), values=c("#0072B2"))
# plot(p1b)
# ggsave(outfile_civic_cnv, p1b, width=10, height=8, units="in", dpi=1200)


## Plot combined version of this plot, showing SNVs and CNVs as facets (same axis, use sample sorting from SNV dataset)

# Sanity check that required input columns can be found
plotting_columns = c("sample_name", "fraction", "sample_order")
df_civic_fraction =  combine_snv_and_cnv_for_plotting(df_civic_snv_fraction, df_civic_cnv_fraction, plotting_columns)
# plotting_columns = c("sample_name", "variable", "value", "label_ypos", "sample_order")
# df_civic_sorted =  combine_snv_and_cnv_for_plotting(df_civic_snv_sorted, df_civic_cnv_sorted, plotting_columns)

# outfile_civic = paste0(opt$outfile_tag, ".variants_in_civic.pdf")
# p1c = ggplot(df_civic_sorted, aes(x=sample_order, y=value, fill=variable)) +
#   geom_bar(stat='identity', position = "identity", width=.7) +
#   ylab("Number of variants") + 
#   xlab("Patient") +
#   # ggtitle("") +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#   theme(legend.justification = c(1, 1),legend.position = c(1,1)) +
#   facet_grid(dataset ~ .) +
#   # facet_grid(dataset ~ ., switch="both") +
#   # coord_cartesian(ylim=c(0, 200)) + 
#   scale_fill_manual(name="Variants", breaks=c("all_variants", "all_civic_variants"), labels=c("all", "matched in CIViC"), values=civic_colors)
# plot(p1c)
# ggsave(outfile_civic, p1c, width=10, height=8, units="in", dpi=1200)

outfile_civic = paste0(opt$outfile_tag, ".fraction_variants_in_civic.pdf")
p1c = ggplot(df_civic_fraction, aes(x=sample_order, y=fraction, fill="#0072B2")) +
  geom_bar(stat='identity', position = "identity", width=.7) +
  ylab("Percent of variants") + 
  xlab("Patient") +
  # ggtitle("") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(legend.justification = c(1, 1),legend.position = c(1,1)) +
  facet_grid(dataset ~ .) +
  scale_fill_manual(name="Variants", labels=c("Matched in CIViC"), values=c("#0072B2"))
# plot(p1c)
ggsave(outfile_civic, p1c, width=10, height=8, units="in", dpi=1200)

# TODO: look into "cutting" a section of the y-axis, to avoid having a very broad range just to show tier4 variants (have a "gap" instead, and show only relevant part of the barplot + max value of the y-axis)



## A2) TOTAL NUMBER OF VARIANTS MATCHED IN CIVIC ACROSS PATIENTS

fraction_civic_snv = aggregate_civic_across_patients(df_civic_snv_sorted)
mean_fraction_civic_snv = (sum(df_civic_snv_fraction$fraction) / nrow(df_civic_snv_fraction))

fraction_civic_cnv = aggregate_civic_across_patients(df_civic_cnv_sorted)
mean_fraction_civic_cnv = (sum(df_civic_cnv_fraction$fraction) / nrow(df_civic_cnv_fraction))

df_civic_pie = data.frame(variant_type = rep(c("No match", "Matched in CIViC"), 2),
                          dataset = c(rep("SNVs", 2), rep("CNVs", 2)),
                          value = c((100 - fraction_civic_snv), fraction_civic_snv, (100 - fraction_civic_cnv), fraction_civic_cnv))

# Assign percentages as labels for the piechart
df_civic_pie$variant_type = factor(df_civic_pie$variant_type, levels=c("No match", "Matched in CIViC"))
df_civic_pie$dataset = factor(df_civic_pie$dataset, levels=c("SNVs", "CNVs"))
df_civic_pie$label = percent(df_civic_pie$value/100, accuracy=1)

outfile_civic_pie = paste0(opt$outfile_tag, ".pie_chart_variants_in_civic_across_patients.pdf")
p1 = ggplot(df_civic_pie, aes(x="", y=value, fill=variant_type)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start=0) +
  # "#CC79A7","#F0E442",
  # scale_fill_manual(values = c("#999999", "#0072B2"), name="") +
  scale_fill_manual(values = c("grey", "#0072B2"), name="") +
  geom_text(aes(label = label, x = 1.05), position = position_stack(vjust = 0.5), size=10) +
  ylab(paste0("Total patients = ", nrow(df_civic_snv_fraction))) +
  theme_minimal() +
  facet_grid(. ~ dataset) +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_text(size=40, hjust = 0.5),#, face="bold"),
        axis.title.y = element_blank(),
        panel.grid=element_blank(),
        legend.text=element_text(size=30),
        strip.text = element_text(size = 35))
# plot(p1)
ggsave(outfile_civic_pie, p1, width=10, height=8, units="in", dpi=1200)

# TODO: improve centering of labels (especially 21% of SNVs is off)



## B1) DISTRIBUTION OF TIER CLASSIFICATION (only for variants matched in CIVIC)

## Plot version showing absolute variant numbers

# Sanity check that required input columns can be found
n_tier_columns = c("n_tier_1", "n_tier_1b", "n_tier_2", "n_tier_3", "all_civic_variants")
df_tier_snv_sorted = prepare_data_for_plotting(data_snv, n_tier_columns, input_colnames, "all_civic_variants")
df_tier_cnv_sorted = prepare_data_for_plotting(data_cnv, n_tier_columns, input_colnames, "all_civic_variants")

# Sort tiers by increasing order of priority
tier_order = c("n_tier_3", "n_tier_2", "n_tier_1b", "n_tier_1")
df_tier_snv_sorted$tier_order = factor(df_tier_snv_sorted$variable, levels=tier_order)
df_tier_cnv_sorted$tier_order = factor(df_tier_cnv_sorted$variable, levels=tier_order)

# Remove number of variants matched in CIViC (used only for sorting samples)
df_tier_snv_sorted = df_tier_snv_sorted[df_tier_snv_sorted$variable != "all_civic_variants",,drop=F]
df_tier_cnv_sorted = df_tier_cnv_sorted[df_tier_cnv_sorted$variable != "all_civic_variants",,drop=F]

# Choose color scheme for plot
# #56B4E9
# #E69F00
tier_cols = c("n_tier_1" = "red", "n_tier_1b" = "#0072B2", "n_tier_2" = "gold", "n_tier_3" = "skyblue")
# tier_cols = c("n_tier_1" = "red", "n_tier_1b" = "#0072B2", "n_tier_2" = "gold", "n_tier_3" = "skyblue", "tier_4" = "grey")

# NOTE: barplots shown below are stacked and not overlaying!

# outfile_tier_snv = paste0(opt$outfile_tag, ".snvs_tier_distribution.pdf")
# p2a = ggplot(df_tier_snv_sorted, aes(x=sample_order, y=value, fill=tier_order)) +
#   geom_bar(stat='identity', position = "stack", width=.7) +
#   ylab("Number of SNVs") + 
#   xlab("Patient") +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#   theme(legend.justification = c(1, 1),legend.position = c(1,1)) +
#   # coord_cartesian(ylim=c(0, 30)) +
#   # scale_y_continuous(breaks=seq(0, 70, by = 10)) +
#   scale_fill_manual(name="Tier class", breaks=c("n_tier_1", "n_tier_1b", "n_tier_2", "n_tier_3"), labels=c("tier1", "tier1b", "tier2", "tier3"), values=tier_cols)
#   # scale_fill_manual(name="Tier class", breaks=c("n_tier_1", "n_tier_1b", "n_tier_2", "n_tier_3", "n_tier4"), labels=c("tier1", "tier1b", "tier2", "tier3", "tier4"), values=tier_cols)
# plot(p2a)
# ggsave(outfile_tier_snv, p2a, width=10, height=8, units="in", dpi=1200)

# outfile_tier_cnv = paste0(opt$outfile_tag, ".cnvs_tier_distribution.pdf")
# p2b = ggplot(df_tier_cnv_sorted, aes(x=sample_order, y=value, fill=tier_order)) +
#   geom_bar(stat='identity', position = "stack", width=.7) +
#   ylab("Number of CNVs") + 
#   xlab("Patient") +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#   theme(legend.justification = c(1, 1),legend.position = c(1,1)) +
#   # coord_cartesian(ylim=c(0, 30)) + 
#   scale_fill_manual(name="Tier class", breaks=c("n_tier_1", "n_tier_1b", "n_tier_2", "n_tier_3"), labels=c("tier1", "tier1b", "tier2", "tier3"), values=tier_cols)
# # scale_fill_manual(name="Tier class", breaks=c("n_tier_1", "n_tier_1b", "n_tier_2", "n_tier_3", "n_tier4"), labels=c("tier1", "tier1b", "tier2", "tier3", "tier4"), values=tier_cols)
# plot(p2b)
# ggsave(outfile_tier_cnv, p2b, width=10, height=8, units="in", dpi=1200)


## Plot version showing fraction of total variants having CIViCutils results

df_tier_snv_fraction = prepare_fraction_tier_data_for_plotting(data_snv, input_colnames)
df_tier_cnv_fraction = prepare_fraction_tier_data_for_plotting(data_cnv, input_colnames)


## Plot combined version of this plot, showing SNVs and CNVs as facets (same axis, use sample sorting from SNV dataset)

# Sanity check that required input columns can be found
plotting_columns_tier = c("sample_name", "variable", "value", "label_ypos", "sample_order", "tier_order")
df_tier_fraction =  combine_snv_and_cnv_for_plotting(df_tier_snv_fraction, df_tier_cnv_fraction, plotting_columns_tier)

# Choose color scheme for plot
tier_cols = c("fraction_tier_1" = "red", "fraction_tier_1b" = "#0072B2", "fraction_tier_2" = "gold", "fraction_tier_3" = "skyblue")

# NOTE: barplot shown below is stacked and not overlaying!

outfile_tier = paste0(opt$outfile_tag, ".fraction_variants_tier_distribution.pdf")
p2c = ggplot(df_tier_fraction, aes(x=sample_order, y=value, fill=tier_order)) +
  geom_bar(stat='identity', position = "stack", width=.7) +
  ylab("Percent of variants matched in CIViC") + 
  xlab("Patient") +
  # ggtitle("") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  # theme(legend.justification = c(1, 1), legend.position = c(1,1)) +
  theme(legend.justification = c(1, 1)) +
  facet_grid(dataset ~ .) +
  # facet_grid(dataset ~ ., switch="both") +
  scale_fill_manual(name="Tier class", breaks=c("fraction_tier_1", "fraction_tier_1b", "fraction_tier_2", "fraction_tier_3"), labels=c("tier1", "tier1b", "tier2", "tier3"), values=tier_cols)
# plot(p2c)
ggsave(outfile_tier, p2c, width=10, height=8, units="in", dpi=1200)


## B2) TOTAL DISTRIBUTION OF TIER CLASSIFICATION ACROSS PATIENTS

fractions_per_tier_snv = aggregate_tier_across_patients(df_tier_snv_sorted)
mean_fractions_per_tier_snv = aggregate_tier_fraction_across_patients(df_tier_snv_fraction)

fractions_per_tier_cnv = aggregate_tier_across_patients(df_tier_cnv_sorted)
mean_fractions_per_tier_cnv = aggregate_tier_fraction_across_patients(df_tier_cnv_fraction)

df_tier_pie = data.frame(tier = rep(c("Tier 1", "Tier 1b", "Tier 2", "Tier 3"), 2),
                         dataset = c(rep("SNVs", 4), rep("CNVs", 4)),
                         value = c(unlist(fractions_per_tier_snv), unlist(fractions_per_tier_cnv)))

# Assign percentages as labels for the piechart
df_tier_pie$dataset = factor(df_tier_pie$dataset, levels=c("SNVs", "CNVs"))
df_tier_pie$label = percent(df_tier_pie$value/100, accuracy=1)
# Avoid showing 0% classes in the pie chart
# df_tier_pie$tier = factor(df_tier_pie$tier, levels=c("Tier 1", "Tier 1b", "Tier 2", "Tier 3"))

# Remove entries of 0% variants to avoid showing them in the pie chart
sub_df_tier_pie = subset(df_tier_pie, subset= label != "0%")
# subset(df_tier_pie, subset= label == "0%")

# Choose color scheme for plot
tier_cols = c("Tier 1" = "red", "Tier 1b" = "#0072B2", "Tier 2" = "gold", "Tier 3" = "skyblue")

outfile_tier_pie = paste0(opt$outfile_tag, ".pie_chart_tiers_across_patients.pdf")
p2 = ggplot(sub_df_tier_pie, aes(x="", y=value, fill=tier)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = tier_cols, name="") +
  geom_text(aes(label = label, x = 1.1), position = position_stack(vjust = 0.5), size=7) +
  ylab(paste0("Total patients = ", nrow(df_civic_snv_fraction))) +
  theme_minimal() +
  facet_grid(. ~ dataset) +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_text(size=40, hjust = 0.5),#, face="bold"),
        axis.title.y = element_blank(),
        panel.grid=element_blank(),
        legend.text=element_text(size=30),
        strip.text = element_text(size = 35))
# plot(p2)
ggsave(outfile_tier_pie, p2, width=10, height=8, units="in", dpi=1200)


# TODO: use functionality or package to have arrow poiting to insignificant percent in the pie chart, and mark 0% (or better, show some decimals when applicable)
 


## C1) DISTRIBUTION OF EVIDENCE TYPES (only for variants matched in CIVIC)

## Plot version showing absolute variant numbers
## (to avoid interpretation bias, since the same variant can be accounted for toward several evidence types)

# Sanity check that required input columns can be found
n_evidence_columns = c("n_predictive_vars", "n_diagnostic_vars", "n_prognostic_vars", "n_predisposing_vars",
                       "n_predictive_vars_no_tier3", "n_diagnostic_vars_no_tier3", "n_prognostic_vars_no_tier3", "n_predisposing_vars_no_tier3")

df_evidence_snv_sorted = prepare_data_for_plotting(data_snv, n_evidence_columns, input_colnames, "n_predictive_vars_no_tier3")
df_evidence_cnv_sorted = prepare_data_for_plotting(data_cnv, n_evidence_columns, input_colnames, "n_predictive_vars_no_tier3")

df_evidence_snv_sorted$evidence_type = gsub("_no_tier3", "", gsub("_vars", "", gsub("n_", "", df_evidence_snv_sorted$variable)))
df_evidence_snv_sorted$tier3_variants = ifelse(grepl("_no_tier3", df_evidence_snv_sorted$variable), "Excluding", "Including")
df_evidence_cnv_sorted$evidence_type = gsub("_no_tier3", "", gsub("_vars", "", gsub("n_", "", df_evidence_cnv_sorted$variable)))
df_evidence_cnv_sorted$tier3_variants = ifelse(grepl("_no_tier3", df_evidence_cnv_sorted$variable), "Excluding", "Including")


## Plot combined version of this plot, showing SNVs and CNVs as facets (same axis, use sample sorting from SNV dataset)

# Sanity check that required input columns can be found
plotting_columns_evidence = c("sample_name", "variable", "value", "label_ypos", "sample_order", "evidence_type", "tier3_variants")
df_evidence_fraction =  combine_snv_and_cnv_for_plotting(df_evidence_snv_sorted, df_evidence_cnv_sorted, plotting_columns_evidence)

# Choose color scheme for plot
tier3_colors = c("Including" = "grey", "Excluding" = "#0072B2")

# TODO: camel case for Evidence Types
# TODO: check Predisposing is never available and if so, entirely remove from plot
# TODO: use different layout -> I believe the X-axis just repeats the same set of Patients 4 times?


# NOTE: barplot shown below is stacked and not overlaying!

outfile_evidence = paste0(opt$outfile_tag, ".variants_evidence_type_distribution.pdf")
p3c = ggplot(df_evidence_fraction, aes(x=sample_order, y=value, fill=tier3_variants)) +
  geom_bar(stat='identity', position = "stack", width=.7) +
  ylab("Number of variants") + 
  xlab("Patient") +
  # ggtitle("") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  # theme(legend.justification = c(1, 1), legend.position = c(1,1)) +
  theme(legend.justification = c(1, 1)) +
  facet_grid(dataset ~ evidence_type) +
  # facet_grid(dataset ~ ., switch="both") +
  scale_fill_manual(name="Tier3 variants", values=tier3_colors)
plot(p3c)
# ggsave(outfile_evidence, p3c, width=10, height=8, units="in", dpi=1200)


## C2) TOTAL DISTRIBUTION OF EVIDENCE TYPES ACROSS PATIENTS




## D) DISEASE DISTRIBUTION OF CT CLASSIFICATION

# Sanity check that required input columns can be found
n_disease_columns = c("n_diseases", "n_diseases_ct", "n_diseases_gt", "n_diseases_nct")
df_disease_snv_sorted = prepare_data_for_plotting(data_snv, n_disease_columns, input_colnames, "n_diseases")
df_disease_cnv_sorted = prepare_data_for_plotting(data_cnv, n_disease_columns, input_colnames, "n_diseases")

# Sort ct classes by increasing order of priority
ct_order = c("n_diseases_nct", "n_diseases_gt", "n_diseases_ct")
df_disease_snv_sorted$ct_order = factor(df_disease_snv_sorted$variable, levels=ct_order)
df_disease_cnv_sorted$ct_order = factor(df_disease_cnv_sorted$variable, levels=ct_order)

# Remove number of diseases matched in CIViC (used only for sorting samples)
df_disease_snv_sorted = df_disease_snv_sorted[df_disease_snv_sorted$variable != "n_diseases",,drop=F]
df_disease_cnv_sorted = df_disease_cnv_sorted[df_disease_cnv_sorted$variable != "n_diseases",,drop=F]

# Choose color scheme for plot
disease_cols = c("n_diseases_ct" = "red", "n_diseases_gt" = "gold", "n_diseases_nct" = "skyblue")

# NOTE: barplots shown below are stacked and not overlaying!

outfile_disease_snv = paste0(opt$outfile_tag, ".snvs_ct_distribution.pdf")
p3a = ggplot(df_disease_snv_sorted, aes(x=sample_order, y=value, fill=ct_order)) +
  geom_bar(stat='identity', position = "stack", width=.7) +
  ylab("Number of disease names (SNVs)") + 
  xlab("Patient") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(legend.justification = c(1, 1),legend.position = c(1,1)) +
  # coord_cartesian(ylim=c(0, 30)) +
  # scale_y_continuous(breaks=seq(0, 70, by = 10)) +
  scale_fill_manual(name="Cancer specificity", breaks=c("n_diseases_ct", "n_diseases_gt", "n_diseases_nct"), labels=c("ct", "gt", "nct"), values=disease_cols)
# plot(p)
ggsave(outfile_disease_snv, p3a, width=10, height=8, units="in", dpi=1200)

outfile_disease_cnv = paste0(opt$outfile_tag, ".cnvs_ct_distribution.pdf")
p3b = ggplot(df_disease_cnv_sorted, aes(x=sample_order, y=value, fill=ct_order)) +
  geom_bar(stat='identity', position = "stack", width=.7) +
  ylab("Number of disease names (CNVs)") + 
  xlab("Patient") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(legend.justification = c(1, 1),legend.position = c(1,1)) +
  scale_fill_manual(name="Cancer specificity", breaks=c("n_diseases_ct", "n_diseases_gt", "n_diseases_nct"), labels=c("ct", "gt", "nct"), values=disease_cols)
# plot(p)
ggsave(outfile_disease_cnv, p3b, width=10, height=8, units="in", dpi=1200)
