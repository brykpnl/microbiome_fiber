library(PECA)
library(imputeLCMD)
library(dplyr)
library(tibble)
set.seed(1)

preprocess <- function(ms_data, log2_transform=TRUE) {
  data <- read.table(ms_data, header=TRUE, sep=',', na.strings=c(""," ", "NA"), stringsAsFactors=FALSE)
  ids <- data[,1:2]
  if (log2_transform == TRUE) {
    intensities <- log2(sapply(data[, -(1:2)], function(x) as.numeric(x)))
  }
  else {
    intensities <- sapply(data[, -(1:2)], function(x) as.numeric(x))
  }
  data <- cbind(ids, intensities)
  split_peptide <- strsplit(as.character(data$Peptide),'.',fixed=TRUE)
  Peptide_mod <- as.character(lapply(split_peptide, function(x) paste(x[c(3:length(x)-1)], collapse = '')))
  data_mod <- add_column(data, Peptide_mod, .after = 2)
  return(data_mod)
}

drop_contaminants <- function(data, ids) {
  non_contaminants <- paste(ids,collapse="|")
  contaminants <- dplyr::filter(data, !grepl(non_contaminants, Protein))
  drop_peptides <- contaminants['Peptide_mod']
  dropped_data <- data[!(data$Peptide_mod %in% drop_peptides$Peptide_mod), ]
  return(dropped_data)
}

remove_missing <- function(samples, min_obs=3) {
  passing_peptides <- row.names(samples[rowSums(!is.na(samples)) >= min_obs, ])
  return(passing_peptides)
}

peptide_select <- function(data) {
  unique_data <- distinct(data)
  data_u_rm <- unique_data[remove_missing(unique_data[c(-1)]),]
  return(data_u_rm)
}

clean <- function(ms_data,
                  extras,
                  grep_ids=c('NODE', 'MOUSE'),
                  log2_transform=TRUE) {
  data <- preprocess(ms_data, log2_transform)
  data <- drop_contaminants(data, grep_ids)
  return(data)
}

prepare <- function(data,
                    condition_A,
                    condition_B) {
  A <- data[,c("Peptide_mod", condition_A)]
  B <- data[,c("Peptide_mod", condition_B)]
  samples <- list(A, B)
  select_peptides <- lapply(samples, function(x) {peptide_select(x)})
  all_peptides = list()
  for (sample in select_peptides) {
    all_peptides <- rbind(all_peptides, sample$Peptide_mod)
  }
  all_labels <- data[data$Peptide_mod %in% all_peptides, c("Protein", "Peptide_mod")]
  all_shared <- data[data$Peptide_mod %in% all_peptides, c("Peptide_mod", condition_A, condition_B)]
  all_distinct <- distinct(all_shared)
  A_imputed <- impute.MinProb(subset(all_distinct, select=c(condition_A)))
  B_imputed <- impute.MinProb(subset(all_distinct, select=c(condition_B)))
  unique_data <- cbind(all_distinct$Peptide_mod, A_imputed, B_imputed)
  colnames(unique_data)[1] <- "Peptide_mod"
  corrected_columns <- c("Protein", "Peptide_mod", condition_A, condition_B)
  full_data <- merge(all_labels, unique_data, by=c("Peptide_mod"))[, corrected_columns]
  single_ids <- full_data[!(duplicated(full_data$Peptide_mod) | duplicated(full_data$Peptide_mod, fromLast = TRUE)), ]$Peptide_mod
  deduplicated_data <- full_data[full_data$Peptide_mod %in% single_ids, ]
  return(list(full_data, deduplicated_data))
}

compare <- function(data, condition_A, condition_B, path) {
  colnames(data)[1:2] <- c("Accession", "Sequence")
  write.table(data, path, row.names=FALSE, quote=FALSE, sep="\t")
  peptide_fold_changes <- data.frame(apply(data, 1, function(x) (mean(as.numeric(x[condition_A]))) - (mean(as.numeric(x[condition_B])))))
  colnames(peptide_fold_changes) <- "FC"
  peptide_fold_changes$Accession <- data$Accession
  protein_fold_change <- aggregate(peptide_fold_changes$FC, by=list(Acc=peptide_fold_changes$Accession), FUN=mean)
  colnames(protein_fold_change) <- c("Accession", "FC")
  rownames(protein_fold_change) <- protein_fold_change$Accession
  AA_results <- PECA_tsv(file=path, samplenames1=condition_A, samplenames2=condition_B, normalize=FALSE,
                         test="modt", type="median", paired=TRUE, progress=FALSE)
  AA_results$pBH <- p.adjust(AA_results$p, method='BH')
  AA_results$FC <- protein_fold_change$FC
  return(AA_results)
}

significant_proteins <- function(AA_results, pBH=.05, FC=2) {
  return(row.names(AA_results[(AA_results$pBH < pBH) & (AA_results$FC > FC | AA_results$FC < -FC),]))
}



#Path to peptide abundances and column names for each sample
ms_data <- "./abpp/inputs/all_peptide_abundances.csv"

hi_probe <- c("A_Hi_GH", "B_Hi_GH", "C_Hi_GH")
lo_probe <- c("A_Lo_GH", "B_Lo_GH", "C_Lo_GH")
hi_noprobe <- c("A_Hi_NP", "B_Hi_NP", "C_Hi_NP")
lo_noprobe <- c("A_Lo_NP", "B_Lo_NP", "C_Lo_NP")

#Clean data
data_clean <- clean(ms_data)
dir.create("./abpp/intermediate_outputs", showWarnings=FALSE)
write.table(data_clean, "./abpp/intermediate_outputs/data_clean.txt", row.names=FALSE, quote=FALSE, sep="\t")

#Prepare clean data for comparisons between conditions
hilo_data <- prepare(data_clean, hi_probe, lo_probe)
hi_data <- prepare(data_clean, hi_probe, hi_noprobe)
lo_data <- prepare(data_clean, lo_probe, lo_noprobe)


#Compare proteins between probe and no probe conditions in the high and low fiber conditions using all peptides or unique peptides only
hi_probe_comp <- compare(hi_data[[1]], 
                         hi_probe, 
                         hi_noprobe, 
                         path="./abpp/intermediate_outputs/hi_probe_comp.txt")
hi_probe_comp_unique <- compare(hi_data[[2]], 
                         hi_probe, 
                         hi_noprobe, 
                         path="./abpp/intermediate_outputs/hi_probe_comp_unique.txt")
lo_probe_comp <- compare(lo_data[[1]], 
                         lo_probe, 
                         lo_noprobe, 
                         path="./abpp/intermediate_outputs/lo_probe_comp.txt")
lo_probe_comp_unique <- compare(lo_data[[2]], 
                         lo_probe, 
                         lo_noprobe, 
                         path="./abpp/intermediate_outputs/lo_probe_comp_unique.txt")

#Obtain significantly different proteins for each comparison
sig_hi <- significant_proteins(hi_probe_comp)
sig_lo <- significant_proteins(lo_probe_comp)
sig_hi_unique <- significant_proteins(hi_probe_comp_unique)
sig_lo_unique <- significant_proteins(lo_probe_comp_unique)

#Combine all unique target identifiers for all probe vs. no-probe results
all_targets <- unique(c(sig_hi, sig_hi_unique, sig_lo, sig_lo_unique))

#Select unique peptides from potential probe targets in the prepared hi vs. lo data
hilo_targets <- hilo_data[[2]][(hilo_data[[2]]$Protein) %in% all_targets,]

#Compare probe-targeted proteins in hi and lo probe conditions
hilo_comp_unique <- compare(hilo_targets, 
                            lo_probe, 
                            hi_probe, 
                            path="./abpp/intermediate_outputs/hilo_comp_unique.csv")


write.table(hilo_comp_unique, "./abpp/outputs/hilo_comp_unique_analysis.txt", row.names=TRUE, quote=FALSE, sep="\t")

# Group protein-level intensities

reg_data <- data.frame(hilo_data[[2]][1:2], 2**hilo_data[[2]][, 3:length(hilo_data[[2]])])

protein_only <- subset(reg_data, select = -c(Peptide_mod))


protein_sums <- protein_only %>% group_by(Protein) %>% summarise(across(everything(), sum))
protein_means <- protein_only %>% group_by(Protein) %>% summarise(across(everything(), mean))
protein_sums_log2 <- data.frame(protein_sums[1], log2(protein_sums[2:length(protein_sums)]))
protein_means_log2 <- data.frame(protein_sums[1], log2(protein_means[2:length(protein_means)]))

# Select unique peptides, then only proteins that are analyzed and then apply to the full data
data_clean_notransform <- data.frame(data_clean[1:3], 2**data_clean[4:length(data_clean)])

full_distinct <- data_clean_notransform %>% distinct(Peptide_mod, .keep_all = T)
full_protein_only <- subset(full_distinct, select = -c(Peptide, Peptide_mod))

select_proteins <- hilo_data[[2]]["Protein"]
select_full_data <- full_protein_only[which(full_protein_only$Protein %in% select_proteins$Protein), colnames(protein_only)]

full_protein_sums <- aggregate(. ~ Protein, select_full_data, FUN=function(x) sum(x, na.rm=TRUE),  na.action=na.pass)
full_protein_means<- aggregate(. ~ Protein, select_full_data, FUN=function(x) mean(x, na.rm=TRUE),  na.action=na.pass)

full_protein_sums_log2 <- data.frame(full_protein_sums[1], log2(full_protein_sums[2:length(full_protein_sums)]))
full_protein_means_log2 <- data.frame(full_protein_means[1], log2(full_protein_means[2:length(full_protein_means)]))

full_protein_sums_log2[full_protein_sums_log2==-Inf] = NA


# Write sums and means out
write.table(full_protein_sums_log2, "./abpp/outputs/full_protein_sums_log2.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(full_protein_means_log2, "./abpp/outputs/full_protein_means_log2.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(protein_sums_log2, "./abpp/outputs/protein_sums_log2.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(protein_means_log2, "./abpp/outputs/protein_means_log2.txt", row.names=FALSE, quote=FALSE, sep="\t")