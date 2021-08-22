library(dplyr)
library(imputeLCMD)

set.seed(1)


preprocess <- function(metabolites_path, annotations_path, annotations=FALSE) {
  metabolites <- data.frame(read.csv(metabolites_path, header=TRUE, row.names='Sample', stringsAsFactors=FALSE))
  log_metabolites <- log2(metabolites)
  log_metabolites[log_metabolites=='-Inf'] <- NaN
  imputed_metabolites <- impute.MinProb(log_metabolites)
  if (annotations==TRUE) {
    annotations <- data.frame(read.csv(annotations_path, header=TRUE, row.names="Metabolite", stringsAsFactors=FALSE))
    return(list(t(imputed_metabolites), annotations))
  }
  else {
    return(list(t(imputed_metabolites)))
  }
}


analyze <- function(design_matrix, data, selector="design_matrixdietLow", annotations=FALSE) {
  df <- as.data.frame(data[[1]])
  reslist <- list()
  for (col in colnames(df)) {
    subdf <- df[col]
    colnames(subdf) <- "intensity"
    mylm <- lm(formula = intensity ~ design_matrix + 0, data=subdf)
    summed <- summary(mylm)
    ind_results <- summed[["coefficients"]][selector, c(1, 4)]
    reslist[[col]] <- ind_results
  }
  newdf <- data.frame(reslist)
  colnames(newdf) <- colnames(df)
  res <- t(newdf)
  colnames(res) <- c("coefficient", "p_value")
  results <- data.frame(res)
  results$bh <- p.adjust(results$p, method = "BH")
  
  if (annotations==TRUE) {
    merged <- merge(results, data[[2]], by="row.names")
    row.names(merged) <- merged$Row.names
    merged <- subset(merged, select = -c(Row.names))
    return(merged)
  }
  
  else {
    return(results)
  }
}


#Import, log2 transform, and impute metabolomics data for both polar and volatile data.
polar_metabolites <- preprocess(metabolites_path="metabolomics/inputs/filtered_polar_metabolites.csv", 
                                annotations_path="metabolomics/inputs/filtered_polar_annotations.csv",
                                annotations=TRUE)

volatile_metabolites <- preprocess(metabolites_path="metabolomics/inputs/filtered_volatile_metabolites.csv", 
                                  annotations_path=NULL,
                                  annotations=FALSE)

#Create design matrix
diet <- c(rep("High", times=6),
          rep("Low", times=6))
litter <- rep(c(rep("A", times=2),
                rep("B", times=2),
                rep("C", times=2)), times=2)
Design <- model.matrix(~diet + litter)

#Compare high and low dietary fiber conditions
polar_df <- analyze(Design, polar_metabolites, annotations=TRUE)
volatile_df <- analyze(Design, volatile_metabolites, annotations=FALSE)

polar_metabolites_analysis <- polar_df[order(polar_df$p_value),]
volatile_metabolites_analysis <- volatile_df[order(volatile_df$p_value),]


#Write analysis outputs
write.table(data.frame("Metabolite"=rownames(polar_metabolites_analysis), polar_metabolites_analysis),
            file="metabolomics/outputs/polar_metabolite_analysis.txt",
            sep ="\t", 
            row.names=FALSE)

write.table(data.frame("Metabolite"=rownames(volatile_metabolites_analysis), volatile_metabolites_analysis),
            file="metabolomics/outputs/volatile_metabolites_analysis.txt",
            sep ="\t", 
            row.names=FALSE)



### PLOTTING OF DATA ###
library(EnhancedVolcano)

plot_data <- function(analysis_df, metabolite_labels, shapes, pcut=.05, pt_size=4) {
  color_labels <- rep('grey', nrow(analysis_df))
  color_labels[which(analysis_df$bh <= .05 & analysis_df$coefficient >= 1)] <- 'green1'
  color_labels[which(analysis_df$bh <= .05 & analysis_df$coefficient <= -1)] <- 'blue1'

  names(color_labels)[which(color_labels == 'green1')] <- 'LF'
  names(color_labels)[which(color_labels == 'blue1')] <- 'HF'
  names(color_labels)[which(color_labels == 'grey')] <- 'Not significant'
  
  min_coef <- min(floor(analysis_df$coefficient))
  max_coef <- max(ceiling(analysis_df$coefficient))
  max_coef_all <- max(c(abs(min_coef), abs(max_coef)))
  max_p <- max(ceiling(-log10(analysis_df$p_value)))

  volcano_plot <- EnhancedVolcano(analysis_df,
                                  lab = rownames(analysis_df),
                                  title = NULL,
                                  subtitle = NULL,
                                  x = 'coefficient',
                                  y = 'p_value',
                                  xlab = bquote(~Log[2]~ "fold-change"),
                                  ylab = bquote(~-Log[10]~italic(P)),
                                  pCutoff = pcut,
                                  FCcutoff = 1,
                                  #pointSize = 4,
                                  pointSize = analysis_df$sizes,
                                  labSize = 4,
                                  boxedLabels = TRUE,
                                  xlim = c(-5, 5),
                                  ylim = c(0, 6),
                                  drawConnectors = TRUE,
                                  widthConnectors = 0.5,
                                  colConnectors = 'grey30',
                                  arrowheads = FALSE,
                                  colCustom = color_labels,
                                  selectLab = metabolite_labels,
                                  shapeCustom = shapes,
                                  colAlpha = .65,
                                  caption=NULL,
                                  cutoffLineType="blank"
  )
  return(volcano_plot)
}



### Selection of metabolites that are mentioned in the manuscript
aas <- polar_df[polar_df$Class == "Amino acid", ] # Amino acids
protein_aas <- aas[!(row.names(aas) %in% c("aminomalonic acid", "dehydroalanine", "N-methylalanine")),] #Protein aas sonly

select_polar_metabolites <- c("D-galacturonic acid",
                              "D-mannose",
                              "fructose",
                              "mannitol phosphate",
                              "scyllo-inositol",
                              "lactate",
                              "N-acetylneuraminic acid")

select_volatile_metabolites <- c("acetate",
                                 "propionate",
                                 "butanoic acid")


# polar_metabolites_analysis$class_id <- values[match(polar_metabolites_analysis$Class, index)]
# 
# Placeholder for volatiles so they can be stacked with polar metabolites
volatile_metabolites_analysis$KEGGid <- ""
volatile_metabolites_analysis$Class <- "Volatile"

full_analysis <- rbind(polar_metabolites_analysis, volatile_metabolites_analysis)

# Label corrections for classes
full_analysis$Class[full_analysis$Class == "Sterol"] <- "Lipid"
full_analysis$Class[full_analysis$Class == "Misc"] <- "Other"
full_analysis$Class[full_analysis$Class == "TCA"] <- "Other"
full_analysis$Class[full_analysis$Class == "Inorganic ion"] <- "Other"
full_analysis$Class[full_analysis$Class == "Volatile"] <- "Other"
full_analysis$Class[full_analysis$Class == "Volatile"] <- "Other"

# Label corrections for SCFA's (n_SCFA due to class number ordering on legend)
full_analysis$Class[rownames(full_analysis) == "acetate"] <- "n_SCFA"
full_analysis$Class[rownames(full_analysis) == "propionate"] <- "n_SCFA"
full_analysis$Class[rownames(full_analysis) == "butanoic acid"] <- "n_SCFA"

# Assigning shapes to different classes
index <- c("Carbohydrate", "Amino acid", "Lipid", "n_SCFA", "Other")
values <- c(16, 15, 17, 18, 1)


# Match values to index
full_analysis$class_id <- values[match(full_analysis$Class, index)]

full_select_metabolites <- c(select_polar_metabolites, select_volatile_metabolites)

# Assigning size based on whether or not metabolite is in full_select_metabolites
sizes <- c(ifelse(rownames(full_analysis) %in% full_select_metabolites, 6, 3))
full_analysis$sizes <- sizes

# Ordering by Class to get the legend output in better order
full_analysis_ordered <- full_analysis[order(as.character(full_analysis$Class)),]

#Correcting n_SCFA (for legend ordering)
full_analysis_ordered$Class[full_analysis_ordered$Class == "n_SCFA"] <- "SCFA"

custom_shapes <- full_analysis_ordered$class_id
names(custom_shapes) <- full_analysis_ordered$Class

full_plot <- plot_data(full_analysis_ordered, full_select_metabolites, custom_shapes)
full_plot


full_out_path="metabolomics/outputs/figures/full_volcano.png"
#metab_plot_nol <- full_plot + theme(legend.position="none")
ggsave(full_out_path, width=6.5, scale=2, dpi=600)












