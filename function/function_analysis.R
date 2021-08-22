library(plyr)
library(edgeR)
library(limma)

preprocess <- function(fn_path, pattern) {
  files <- list.files(path=fn_path, pattern=pattern, full.names=TRUE, recursive=FALSE)
  freq_dfs <- lapply(files, function(x) {
    fn_df <- as.data.frame(read.table(x, header=FALSE, sep="\n", quote="", col.names="fn", stringsAsFactors=FALSE))
    new_fn <- plyr::count(fn_df, "fn")
    colnames(new_fn)[2] <- basename(x)
    row.names(new_fn) <- new_fn$fn
    return(new_fn)
  })
  combined <- as.data.frame(freq_dfs[[1]]$fn)
  colnames(combined) <- c("fn")
  for (df in freq_dfs) {
    combined <- merge(combined, df, by=c("fn"), all=TRUE)
  }
  row.names(combined) <- combined$fn
  final_fn <- subset(combined, select=-c(fn))
  colnames(final_fn) <- c("A1HF", "A1LF", "A2HF", "A2LF", "B1HF", "B1LF", "B2HF", "B2LF", "C1HF", "C1LF", "C2HF", "C2LF")
  final_fn[is.na(final_fn)] <- 0
  final_fn <- as.matrix(final_fn)
  return(final_fn)
}

# Differential abundance analysis
analyze <- function(data, design) {
  y <- DGEList(counts=data,)
  keep <- filterByExpr(y, design)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  results <- glmLRT(fit)
  output <- data.frame(topTags(results, n=Inf))
  return(output)
}

# Run a single analysis
run_analysis <- function(id_folder, pattern, results_file, design) {
  processed_data_file <- paste(results_file, "_data.txt", sep="")
  processed_data <- preprocess(fn_path=id_folder, pattern=pattern)
  write.table(data.frame("Function"=rownames(processed_data), processed_data),
              file=processed_data_file,
              sep ="\t", 
              row.names=FALSE)
  
  results <- analyze(processed_data, design=design)
  write.table(data.frame("Function"=rownames(results), results),
              file=results_file, 
              sep ="\t", 
              row.names=FALSE)
  return(list(processed_data, results))
}

# Run a series of analyses
multi_analysis <- function(id_folder, output_prefix, patterns, design) {
  analysis_results <- list()
  for (pattern in patterns) {
    search_pattern <- paste("*.", pattern, ".out", sep="")
    out_file <- paste(output_prefix, pattern, "_analysis.txt", sep="")
    analysis_results[[pattern]] <- run_analysis(id_folder=id_folder,
                                                pattern=search_pattern,
                                                results_file=out_file,
                                                design=design)
  }
  return(analysis_results)
}

# Create design matrix
Litter <- c(rep("A", times=4),
            rep("B", times=4),
            rep("C", times=4))
Diet <- rep(c("High", "Low"), times=6)
coldata <- cbind(Litter, Diet)
Design <- model.matrix(~Litter + Diet)


# Run all analyses
runner <- multi_analysis(id_folder='function/inputs/', 
                         output_prefix="function/outputs/metagenome_",
                         patterns="go", 
                         design=Design)
