library(ggplot2)
library(dplyr)
library(tidyr)
library(vcd)
library(reshape2)
library(survival)
library(car)
library(glmnet)
library(timeROC)
library(Hmisc)
library(colorspace)

# -------------------------- Read data -------------------------------------

# Training data
t <- as.data.frame(readRDS("./data/UROMOL_TaLG.teachingcohort.rds"))
# Validation data
v <- as.data.frame(readRDS("./data/knowles_matched_TaLG_final.rds"))

# ------------------------ Exploratory Data Analysis ----------------------------

# Compare if the two datasets have the same column names

if (setequal(names(t),names(v)) == FALSE ){
  print("The training and validation datasets 
        don't contain the same variables")
}else{
    print("The training and validation 
          datasets contain the same variables")
  }

t_no_v <- setdiff(names(t), names(v)) 
cat("Training dataset variables not available in the validation dataset:", t_no_v)
# "UROMOL.ID"  "Smoking"  "Tumor.size"  "Incident.tumor" low grade "EAU.risk"   

v_no_t <- setdiff(names(v), names(t))
cat("Validation dataset variables not available in the training dataset:", v_no_t)
# "knowles_ID"

#-----------Standarized values in columns between the two datasets ---------

# For the training dataset
unique(t$UROMOL2021.classification)

t <- t %>%
  mutate(UROMOL2021.classification = recode_factor(UROMOL2021.classification,
                                 "Class 1" = "Class_1",
                                 "Class 2b" = "Class_2b",
                                 "Class 2a" = "Class_2a",
                                 "Class 3" = "Class_3"))
unique(t$UROMOL2021.classification)

# For the validation dataset
unique(v$UROMOL2021.classification)
# Now they both have the same labels

# ---------------------------Modify columns types --------------------------
#View(t) # Identify the columns types
names(t) # Retrieve names to mdoofy format

cols_to_factor <- c("Progression", 
                    "Recurrence", 
                    "Sex",
                    "Smoking",
                    "Tumor.stage",
                    "Tumor.grade",
                    "Concomitant.CIS",
                    "Tumor.size",
                    "Incident.tumor",
                    "EAU.risk",
                    "BCG",
                    "UROMOL2021.classification" 
                    )

FactorConvert <- function(t, cols_to_factor){
  for (col in cols_to_factor) {
    if (col %in% names(t)){
      t[[col]] <- as.factor(t[[col]])
    }
  }
  return(t)
}

# Transform variables to categorical
t <- FactorConvert(t, cols_to_factor)
v <- FactorConvert(v, cols_to_factor)

# Constant values
unique(t$Tumor.stage) 
unique(t$Tumor.grade)

# ---------------------Summary Table ----------------------------
# Function to obtain variables types, number of nans, levels or range
SummaryTable <- function(df) {
  variable_names <- names(df)
  
  variable_types <- sapply(df, function(x) {
    if (is.numeric(x)) {
      "Numeric"
    } else if (is.factor(x) || is.character(x)) {
      "Categorical"
    } else if (is.logical(x)) {
      "Logical"
    } else {
      class(x)[1]
    }
  })
  
  na_counts <- sapply(df, function(x) sum(is.na(x)))
  
  range_or_levels <- sapply(seq_along(df), function(i) {
    x <- df[[i]]
    if (is.numeric(x)) {
      if (all(is.na(x))) {
        return(NA)
      } else {
        return(round(max(x, na.rm = TRUE) - min(x, na.rm = TRUE), 3))
      }
    } else {
      return(length(unique(na.omit(x))))
    }
  })
  
  min_max_or_levels <- sapply(seq_along(df), function(i) {
    x <- df[[i]]
    if (is.numeric(x)) {
      if (all(is.na(x))) {
        return("[NA, NA]")
      } else {
        return(paste0("[", round(min(x, na.rm = TRUE), 3), ", ", round(max(x, na.rm = TRUE), 3), "]"))
      }
    } else {
      lvls <- unique(na.omit(x))
      if (length(lvls) > 5) {
        return(paste0(paste0(head(lvls, 5), collapse = ", "), ", ..."))
      } else {
        return(paste(lvls, collapse = ", "))
      }
    }
  })
  
  summary_df <- data.frame(
    Variable = variable_names,
    Type = variable_types,
    NA_Count = na_counts,
    Range_or_Levels = range_or_levels,
    Min_Max_or_Levels = min_max_or_levels,
    stringsAsFactors = FALSE
  )
  
  return(summary_df)
}

summary_t <- SummaryTable(t)
summary_v <- SummaryTable(v)

write.csv(summary_t, file = "./data/summary_t.csv", row.names = FALSE)
write.csv(summary_v, file = "./data/summary_v.csv", row.names = FALSE)

# ------------------------ Histogram --------------------------
# Function to plot categorical variables
PlotCategorical <- function(t, cat_vars, colors) {
  # One plot per variable
  par(mfrow = c(2, 3),
      mar   = c(4, 2, 2, 2),  # Shrink margins around each subplot
      oma   = c(1, 1, 1, 1)   # No outer margins
  )
  
  for (i in 1:length(cat_vars)) {
    var <- cat_vars[i]
    
    if (var == "Age") {
      # Get histogram info without plotting
      hist_info <- hist(t[[var]], breaks = 4, plot = FALSE)
      percentages <- round(hist_info$counts / sum(hist_info$counts) * 100, 1)
      # Plot histogram with extended ylim to accommodate labels
      hist(t[[var]], breaks = 4, 
           main = paste(var),
           col = colors[i],
           xlab = var,
           cex.main = 0.8,  # Smaller title
           cex.axis = 0.6,  # Smaller axis labels
           ylim = c(0, max(hist_info$counts) * 1.2))
      # Add percentages at the top of each bin
      text(hist_info$mids, hist_info$counts, 
           labels = paste0(percentages, "%"), 
           pos = 3, offset = 0.5, cex = 0.5, col = "black")
    } else {
      bar_data <- table(t[[var]])
      # Set ylim to extend above the highest bar
      bar_pos <- barplot(bar_data,
                         main = paste(var),
                         col = colors[i],
                         las = 2,         # Rotate axis labels
                         cex.names = 0.8, # Smaller category labels
                         cex.main = 0.8,  # Title smaller
                         cex.axis = 0.6,
                         ylim = c(0, max(bar_data) * 1.2))
      percentages <- round((bar_data / sum(bar_data)) * 100, 1)
      # Add percentages on top of each bar
      text(bar_pos, bar_data, 
           labels = paste0(percentages, "%"), 
           pos = 3, cex = 0.5, col = "black")
    }
  }
  
  par(mfrow = c(1, 1))  # Reset layout
}

# Get categorical columns
omit_vars <- c("UROMOL.ID", "knowles_ID", "Progression", "Smoking", "Tumor.stage",
               "Tumor.grade", "Tumor.size", "Incident.tumor", "EAU.risk")

# Categorical variables (shared in training and validation sets)
cat_vars <- names(t)[sapply(t, function(x) is.factor(x) || is.character(x))]
cat_vars <- cat_vars[!cat_vars %in% omit_vars]
cat_vars <- append(cat_vars, "Age")
cat_vars

# Colors
colors1 <- rainbow(length(cat_vars))
colors2 <- darken(rainbow(length(cat_vars)), amount = 0.3)  # 30% darker

# Plot training 
png(filename = "./images/t_barplots.png", width = 1200, height = 1100, res = 300)
PlotCategorical(t, cat_vars, colors1)
dev.off()

# Plot validation
png(filename = "./images/v_barplots.png", width = 1200, height = 1000, res = 300)
PlotCategorical(v, cat_vars, colors2)
dev.off()

# Counts
#lapply(cat_vars, function(var) table(t[[var]], useNA = "ifany"))

# -------------------------------- Correlation ------------------
# Function to compute correlation between the different variables 
CorrVariables <- function(df, cols_cor) {
  df_subset <- df[cols_cor]
  n <- length(cols_cor)
  corr_matrix <- matrix(NA, nrow = n, ncol = n,
                        dimnames = list(cols_cor, cols_cor))
  
  for (i in 1:n) {
    for (j in 1:n) {
      
      # Set self-correlation to 1
      if (i == j) {
        corr_matrix[i, j] <- 1
        next
      }
      
      x <- df_subset[[i]]
      y <- df_subset[[j]]
      
      # Remove rows with NA in either variable
      complete_cases <- complete.cases(x, y)
      x <- x[complete_cases]
      y <- y[complete_cases]
      
      if (is.numeric(x) && is.numeric(y)) {
        # Pearson correlation for numeric vs numeric
        corr_matrix[i, j] <- cor(x, y, method = "pearson")
      } else if ((is.factor(x) || is.character(x)) &&
                 (is.factor(y) || is.character(y))) {
        # Cramér's V for categorical vs categorical
        tab <- table(x, y)
        corr_matrix[i, j] <- suppressWarnings(assocstats(tab)$cramer)
      } else if (is.numeric(x) && length(unique(y)) == 2) {
        # Point-biserial for numeric vs binary categorical
        corr_matrix[i, j] <- cor(x, as.numeric(as.factor(y)), method = "pearson")
      } else if (is.numeric(y) && length(unique(x)) == 2) {
        corr_matrix[i, j] <- cor(as.numeric(as.factor(x)), y, method = "pearson")
      } else if (is.numeric(x) && (is.factor(y) || is.character(y))) {
        # Correlation ratio (eta) for numeric vs categorical (with >2 levels)
        y <- as.factor(y)
        eta <- sqrt(sum(tapply(x, y, function(sub) {
          n_i <- length(sub)
          (mean(sub) - mean(x))^2 * n_i
        })) / sum((x - mean(x))^2))
        corr_matrix[i, j] <- eta
      } else if (is.numeric(y) && (is.factor(x) || is.character(x))) {
        # Correlation ratio (eta) for categorical vs numeric
        x <- as.factor(x)
        eta <- sqrt(sum(tapply(y, x, function(sub) {
          n_i <- length(sub)
          (mean(sub) - mean(y))^2 * n_i
        })) / sum((y - mean(y))^2))
        corr_matrix[i, j] <- eta
      } else {
        corr_matrix[i, j] <- NA  # not supported
      }
    }
  }
  return(corr_matrix)
}

# Function to plot the heatmap of the correlation
CorrHeatmap <- function(corr_matrix) {
  melted_corr <- melt(corr_matrix, na.rm = TRUE)
  ggplot(melted_corr, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 9, hjust = 1)) +
    coord_fixed()
}

# Columns I don't want to correlate 
remove_cols <- c("UROMOL.ID", 
                 "exprs", 
                 "Tumor.stage", 
                 "Tumor.grade",
                 "Smoking" ,    # Not in the validation   
                 "Tumor.size", # Constant variable: Ta     
                 "Incident.tumor", # Constant variable: low grade,
                 "EAU.risk" # not in the validation, also compare to it
                 )

# Remove the columns I don't want
cols_cor <- setdiff(names(t), remove_cols)
cols_cor

#cols_cor[9]
# Compute the correlation among the different variables
t_corr <- CorrVariables(t,cols_cor)

# Plot correlations
png(filename = "./images/t_correlation.png", width = 1700, height = 1700, res = 300)
CorrHeatmap(t_corr)
dev.off()

# For computing the correlation, the NA values are omitted automatically right? 

# ----------------------- Time dependent variables: PFS, RFS, FU 

PlotTimeVariables <- function(t){
  par(mfrow = c(1, 3))  # Show 3 plots side by side
  hist(t$RFS_time, main = "RFS_time", xlab = "Time (months)", col = "skyblue", breaks = 15)
  hist(t$PFS_time., main = "PFS_time", xlab = "Time (months)", col = "salmon", breaks = 15)
  hist(t$FUtime_days., main = "FUtime_days", xlab = "Time (days)", col = "lightgreen", breaks = 15)
  par(mfrow = c(1, 1))  # Reset layout
  
}

png(filename = "./images/t_SurvivalTimes.png", width = 1400, height = 1300, res = 300)
PlotTimeVariables(t)
dev.off()

hist(v$RFS_time, main = "RFS_time", xlab = "Time (days)", col = "skyblue", breaks = 15)
hist(v$FUtime_days., main = "Follow-Up Time", xlab = "Time (days)", col = "skyblue", breaks = 15)


PlotRecurrence <- function(t){
  par(mfrow = c(1, 2))
  
  # Histogram for RFS_time when Recurrence is 0
  hist(t$RFS_time[t$Recurrence == 0],
       main = "RFS_time (Recurrence = 0)",
       xlab = "RFS_time",
       col = "lightblue")
  
  # Histogram for RFS_time when Recurrence is 1
  hist(t$RFS_time[t$Recurrence == 1],
       main = "RFS_time (Recurrence = 1)",
       xlab = "RFS_time",
       col = "salmon")
  
  # Reset layout to default
  par(mfrow = c(1, 1))
}

PlotRecurrence(t)

#PlotRecurrence(v)
# ---------------------------------Remove NA  -------------------------------

## Drop Na in Recurrence and RFS_time
# Training dataset
t_new <- t[!is.na(t$Recurrence) | !is.na(t$RFS_time), ]
all(is.na(t$Recurrence) == is.na(t$RFS_time))
# Validation dataset
v_new <- v[!is.na(v$Recurrence) | !is.na(v$RFS_time), ]
all(is.na(v$Recurrence) == is.na(v$RFS_time))
sum(is.na(v_new$RFS_time))
v_rec_0 <- v_new[v_new$Recurrence == 0, ]
nrow(v_rec_0)
sum(is.na(v_rec_0$RFS_time))

cat("Before removing rows with NA in Recurrence or RFS_time in t: ", 
    nrow(t),
    " \n After removing the rows with NA in Recurrence or RFS_time in t: ",
    nrow(t_new))

cat("Before removing rows with NA in Recurrence or RFS_time in v: ", 
    nrow(v),
    " \n After removing the rows with NA in Recurrence or RFS_time in v: ",
    nrow(v_new))


# Check that there are no NA in Recurrence and RFS_time
## Training dataset
sum(is.na(t_new$Recurrence)) == 0
sum(is.na(t_new$RFS_time)) == 0
## Validation dataset
sum(is.na(v_new$Recurrence)) == 0
sum(is.na(v_new$RFS_time)) == 0
sum(is.na(v_new$FUtime_days.)) == 0
sum(is.na(v_new$RFS_time)) 

t_new_recurrence_0 <- t_new[t_new$Recurrence == 0, ]
t_new_recurrence_0$RFS_time_days <- round(t_new_recurrence_0$RFS_time*30, 0)
cor.test(t_new_recurrence_0$RFS_time_days, t_new_recurrence_0$FUtime_days.)

#png(filename = "./images/t_RFS_FU_timedays.png", width = 1400, height = 1700, res = 300)
plot(t_new_recurrence_0$RFS_time, t_new_recurrence_0$FUtime_days, 
     xlab = "RFS time (days)", 
     ylab = "FU time (days)", 
     main = "Recurrence Free Survival and Follow-Up time \n when Recurrence = 0", 
     col = "purple",
     cex.main = 0.95,
     pch = 19)
#dev.off()

#View(t_new_recurrence_0[, c("RFS_time_days","FUtime_days.")])

# Fill values of v_new 
# Run this if imputing the values of RFS_time
#v_new$RFS_time[v_new$Recurrence == 0] <- v_new$FUtime_days.[v_new$Recurrence == 0] / 30

# Convert Recurrence to numerical variable
## Training dataset
t_new$Recurrence <- as.numeric(as.character(t_new$Recurrence))
# Validation dataset
v_new$Recurrence <- as.numeric(as.character(v_new$Recurrence))

# Check that there are no errors while converting the variable
sum(is.na(t_new$Recurrence)) == 0
sum(is.na(v_new$Recurrence)) == 0
sum(is.na(v_new$RFS_time))

# ---------------------------Expression Data --------------------------
# Dimension of the expression data
## Training dataset
str(t_new$exprs)
dim(t_new$exprs)
## Validation dataset
str(v_new$exprs)
dim(v_new$exprs)

# Find common genes between training and validation datasets
common_genes <- intersect(colnames(t_new$exprs), 
                          colnames(v_new$exprs))

# Keep only common genes
t_new$exprs <- t_new$exprs[, common_genes, drop = FALSE]
v_new$exprs <- v_new$exprs[, common_genes, drop = FALSE]

# Check that dimensions are the same (number of genes only)
dim(t_new$exprs)[2] == dim(v_new$exprs)[2]
dim(t_new$exprs) # I have now 19087 genes

# Important genes for recurrence
genes <- c("TERT","EGFR","ATM", "NF1","NOTCH1","MYCL")

# Check if any of the genes are in the expression data (column names)
## Training dataset
genes %in% colnames(t_new$exprs)
## Validation dataset
genes %in% colnames(v_new$exprs)


# Function to detect the most variable genes
# genes with no variance and to subset the
# dataset only to the top variable genes
# Second try of variable genes
VariableGenes <- function(t_new, genes, alpha) {
  # Calculate variance for each gene
  t_var <- apply(t_new$exprs, 2, var)
  
  # Use the median variance as a baseline reference
  baseline_var <- median(t_var)
  
  n <- nrow(t_new$exprs)  # number of samples
  
  # For each gene, compute a p-value testing whether its variance is significantly greater than baseline_var.
  # Under the null, (n-1)*S^2 / baseline_var ~ Chi-square(df = n-1)
  pvals <- sapply(t_var, function(S2) {
    pchisq((n - 1) * S2 / baseline_var, df = n - 1, lower.tail = FALSE)
  })
  
  # Select genes with p-value less than alpha (i.e., significantly more variable than the baseline)
  variable_genes <- names(pvals)[pvals < alpha]
  
  # Always include the genes in your provided vector, even if they don't pass the test
  total_genes <- unique(c(variable_genes, genes))
  
  # Also, detect genes with zero variance (if any)
  t_zero_var_genes <- names(t_var)[t_var == 0]
  
  # Subset the expression matrix to only include the selected genes
  t_new_genes <- t_new$exprs[, total_genes]
  
  result <- list(t_zero_var_genes = t_zero_var_genes, 
                 selected_genes = total_genes, 
                 t_new_genes = t_new_genes)
  return(result)
}


# Call funciton to identify top genes, genes with no variance and the subset of the dataframe with the top genes
#list_var_genes <- VariableGenes(t_new, 1000, genes)
list_var_genes <- VariableGenes(t_new, genes, alpha = 0.05)
t_zero_var_genes <- list_var_genes[[1]]
total_genes <- list_var_genes[[2]]
t_new_genes <- list_var_genes[[3]]

# Now I have 8643 variable genes
length(total_genes)
# View genes
#View(t_new_genes)

# Plot
png(filename = "./images/t_TopVariableGenes.png", width = 1400, height = 1400, res = 300)
boxplot(t_new_genes[,1:20], 
        las = 2, main = "Top 20 Most Variable Genes", outline = FALSE, col = rainbow(20), cex=0.9)
dev.off()

genes_interes <- c("TAF4", "ZBTB6", "DSC2", "CDKN2B", "TNFRSF10D", "ARHGEF25",  "BEST1",
                   "MYO15A","BRSK1", "AATK", "PRR19",  "CHKB-DT", "AGGF1", "C10orf62", 
                   "LINC02610", "RSPO2") 
setdiff(genes_interes, colnames(t_new_genes))
# Plot genes of interest

png(filename = "./images/t_PredictorGenes.png", width = 1400, height = 1400, res = 300)
par(oma = c(1, 1, 1, 1)) 
boxplot(t_new_genes[,genes_interes], 
        las = 2, main = "Predictor Genes", outline = FALSE, col = rainbow(20), cex=0.7)
dev.off()
# Kepp only variable genes in t
t_new$exprs <- t_new_genes
# Save z-score
t_new$exprs_z <-  scale(t_new_genes)
dim(t_new$exprs)
dim(t_new$exprs_z)

# Keep the genes in v
v_new$exprs <- v_new$exprs[, total_genes, drop = FALSE]
# Save z-scores
v_new$exprs_z <- scale(v_new$exprs)
dim(v_new$exprs)
dim(v_new$exprs_z)
# --------------------------------Model--------------------------------
# Select columns for modelling 
cols_keep <- c("Recurrence",
               "RFS_time",
               "Sex",
               "Age",
               "Concomitant.CIS",
               "UROMOL2021.classification",
               "BCG", 
               "exprs_z"
               )

## Filter columns
# Training dataset
t_m <- t_new[,cols_keep]

# Validation dataset
v_m <- v_new[,cols_keep]

# unify factor variables that in v are only one level but in t are two
v_m$Sex <- factor(v_m$Sex, levels = levels(t_m$Sex))
v_m$BCG <- factor(v_m$BCG, levels = levels(t_m$BCG))
v_m$Concomitant.CIS <- factor(v_m$Concomitant.CIS, levels = levels(t_m$Concomitant.CIS))
v_m$UROMOL2021.classification <- factor(v_m$UROMOL2021.classification, levels = levels(t_m$UROMOL2021.classification))

# ------------------------- LASSO with COX ------------------------

# Create a design matrix that includes all predictors.
# The model.matrix function will automatically create dummy variables for factors.
# Remove the intercept (the first column) by selecting -1.
predictor_formula <- ~ Age +  
  Sex + 
  Concomitant.CIS + 
  BCG + 
  UROMOL2021.classification +
  exprs_z 

# Define all relevant columns: survival + predictors
cols_required <- c("RFS_time", "Recurrence", all.vars(predictor_formula))

# Drop rows with NA in any of them
t_model <- t_m[complete.cases(t_m[, cols_required]), ]

# Create the model matrix (without the intercept)
x <- model.matrix(predictor_formula, data = t_model)[, -1]

# Create the survival response
y <- with(t_model, Surv(RFS_time, Recurrence))

# Use cross-validation to select the optimal lambda for LASSO (alpha = 1)
cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1)

# Get the best lambda value
best_lambda <- cv_fit$lambda.min
best_lambda

# Fit the final LASSO Cox model using the selected lambda
final_model <- glmnet(x, y, family = "cox", alpha = 1, lambda = best_lambda)

# View a summary of the final model
print(final_model)

coef_final <- coef(final_model)
# Convert coefficients to a matrix and exponentiate them
hr <- exp(as.matrix(coef_final))
print(hr[hr!=1,])
      
selected_vars <- rownames(coef_final)[which(coef_final != 0)]
selected_vars

# Obtain linear predictors from the final model (using type = "link")
lp <- predict(final_model, newx = x, type = "link")

# -------------------------------Training: ROC -----------------------------

# Define time points of interest
times <- c(3, 6, 12, 24, 36, 60) 

# Calculate time-dependent ROC curves using your linear predictors
roc_obj_t <- timeROC(T = t_model$RFS_time,
                     delta = t_model$Recurrence,
                     marker = as.vector(lp),
                     cause = 1,         # specify the event of interest (usually coded as 1)
                     weighting = "marginal",
                     times = times,
                     ROC = TRUE,
                     iid = TRUE)

# Print the ROC results (this includes the AUC values at the specified times)
print(roc_obj_t$AUC)

png(filename = "./images/t_ROC.png", width = 1400, height = 1400, res = 300)
# Plot ROC curve for 3 months
plot(roc_obj_t, time = times[1], col = "red", title = "", add = FALSE)
# Add ROC curve for 6 months
plot(roc_obj_t, time = times[2], col = "blue", add = TRUE)
# Add ROC curve for 12 months
plot(roc_obj_t, time = times[3], col = "darkgreen", add = TRUE)
# Add ROC curve for 24 months
plot(roc_obj_t, time = times[4], col = "purple", add = TRUE)
# Add ROC curve for 36 months
plot(roc_obj_t, time = times[5], col = "skyblue", add = TRUE)
# Add ROC curve for 60 months
plot(roc_obj_t, time = times[6], col = "brown", add = TRUE)
title("Training Time-dependent ROC Curves")
legend("bottomright",
       legend = c("3 months", "6 months", "1 year", "2 years", "3 years", "5 years"),
       col = c("red", "blue", "darkgreen", "purple", "skyblue", "brown"),
       lwd = 2, cex = 0.8)
dev.off()
# -------------------------Training: Kaplan Meier -------------------------------

# Obtain risk scores from your final model
risk_score <- as.vector(predict(final_model, newx = x, type = "link"))

# Create risk groups
t_model$risk_group <- ifelse(risk_score > median(risk_score), "High Risk", "Low Risk")

# Fit the Kaplan-Meier survival curves based on risk group
km_fit <- survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = t_model)

# Plot the Kaplan-Meier curves
png(filename = "./images/t_KaplanMeier.png", width = 1400, height = 1400, res = 300)
plot(km_fit, col = c("red", "blue"), lwd = 2,
     xlab = "Time", ylab = "Recurrence-Free Survival Probability",
     main = "Training Kaplan-Meier Curves")
legend("bottomleft", legend = c("High Risk", "Low Risk"), col = c("red", "blue"), lwd = 2, cex = 0.8)
dev.off()

# --------------------------Validation -----------------------

# Define the columns required for the predictors and survival outcome
cols_required <- c("RFS_time", "Recurrence", all.vars(predictor_formula))

# Subset the validation data to keep only rows with complete cases
v_model <- v_m[complete.cases(v_m[, cols_required]), ]

# Create the model matrix for the validation set (drop the intercept column)
x_val <- model.matrix(predictor_formula, data = v_model)[, -1]

# Create the survival response for validation data
y_val <- with(v_model, Surv(RFS_time, Recurrence))

# Obtain linear predictors from the final model
lp_val <- predict(final_model, newx = x_val, type = "link")

# -------------------------------- Validation: ROC ------------------------
# Define time points of interest
times <- c(3, 6, 12, 24, 36, 60) 

# Calculate time-dependent ROC curves using your linear predictors
roc_obj_v <- timeROC(T = v_model$RFS_time,
                     delta = v_model$Recurrence,
                     marker = as.vector(lp_val),
                     cause = 1,       
                     weighting = "marginal",
                     times = times,
                     ROC = TRUE,
                     iid = TRUE)

# Print the ROC results (this includes the AUC values at the specified times)
print(roc_obj_v$AUC)

#png(filename = "./images/v_ROC_nofillofRFS_time.png", width = 1400, height = 1400, res = 300)
png(filename = "./images/v_ROC.png", width = 1400, height = 1400, res = 300)
# Plot ROC curve for 3 months
plot(roc_obj_v, time = times[1], col = "red", title = "", add = FALSE)
# Add ROC curve for 6 months
plot(roc_obj_v, time = times[2], col = "blue", add = TRUE)
# Add ROC curve for 12 months
plot(roc_obj_v, time = times[3], col = "darkgreen", add = TRUE)
# Add ROC curve for 24 months
plot(roc_obj_v, time = times[4], col = "purple", add = TRUE)
# Add ROC curve for 36 months
plot(roc_obj_v, time = times[5], col = "skyblue", add = TRUE)
# Add ROC curve for 60 months
plot(roc_obj_v, time = times[6], col = "brown", add = TRUE)
title("Validation Time-dependent ROC Curves")
legend("bottomright",
       legend = c("3 months", "6 months", "1 year", "2 years", "3 years", "5 years"),
       col = c("red", "blue", "darkgreen", "purple", "skyblue", "brown"),
       lwd = 2, cex = 0.8)
dev.off()


# ------------------------------ Validation: Kaplan Meier --------------------------
# Obtain risk scores from your final model
risk_score_v <- -as.vector(predict(final_model, newx = x_val, type = "link"))

# Create risk groups
v_model$risk_group <- ifelse(risk_score_v > median(risk_score_v), "High Risk", "Low Risk")

# Fit the Kaplan-Meier survival curves based on risk group
v_km_fit <- survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = v_model)

# Plot the Kaplan-Meier curves
png(filename = "./images/v_KaplanMeier.png", width = 1400, height = 1400, res = 300)
plot(v_km_fit, col = c("red", "blue"), lwd = 2,
     xlab = "Time", ylab = "Recurrence-Free Survival Probability",
     main = "Validation Kaplan-Meier Curves ")
legend("bottomleft", legend = c("High Risk", "Low Risk"), col = c("red", "blue"), lwd = 2, cex = 0.8)
dev.off()


