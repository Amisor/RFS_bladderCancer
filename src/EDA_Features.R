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

# ------------------------ Histogram --------------------------
# Function to plot categorical variables
PlotCategorical <- function(t, cat_vars, colors){
  # One plot per variable
  par(mfrow = c(2, 3),
      mar   = c(4, 2, 2, 2),  # shrink margins around each subplot
      oma   = c(1, 1, 1, 1)   # no outer margins
  )
  
  for (i in 1:length(cat_vars)) {
    var <- cat_vars[i]
    bar_data <- table(t[[var]])
    barplot(bar_data,
            main = paste(var),
            col = colors[i],
            las = 2,  # rotate axis labels
            cex.names = 0.8, # smaller category labels
            cex.main = 0.8, # title smaller
            cex.axis = 0.8)  
  }
  par(mfrow = c(1, 1))  # Reset layout
}

PlotCategorical <- function(t, cat_vars, colors) {
  # One plot per variable
  par(mfrow = c(2, 3),
      mar   = c(4, 2, 2, 2),  # Shrink margins around each subplot
      oma   = c(1, 1, 1, 1)   # No outer margins
  )
  
  for (i in 1:length(cat_vars)) {
    var <- cat_vars[i]
    
    # Check if the variable is "Age"
    if (var == "Age") {
      hist(t[[var]], breaks = 4, 
           main = paste(var),
           col = colors[i],
           xlab = var,
           cex.main = 0.8,  # Smaller title
           cex.axis = 0.8)  # Smaller axis labels
    } else {
      bar_data <- table(t[[var]])
      barplot(bar_data,
              main = paste(var),
              col = colors[i],
              las = 2,  # Rotate axis labels
              cex.names = 0.8,  # Smaller category labels
              cex.main = 0.8,  # Title smaller
              cex.axis = 0.8)  
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
PlotCategorical(t, cat_vars, colors1)

# Plot validation
PlotCategorical(v, cat_vars, colors2)

unique(v$Tumor.size)
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
      x <- df_subset[[i]]
      y <- df_subset[[j]]
      
      # Remove rows with NA in either variable
      complete_cases <- complete.cases(x, y)
      x <- x[complete_cases]
      y <- y[complete_cases]
      
      print("Paso esta parte")
      cat("i es", i, "j es ", j)
      
      if (is.numeric(x) && is.numeric(y)) {
        corr_matrix[i, j] <- cor(x, y, method = "pearson")
      } else if ((is.factor(x) || is.character(x)) && (is.factor(y) || is.character(y))) {
        # Cramér's V
        tab <- table(x, y)
        corr_matrix[i, j] <- suppressWarnings(assocstats(tab)$cramer)
        print("Paso Cramer")
      } else if (is.numeric(x) && length(unique(y)) == 2) {
        # Point-biserial (numeric vs binary categorical)
        corr_matrix[i, j] <- cor(x, as.numeric(as.factor(y)), method = "pearson")
        print("Paso Point-biserial")
      } else if (is.numeric(y) && length(unique(x)) == 2) {
        corr_matrix[i, j] <- cor(as.numeric(as.factor(x)), y, method = "pearson")
        print("Paso pearson")
      } else {
        corr_matrix[i, j] <- NA  # not supported
        print("Paso nans")
      }
    }
  }
  return(corr_matrix)
}

# Version 2 of correlation of variables
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
CorrHeatmap(t_corr)
names(t)

# For computing the correlation, the NA values are omitted automatically right? 

# ----------------------- Time dependent variables: PFS, RFS, FU 

par(mfrow = c(1, 3))  # Show 3 plots side by side
hist(t$RFS_time, main = "RFS_time", xlab = "Time (days)", col = "skyblue", breaks = 30)
hist(t$PFS_time., main = "PFS_time", xlab = "Time (days)", col = "salmon", breaks = 30)
hist(t$FUtime_days., main = "FUtime_days", xlab = "Time (days)", col = "lightgreen", breaks = 30)
par(mfrow = c(1, 1))  # Reset layout

#plot(density(t$FUtime_days., na.rm = TRUE), col = "darkgreen", main = "Density of Follow-up, RFS, PFS", xlab = "Days")
#lines(density(t$RFS_time, na.rm = TRUE), col = "blue")
#lines(density(t$PFS_time., na.rm = TRUE), col = "red")
#legend("topright", legend = c("FUtime", "RFS", "PFS"), col = c("darkgreen", "blue", "red"), lty = 1)

par(mfrow = c(1, 3),
    mar   = c(4, 2, 2, 2),  # Shrink margins around each subplot
    oma   = c(1, 1, 1, 1))  # Show 3 plots side by side
plot(t$RFS_time, t$PFS_time., xlab = "RFS_time", ylab = "PFS_time", 
     main = "Recurrence vs Progression", col = "purple", pch = 19)
plot(t$RFS_time, t$FUtime_days., xlab = "RFS_time", ylab = "FUtime_days", 
     main = "Recurrence vs Follow-Up", col = "darkorange", pch = 19)
plot(t$PFS_time., t$FUtime_days., xlab = "PFS_time", ylab = "FUtime_days", 
     main = "Progression vs Follow-Up", col = "darkred", pch = 19)
par(mfrow = c(1, 1))  # Reset layout


# Recurrence before progression
table(t$RFS_time < t$PFS_time., useNA = "ifany")

# FU time longer than RFS or PFS
table(t$FUtime_days. < t$RFS_time, useNA = "ifany")
table(t$FUtime_days. < t$PFS_time., useNA = "ifany")

# Is PFS_time greater than RFS_time ? 

# What periods do we have

t$recur_before_prog <- with(t, ifelse(!is.na(RFS_time) & !is.na(PFS_time.), RFS_time < PFS_time., NA))
t$fu_before_events <- with(t, ifelse(!is.na(FUtime_days.), FUtime_days. < pmin(RFS_time, PFS_time., na.rm = TRUE), NA))

table(t$recur_before_prog, useNA = "ifany")
table(t$fu_before_events, useNA = "ifany")



# ---------------------------------Remove NA  -------------------------------

# Drop Na in Recurrence and RFS_time
before_t <- nrow(t)
t_new <- t[!is.na(t$Recurrence) | !is.na(t$RFS_time), ]
after_t <- nrow(t_new)

cat("Before removing rows with NA in Recurrence or RFS_time: ", 
    before_t,
    " \n After removing the rows with NA in Recurrence or RFS_time: ",
    after_t)

# Check that there are no NA in Recurrence and RFS_time
sum(is.na(t_new$Recurrence)) == 0
sum(is.na(t_new$RFS_time)) == 0

# Convert Recurrence to numerical variable
t_new$Recurrence <- as.numeric(as.character(t_new$Recurrence))

# Check that there are no error while converting the variable
sum(is.na(t_new$Recurrence)) == 0

str(t_new$exprs)
dim(t_new$exprs)

# ---------------------------Expression Data --------------------------
# Gene with no variance 
sd_vals <- apply(t_new$exprs, 2, sd)
which(sd_vals == 0)  # Gives column indices with zero variance

zero_var_genes <- names(sd_vals)[which(sd_vals == 0)]
zero_var_genes

View(t_new$exprs[, "LncRNA2747_ENSG00000256494"])

# Remove zero-variance columns 
t_new <- t_new$exprs[, sd_vals != 0]

# Detect most variable genes 
gene_vars <- apply(t$exprs, 2, var)
sorted_gene_vars <- sort(gene_vars, decreasing = TRUE)
head(sorted_gene_vars, 10)  # top 10 most variable genes
names(head(sorted_gene_vars, 10))

top_n <- 500
top_genes <- names(sorted_gene_vars)[1:top_n]


# Subset expression matrix to top genes
exprs_top <- t_new$exprs[, top_genes]

boxplot(exprs_top[, 1:50], 
        las = 2, main = "Top 50 Most Variable Genes", outline = FALSE)

# Scaling 
exprs_top_z <- scale(exprs_top)

names(t)

# Now safely z-score
z_exprs <- scale(t_exprs_filtered)

summary(z_exprs[, 1])  # mean ~0 and SD ~1 for a gene
unique(apply(z_exprs, 2, sd)) # SDs should be ~1, Por qué tengo un cero?

#t$exprs_z <- z_exprs

# -------------------------------PCA ----------------------------------
str(t_new$exprs)
dim(t_new$exprs)

# PCA over the expression data 
#pca <- prcomp(t_new$exprs, scale. = TRUE)

# Proportion of variance explained by each PC
#var_explained <- pca$sdev^2 / sum(pca$sdev^2)

# Cumulative variance
#cum_var_explained <- cumsum(var_explained)

# At leas 80% of the variance
#threshold <- 0.8
#num_pc <- which(cum_var_explained >= threshold)[1]

#cat("Number of PCs to reach 80% variance:", num_pc, "\n")

# Elbow plot
#plot(var_explained, type = "b", pch = 19,  xlab = "Principal Component", ylab = "Proportion of Variance Explained", main = "Elbow Plot (Scree Plot) of PCA", col = "darkgreen")

# Plot PCs cumulative variance
#plot(cum_var_explained, type = "o", pch = 1, xlab = "Number of Principal Components", ylab = "Cumulative Variance Explained", main = "Cumulative Variance Explained by PCA")
#abline(h = threshold, col = "red", lty = 2)
#abline(v = num_pc, col = "blue", lty = 2)
#text(x = num_pc, y = threshold + 0.02, 
#     labels = paste0("PC-", num_pc), col = "blue", cex = 0.9)

# Assign colors to each level
risk_colors <- as.numeric(t_new$EAU.risk)

# Plot PCA
plot(pca$x[, 1], pca$x[, 2],
     col = risk_colors,
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "PCA: Samples colored by EAU Risk")

# HEre 

#pca <- prcomp(t_new$exprs, scale. = TRUE)

#pc_scores <- pca$x  # Save scores

# Combine PC scores with metadata
#t_new <- cbind(t_new, pc_scores)

# --------------------------------Model--------------------------------
# Select columns for modelling 
cols_no_mod <- c("UROMOL.ID",
                 "Tumor.stage", # It's not in the validation
                 "Tumor.grade", # It's not in the validation
                 "Smoking" ,    # Not in the validation   
                 "Tumor.size", # Constant variable: Ta     
                 "Incident.tumor", # Constant variable: low grade,
                 "EAU.risk" # not in the validation, also compare to it
                 )

cols_keep <- setdiff(names(t_new), cols_no_mod)
cols_keep

# Filter columns
t_m <- t_new[,cols_keep]

# ------------------------- LASSO with COX ------------------------

# Create a design matrix that includes all predictors.
# The model.matrix function will automatically create dummy variables for factors.
# Remove the intercept (the first column) by selecting -1.
predictor_formula <- ~ Progression + 
  #PFS_time. + 
  #FUtime_days. + 
  Incident.tumor +
  Tumor.size + 
  Age +  
  Sex + 
  Concomitant.CIS + 
  BCG + 
  UROMOL2021.classification +
  exprs 

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

# Fit the final LASSO Cox model using the selected lambda
final_model <- glmnet(x, y, family = "cox", alpha = 1, lambda = best_lambda)

# View a summary of the final model
print(final_model)

coef_final <- coef(final_model)
selected_vars <- rownames(coef_final)[which(coef_final != 0)]
selected_vars

# ---------------------- Evaluation metrics: C- index ----------------------

# Obtain linear predictors from the final model (using type = "link")
lp <- predict(final_model, newx = x, type = "link")

# Compute the C-index using rcorr.cens:
c_index <- rcorr.cens(lp, Surv(t_model$RFS_time, t_model$Recurrence))
print(c_index)

c_index_neg <- rcorr.cens(-lp, Surv(t_model$RFS_time, t_model$Recurrence))
c_index_neg

# ------------------------------------ ROC -----------------------------


# Define time points of interest (for example, in days: 1-year, 3-year, 5-year)
times <- c(30, 45, 60) 

# Calculate time-dependent ROC curves using your linear predictors
roc_obj <- timeROC(T = t_model$RFS_time,
                   delta = t_model$Recurrence,
                   marker = as.vector(lp),
                   cause = 1,         # specify the event of interest (usually coded as 1)
                   weighting = "marginal",
                   times = times,
                   ROC = TRUE,
                   iid = TRUE)

# Print the ROC results (this includes the AUC values at the specified times)
print(roc_obj$AUC)

# Plot ROC curve for 3 months
plot(roc_obj, time = times[1], col = "red", title = "Time-dependent ROC Curves")
# Add ROC curve for 6 months
plot(roc_obj, time = times[2], col = "blue", add = TRUE)
# Add ROC curve for 9 months
plot(roc_obj, time = times[3], col = "green", add = TRUE)

legend("bottom",
       legend = c("10 days", "20 days", "30 days"),
       col = c("red", "blue", "green"),
       lwd = 2)



# -------------------------Kaplan Meier 

# Obtain risk scores from your final model
# (Depending on the direction, you might want to negate the risk score if higher scores imply lower risk)
risk_score <- as.vector(predict(final_model, newx = x, type = "link"))
# Optionally, if needed (see previous discussion about inversion):
risk_score <- -risk_score

# Create risk groups (using the median as a cutoff)
t_model$risk_group <- ifelse(risk_score > median(risk_score), "High Risk", "Low Risk")

# Fit the Kaplan-Meier survival curves based on risk group
km_fit <- survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = t_model)

# Plot the Kaplan-Meier curves
plot(km_fit, col = c("red", "blue"), lwd = 2,
     xlab = "Time", ylab = "Recurrence-Free Survival Probability",
     main = "Kaplan-Meier Curves by Risk Group")
legend("bottomleft", legend = c("High Risk", "Low Risk"), col = c("red", "blue"), lwd = 2)


# -------------------------------- Extra -----------------------





# Other cox 


# Cox model
cox_model <- coxph(Surv(RFS_time, Recurrence) ~ Age + Sex + Smoking + Concomitant.CIS +
                     Tumor.size + Incident.tumor + EAU.risk + BCG +
                     UROMOL2021.classification + exprs, data = t_m)

# Check the model summary
summary(cox_model)

# Compute the variance inflation factors (VIFs)
vif_values <- vif(cox_model)
print(vif_values)


