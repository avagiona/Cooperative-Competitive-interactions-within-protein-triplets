#REQUIREMENTS
library("dplyr")
library("tidyr")
library("tibble")
library("writexl")
library("ggplot2")
library('future') #parallezation 
library("furrr") #parallezation
#set pathway
setwd("~/R-4.2.1/predict_complexes_ML")

# Setup parallel backend
plan(multisession, workers = parallel::detectCores() - 1)
# PARALLELIZED SCRIPT: IDENTIFY OPEN TRIANGLES + ANNOTATE INTERACTION REGIONS

# Step 1: Load structural interactions data
struct_intr <- read.delim("interactions.dat", header = TRUE, sep = "")

# Step 2: Clean and filter interactions
struct_intr_clean <- struct_intr %>%
  filter(TYPE == "Structure", PROT1 != PROT2) %>%
  distinct(PDB_ID, PROT1, PROT2, CHAIN1, CHAIN2, SEQ_BEGIN1, SEQ_END1, SEQ_BEGIN2, SEQ_END2)

# Step 3: List unique proteins per PDB
proteins_per_pdb <- struct_intr_clean %>%
  pivot_longer(cols = c(PROT1, PROT2), names_to = "which_prot", values_to = "protein") %>%
  group_by(PDB_ID) %>%
  summarise(proteins = paste(unique(protein), collapse = " "), .groups = "drop")

# Step 4: Generate triplets per PDB (in parallel)
triplets_input <- proteins_per_pdb %>%
  mutate(
    protein_list = strsplit(proteins, " "),
    num_proteins = lengths(protein_list)
  ) %>%
  filter(num_proteins >= 3)

triplet_lists <- future_map(triplets_input$protein_list, ~combn(.x, 3, simplify = FALSE))

triplets_per_pdb <- tibble(
  PDB_ID = rep(triplets_input$PDB_ID, lengths(triplet_lists)),
  triplets = unlist(triplet_lists, recursive = FALSE)
)

# Step 5: Prepare validated interaction pairs with region info
validated_pairs <- struct_intr_clean %>%
  mutate(
    pair = paste(pmin(PROT1, PROT2), pmax(PROT1, PROT2), sep = "_"),
    region_info = paste0(
      PROT1, ":", CHAIN1, "(", SEQ_BEGIN1, "-", SEQ_END1, ") | ",
      PROT2, ":", CHAIN2, "(", SEQ_BEGIN2, "-", SEQ_END2, ")"
    )
  ) %>%
  select(PDB_ID, pair, region_info)

validated_lookup <- split(validated_pairs, validated_pairs$PDB_ID)

# Step 6: Identify open triangles (2 validated interactions out of 3)
open_triangles <- future_map2_dfr(
  triplets_per_pdb$PDB_ID,
  triplets_per_pdb$triplets,
  function(pdb, triplet) {
    triplet <- sort(triplet)
    A0 <- triplet[1]; B0 <- triplet[2]; C0 <- triplet[3]
    pairs0 <- list(
      AB = paste(sort(c(A0, B0)), collapse = "_"),
      AC = paste(sort(c(A0, C0)), collapse = "_"),
      BC = paste(sort(c(B0, C0)), collapse = "_")
    )
    
    if (!(pdb %in% names(validated_lookup))) return(NULL)
    
    validated <- validated_lookup[[pdb]]$pair
    regions <- setNames(validated_lookup[[pdb]]$region_info, validated_lookup[[pdb]]$pair)
    
    is_valid <- sapply(pairs0, function(p) p %in% validated)
    edge_count <- sum(is_valid)
    
    if (edge_count == 2) {
      # Determine common interactor
      common <- NULL
      if (is_valid["AB"] & is_valid["AC"]) common <- A0
      else if (is_valid["AB"] & is_valid["BC"]) common <- B0
      else if (is_valid["AC"] & is_valid["BC"]) common <- C0
      
      others <- setdiff(c(A0, B0, C0), common)
      A <- common; B <- others[1]; C <- others[2]
      
      # Rebuild pair names
      AB_name <- paste(sort(c(A, B)), collapse = "_")
      AC_name <- paste(sort(c(A, C)), collapse = "_")
      BC_name <- paste(sort(c(B, C)), collapse = "_")
      
      # Recalculate validations
      AB_valid <- AB_name %in% validated
      AC_valid <- AC_name %in% validated
      BC_valid <- BC_name %in% validated
      
      AB_region <- if (AB_valid) regions[[AB_name]] else NA
      AC_region <- if (AC_valid) regions[[AC_name]] else NA
      
      return(tibble(
        PDB_ID = pdb,
        A = A, B = B, C = C,
        AB_valid = AB_valid, AB_region = AB_region,
        AC_valid = AC_valid, AC_region = AC_region,
        BC_valid = BC_valid
      ))
    } else {
      return(NULL)
    }
  }
)

# Save open triangles to file
#write.csv(open_triangles, "open_triangles.csv", row.names = FALSE)

# =======================
# PLOT 1: Proteins per PDB
# =======================
protein_count_plot <- triplets_input %>%
  ggplot(aes(x = num_proteins)) +
  geom_histogram(bins = 40, fill = "#C90DFF", color = "black", alpha = 0.7) +
  theme_bw() +
  labs(x = "Number of proteins/PDB id",
       y = "Frequency")

print(protein_count_plot)

# ================================
# PLOT 2: Open triangles per PDB
# ================================
open_triangle_counts <- open_triangles %>%
  count(PDB_ID, name = "num_open_triangles")

open_triangle_plot <- ggplot(open_triangle_counts, aes(x = num_open_triangles)) +
  geom_histogram(binwidth = 1, fill = "#C90DFF", color = "black", alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_bw() +
  labs(
    x = "Number of cooperative triplet complexes / PDB id",
    y = "Frequency"
  )

print(open_triangle_plot)

library(ggpubr)
pdb_stats<-ggarrange(protein_count_plot,open_triangle_plot,
                     ncol = 2, nrow = 1,common.legend = T)
pdb_stats
ggexport(pdb_stats,filename = "pdb_stats.jpeg",width = 3000,height = 1500,
         pointsize = 100,res = 300)

#setwd("C:/Users/avagiona/Desktop/paper_triplets_figures_2025")
#ggexport(pdb_stats,filename = "pdb_stats.pdf",width = 3000,height = 1500,
#         pointsize = 100,res = 300)

# ====================================================
# SETUP for randomization and unique common 
# ====================================================
library(dplyr)
library(data.table)
library(writexl)

# Set working directory
setwd("~/R-4.2.1/predict_complexes_ML")

# Load feature-enriched ML dataset
load("ML_data_new_features.RData")
ML_data <- as.data.table(ML_data_new_features)[, -1, with = FALSE]

# Load structurally supported open triangles
#load("open_triangles.RData")  # or read.csv() if saved differently
open_triangles <- as.data.table(open_triangles)

# ====================================================
# STEP 1: Randomize V1 and V2 to avoid ordering bias
# ====================================================
set.seed(42)
ML_data[, c("V1", "V2") := {
  swap <- runif(.N) < 0.5
  v1_new <- ifelse(swap, V2, V1)
  v2_new <- ifelse(swap, V1, V2)
  list(v1_new, v2_new)
}]

# ====================================================
# STEP 2: Create canonical triplet IDs
# ====================================================
ML_data[, triplet_id := paste(common, pmin(V1, V2), pmax(V1, V2), sep = "__")]
open_triangles[, triplet_id := paste(A, pmin(B, C), pmax(B, C), sep = "__")]

# ====================================================
# STEP 3: Annotate positive (cooperative) cases
# ====================================================
ML_data[, positive := as.integer(triplet_id %in% open_triangles$triplet_id)]

# Check counts
table(ML_data$positive)

# ====================================================
# STEP 4: Generate Supplementary Table (cooperative cases)
# ====================================================
# Join structural info
cooperative_struct <- merge(
  ML_data[positive == 1, .(triplet_id, common, V1, V2)],
  open_triangles[, .(triplet_id, PDB_ID, AB_region, AC_region)],
  by = "triplet_id"
)

# Format for output
supp_table <- cooperative_struct[, .(
  PDB_ID,
  Common = common,
  V1,
  V2,
  AB_region,
  AC_region
)]

#------------------- select only one unique common interactor---------------
# Convert to data.table if needed
supp_table <- as.data.table(supp_table)

# Create triplet_id to identify unique triplets
supp_table[, triplet_id := paste(Common, pmin(V1, V2), pmax(V1, V2), sep = "__")]

# Randomly select one triplet per Common protein
set.seed(42)
selected_triplets <- supp_table[, .SD[sample(.N, 1)], by = Common]

#order correctly the interactin regions
# Function to reorder region string to put Common first
flip_region_if_needed <- function(region_string, common) {
  parts <- strsplit(region_string, " \\| ")[[1]]
  if (startsWith(parts[1], common)) {
    return(region_string)
  } else {
    return(paste(parts[2], parts[1], sep = " | "))
  }
}
#add triplet ID 
# Apply to region columns
selected_triplets[, Common_V1_region := mapply(flip_region_if_needed, AB_region, Common)]
selected_triplets[, Common_V2_region := mapply(flip_region_if_needed, AC_region, Common)]

selected_triplets[, Triplet_ID := .I]

#add triplet ID and final formating
supp_table_final <- selected_triplets[, .(
  Triplet_ID,
  PDB_ID,
  Common,
  V1,
  V2,
  Common_V1_region,
  Common_V2_region
)]

#export
#writexl::write_xlsx(supp_table_final, "Supplementary_Cooperative_Triplets_Cleaned.xlsx")

# ====================================================
# Assign 'cooperative' label in ML_data
#         using only 1 de-biased positive per Common protein
# ====================================================
ML_data[, positive := NULL]
# Use only selected cooperative triplets (1 per common)
ML_data[, cooperative := as.integer(triplet_id %in% selected_triplets$triplet_id)]
rm(ML_data_new_features)
# ====================================================================
# Machine Learning part
# ====================================================================
# ========================
# SETUP
# ========================
library(data.table)
library(dplyr)
library(tidyr)
library(caret)
library(ROSE)
library(pROC)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(forcats)
library(randomForest)
require(caTools)

# Parallel setup
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# ========================
# DATA PREP
# ========================
df_ml <- ML_data %>%
  dplyr::select(4:45, cooperative) %>%
  tidyr::drop_na()

df_ml$cooperative <- as.factor(df_ml$cooperative)

set.seed(42)
sample_df_ml <- sample.split(df_ml$cooperative, SplitRatio = 0.70)
train_df_ml <- subset(df_ml, sample_df_ml == TRUE)
test_df_ml <- subset(df_ml, sample_df_ml == FALSE)

# Undersample to balance
data_balanced_under <- ovun.sample(cooperative ~ ., data = train_df_ml, method = "under", N = 296)$data
data_balanced_under_replace <- data_balanced_under %>%
  mutate(cooperative = ifelse(cooperative == 0, "No", "Yes"))
test_df_ml_replace <- test_df_ml %>%
  mutate(cooperative = ifelse(cooperative == 0, "No", "Yes"))

# ========================
# TRAIN MODELS IN PARALLEL
# ========================
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 10, 
                     returnResamp = "final", savePredictions = TRUE)

set.seed(42)
rf_under <- train(cooperative ~ ., data = data_balanced_under, method = 'rf', metric = "Accuracy", trControl = ctrl)
svm_under <- train(cooperative ~ ., data = data_balanced_under_replace, method = "svmLinear", metric = "Accuracy", trControl = ctrl)
lg_under <- train(cooperative ~ ., data = data_balanced_under, method = "glmnet", metric = "Accuracy", trControl = ctrl)
dt_under <- train(cooperative ~ ., data = data_balanced_under, method = "rpart", metric = "Accuracy", trControl = ctrl)
knn_under <- train(cooperative ~ ., data = data_balanced_under, method = "knn", metric = "Accuracy", trControl = ctrl)

# ========================
# EVALUATION
# ========================
# Define thinning function (if not already defined)
thin_roc <- function(roc_obj, step = 10000) {
  if (length(roc_obj$sensitivities) > step) {
    indices <- seq(1, length(roc_obj$sensitivities), by = step)
    roc_obj$sensitivities <- roc_obj$sensitivities[indices]
    roc_obj$specificities <- roc_obj$specificities[indices]
    roc_obj$thresholds <- roc_obj$thresholds[indices]
  }
  return(roc_obj)
}

# Updated evaluation function to optionally return the full ROC
evaluate_model <- function(model, test_data, label_col) {
  pred_class <- predict(model, newdata = test_data[, !label_col, with = FALSE], type = "raw")
  pred_prob <- tryCatch(predict(model, newdata = test_data[, !label_col, with = FALSE], type = "prob"), error = function(e) NULL)
  
  cm <- confusionMatrix(pred_class, as.factor(test_data[[label_col]]), mode = "everything")
  roc_obj <- if (!is.null(pred_prob)) {
    roc(as.numeric(test_data[[label_col]]), as.numeric(pred_prob[[2]]))
  } else {
    roc(as.numeric(test_data[[label_col]]), as.numeric(pred_class))
  }
  
  roc_obj <- thin_roc(roc_obj)  # Thin right away to save memory
  
  list(confusion = cm, auc = roc_obj$auc, roc = roc_obj)
}

# Evaluate one-by-one, thinning immediately
results_rf <- evaluate_model(rf_under, test_df_ml, "cooperative")
results_svm <- evaluate_model(svm_under, test_df_ml_replace, "cooperative")
results_glm <- evaluate_model(lg_under, test_df_ml, "cooperative")
results_dt <- evaluate_model(dt_under, test_df_ml, "cooperative")
results_knn <- evaluate_model(knn_under, test_df_ml, "cooperative")

# Combine into one list after everything is safe
results <- list(
  RF = results_rf,
  SVM = results_svm,
  GLM = results_glm,
  DT = results_dt,
  KNN = results_knn
)


# ========================
# FEATURE IMPORTANCE
# ========================
i_scores <- varImp(rf_under, scale = FALSE)
importance_df <- as.data.frame(i_scores$importance)
importance_df$feature <- rownames(importance_df)
importance_df_or <- importance_df[order(-importance_df$Overall), ]

feature_plot <- ggplot(importance_df_or, aes(x = fct_reorder(feature, Overall), y = Overall)) +
  geom_col(fill = "#44AA99") +
  coord_flip() +
  labs(x = "Feature", y = "Importance") +
  theme_bw()

ggexport(feature_plot, filename = "feature_importance_cooperative.png", width = 4000, height = 2000, pointsize = 100, res = 250)

# ========================
# CLEANUP PARALLEL
# ========================
stopCluster(cl)
registerDoSEQ()

#====================================================
# Machine learning part
#====================================================
#set pathway
setwd("~/R-4.2.1/predict_complexes_ML")
# Step 7: Load HIPPIE interactions and features for modeling
load("ML_data_new_features.RData")
ML_data<-ML_data_new_features%>%
 dplyr::select(-1)
#data.table::fwrite(ML_data[, .(common, V1, V2)], "ML_data_triplets_2025.csv")
#data.table::fwrite(open_triangles[, .(A, B, C)], "open_triangles_triplets_2025.csv")
#=================add the positive caese to the ML_data
library("data.table")
labels <- fread("ML_data_with_labels.csv")
# Make sure both are data.tables
ML_data <- as.data.table(ML_data)
labels <- as.data.table(labels)

# Merge the 'positive' column into ML_data using common + V1 + V2
ML_data <- merge(
  ML_data,
  labels[, .(common, V1, V2, positive)],
  by = c("common", "V1", "V2"),
  all.x = TRUE
)

#===========randomize V1 and V2
set.seed(42)  # for reproducibility

#randomize V1 and V2 
ML_data[, c("V1", "V2") := {
  swap <- runif(.N) < 0.5  # TRUE for half of the rows
  v1_new <- ifelse(swap, V2, V1)
  v2_new <- ifelse(swap, V1, V2)
  list(v1_new, v2_new)
}]

#=======================================================================
#cooperative cases with structural annotations/ Supplementary table 4(?)
#=======================================================================
# Ensure data.table is loaded

ML_data <- as.data.table(ML_data)
open_triangles <- as.data.table(open_triangles)
# Step 1: Create triplet_id in both datasets
ML_data[, triplet_id := paste(common, pmin(V1, V2), pmax(V1, V2), sep = "__")]
open_triangles[, triplet_id := paste(A, pmin(B, C), pmax(B, C), sep = "__")]

# Step 2: Merge cooperative cases with structural annotations
cooperative_struct <- merge(
  ML_data[positive == 1, .(triplet_id, common, V1, V2)],
  open_triangles[, .(triplet_id, PDB_ID, AB_region, AC_region)],
  by = "triplet_id"
)

# Step 3: Reorder columns
supp_table <- cooperative_struct[, .(
  PDB_ID,
  Common = common,
  V1,
  V2,
  AB_region,
  AC_region
)]

# Step 4: Export
#writexl::write_xlsx(supp_table, "Supplementary_Cooperative_Triplets.xlsx")

























