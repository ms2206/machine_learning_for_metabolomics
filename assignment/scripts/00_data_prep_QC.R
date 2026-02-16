# script.R
# load libraries

library(dplyr)
library(tidyr)
library(mixOmics)
library(dotenv)

################################################
################### Load .env ##################
################################################
cat("Loading environment variables from .env file...\n")
tryCatch(
  {
    dotenv::load_dot_env()
    cat("Environment variables loaded successfully.\n")
  },
  error = function(e) {
    cat("Error loading environment variables:", e$message, "\n")
    cat("Loading from hardcoded path as fallback...\n")
    dotenv::load_dot_env("/Users/mspriggs/Library/CloudStorage/OneDrive-Illumina,Inc./Documents/Applied_Bioinformatics/modules/machine_learning_for_metabolomics/assignment/scripts/.env") # nolint
    cat("Loaded .env from hardcoded path\n")
  }
)


################################################
############# Set up Global Scope ##############
################################################
cat("Setting up global variables and sourcing utility functions...\n")
WORKING_DIR <- Sys.getenv("WORKING_DIR")

# define file paths
UTILS_FILEPATH <- file.path(WORKING_DIR, "scripts", "00_data_prep_utils.R")
DATA_DIR <- file.path(WORKING_DIR, "Assignment_Data")
OUTPUT_DIR <- file.path(WORKING_DIR, "outputs")
PLOTS_DIR <- file.path(OUTPUT_DIR, "plots")

date_str <- format(Sys.Date(), "%Y-%m-%d")



# source utility functions
source(UTILS_FILEPATH)
print(paste("Utility functions sourced from:", UTILS_FILEPATH))

################################################
################### Main #######################
################################################


main <- function(script_dir) {
  cat("\n")
  cat("##################################################", "\n")
  cat("Running data preparation and quality control...\n")
  cat("##################################################", "\n\n")

  ################################################
  ###### data prep & exploration [15marks] #######
  ################################################

  # load data
  ca_year_2 <- load_data(file.path(DATA_DIR, "CA_Year2.csv"))
  dca_year_2 <- load_data(file.path(DATA_DIR, "DCA_Year2.csv"))
  dca_year_1 <- load_data(file.path(DATA_DIR, "DCA_Year1.csv"))

  # check dimensions
  cat("Checking dimensions of the datasets...\n")
  cat("CA Year 2 dimensions:", dim(ca_year_2), "\n")
  cat("DCA Year 2 dimensions:", dim(dca_year_2), "\n")
  cat("DCA Year 1 dimensions:", dim(dca_year_1), "\n\n")

  # check for missing values
  cat("Checking for missing values in the datasets...\n")
  cat("Missing values in CA Year 2:", sum(is.na(ca_year_2)), "\n")
  cat("Missing values in DCA Year 2:", sum(is.na(dca_year_2)), "\n")
  cat("Missing values in DCA Year 1:", sum(is.na(dca_year_1)), "\n\n")

  # infer mean for na values in DCA Year 2
  cat("Inferring mean values for missing entries in DCA Year 2...\n")
  dca_year_2 <- infer_mean_for_nan(dca_year_2)
  cat("Missing values in DCA Year 2:", sum(is.na(dca_year_2)), "\n\n")

  # remove missing values from DCA Year 1
  cat("Removing rows with missing values from DCA Year 1...\n")
  dca_year_1 <- dca_year_1 %>% drop_na()
  cat("Missing values in DCA Year 1:", sum(is.na(dca_year_1)), "\n\n")

  # extract features and labels
  ca_year_2_samples <- extract_sample_names(ca_year_2)
  rownames(ca_year_2) <- ca_year_2_samples
  ca_year_2 <- ca_year_2[, -1]

  dca_year_2_samples <- extract_sample_names(dca_year_2)
  rownames(dca_year_2) <- dca_year_2_samples
  dca_year_2 <- dca_year_2[, -1]

  dca_year_1_samples <- extract_sample_names(dca_year_1)
  rownames(dca_year_1) <- dca_year_1_samples
  dca_year_1 <- dca_year_1[, -1]

  # convert storage.week and ripening.stage to factors
  ca_year_2 <- mutate_df_factors(ca_year_2)
  dca_year_2 <- mutate_df_factors(dca_year_2)
  dca_year_1 <- mutate_df_factors(dca_year_1)

  # principal component analysis (PCA) for data exploration

  #####################################
  ############### PCA ################
  #####################################


  #####################################
  ######## Scaled and Centered ########
  #####################################
  cat("Performing PCA for data exploration...\n")
  cat("PCA with scaling and centering...\n")
  ca_year_2_pca <- pca_wrapper(ca_year_2, scale = TRUE, center = TRUE)
  dca_year_1_pca <- pca_wrapper(dca_year_1, scale = TRUE, center = TRUE)
  dca_year_2_pca <- pca_wrapper(dca_year_2, scale = TRUE, center = TRUE)

  # plots for PCA

  png(filename = file.path(PLOTS_DIR, paste0(date_str, "ca_year_2_pca_pca_scree_plots.png")), width = 1200, height = 800)
  plot(ca_year_2_pca, type = "screeplot", main = "Scree Plot for CA Year 2 PCA | Scaled & Centered")
  dev.off()

  png(filename = file.path(PLOTS_DIR, paste0(date_str, "dca_year_1_pca_scree_plots.png")), width = 1200, height = 800)
  plot(dca_year_1_pca, type = "screeplot", main = "Scree Plot for DCA Year 1 PCA | Scaled & Centered")
  dev.off()

  png(filename = file.path(PLOTS_DIR, paste0(date_str, "dca_year_2_pca_scree_plots.png")), width = 1200, height = 800)
  plot(dca_year_2_pca, type = "screeplot", main = "Scree Plot for DCA Year 2 PCA | Scaled & Centered")
  dev.off()

  #############################
  ######## Scaled only ########
  #############################

  cat("PCA with scaling only (no centering)...\n")
  ca_year_2_pca_scaled <- pca_wrapper(ca_year_2, scale = TRUE, center = FALSE)
  dca_year_1_pca_scaled <- pca_wrapper(dca_year_1, scale = TRUE, center = FALSE)
  dca_year_2_pca_scaled <- pca_wrapper(dca_year_2, scale = TRUE, center = FALSE)

  plot(ca_year_2_pca_scaled, type = "screeplot", main = "Scree Plot for CA Year 2 PCA | Scaled Only")
  plot(dca_year_1_pca_scaled, type = "screeplot", main = "Scree Plot for DCA Year 1 PCA | Scaled Only")
  plot(dca_year_2_pca_scaled, type = "screeplot", main = "Scree Plot for DCA Year 2 PCA | Scaled Only")


  ###############################
  ######## Centered only ########
  ###############################

  cat("PCA with centering only (no scaling)...\n")
  ca_year_2_pca_centered <- pca_wrapper(ca_year_2, scale = FALSE, center = TRUE)
  dca_year_1_pca_centered <- pca_wrapper(dca_year_1, scale = FALSE, center = TRUE)
  dca_year_2_pca_centered <- pca_wrapper(dca_year_2, scale = FALSE, center = TRUE)

  plot(ca_year_2_pca_centered, type = "screeplot", main = "Scree Plot for CA Year 2 PCA | Centered Only")
  plot(dca_year_1_pca_centered, type = "screeplot", main = "Scree Plot for DCA Year 1 PCA | Centered Only")
  plot(dca_year_2_pca_centered, type = "screeplot", main = "Scree Plot for DCA Year 2 PCA | Centered Only")


  #########################################
  ######## Unscaled and Uncentered ########
  #########################################

  cat("PCA without scaling or centering...\n")
  ca_year_2_pca_none <- pca_wrapper(ca_year_2, scale = FALSE, center = FALSE)
  dca_year_1_pca_none <- pca_wrapper(dca_year_1, scale = FALSE, center = FALSE)
  dca_year_2_pca_none <- pca_wrapper(dca_year_2, scale = FALSE, center = FALSE)


  plot(ca_year_2_pca_none, type = "screeplot", main = "Scree Plot for CA Year 2 PCA | Unscaled & Uncentered")
  plot(dca_year_1_pca_none, type = "screeplot", main = "Scree Plot for DCA Year 1 PCA | Unscaled & Uncentered")
  plot(dca_year_2_pca_none, type = "screeplot", main = "Scree Plot for DCA Year 2 PCA | Unscaled & Uncentered")


  ## create a bioplots for the three datasets using the scaled and centered PCA results, colored by storage week and ripening

  cat("Creating biplots for PCA results using scaled and centered data...\n")
  # biplots for CA Year 2
  save_plot(
    biplot_wrapper(ca_year_2_pca, ca_year_2$Storage.Week, ca_year_2_samples),
    "biplot_ca_year_2_storage_week.png",
    PLOTS_DIR
  )

  save_plot(
    biplot_wrapper(ca_year_2_pca, ca_year_2$Ripening.Stage, ca_year_2_samples),
    "biplot_ca_year_2_ripening_stage.png",
    PLOTS_DIR
  )

  # biplots for DCA Year 1
  save_plot(
    biplot_wrapper(dca_year_1_pca, dca_year_1$Storage.Week, dca_year_1_samples),
    "biplot_dca_year_1_storage_week.png",
    PLOTS_DIR
  )

  save_plot(
    biplot_wrapper(dca_year_1_pca, dca_year_1$Ripening.Stage, dca_year_1_samples),
    "biplot_dca_year_1_ripening_stage.png",
    PLOTS_DIR
  )

  # biplots for DCA Year 2
  save_plot(
    biplot_wrapper(dca_year_2_pca, dca_year_2$Storage.Week, dca_year_2_samples),
    "biplot_dca_year_2_storage_week.png",
    PLOTS_DIR
  )

  save_plot(
    biplot_wrapper(dca_year_2_pca, dca_year_2$Ripening.Stage, dca_year_2_samples),
    "biplot_dca_year_2_ripening_stage.png",
    PLOTS_DIR
  )

  # Can I merge dca_year_1 and dca_year_2 and do a PCA on the combined dataset?
  cat("Merging DCA Year 1 and Year 2 datasets for combined PCA...\n")

  # Fix column name mismatch
  colnames(dca_year_1)[colnames(dca_year_1) == "X7OH.ABA"] <- "X7.OH.ABA"

  dca_combined <- rbind(dca_year_1, dca_year_2)
  dca_combined_samples <- c(dca_year_1_samples, dca_year_2_samples)
  dca_combined_pca <- pca_wrapper(dca_combined, scale = TRUE, center = TRUE)

  # add column to indicate year
  dca_combined$Year <- c(rep("Year 1", nrow(dca_year_1)), rep("Year 2", nrow(dca_year_2)))
  dca_combined$Year <- factor(dca_combined$Year)

  save_plot(
    biplot_wrapper(dca_combined_pca, dca_combined$Storage.Week, dca_combined_samples),
    "biplot_dca_combined_storage_week.png",
    PLOTS_DIR
  )

  save_plot(
    biplot_wrapper(dca_combined_pca, dca_combined$Ripening.Stage, dca_combined_samples),
    "biplot_dca_combined_ripening_stage.png",
    PLOTS_DIR
  )

  save_plot(
    biplot_wrapper(dca_combined_pca, dca_combined$Year, dca_combined_samples),
    "biplot_dca_combined_year.png",
    PLOTS_DIR
  )

  # PCA2 on the virtical axis looks to separate year, is this a batch effect?
  # Extract PC2 scores
  pc2_scores <- dca_combined_pca$x[, 2]

  # Add metadata
  combined_pc2_df <- data.frame(
    PC2 = pc2_scores,
    Storage_Week = dca_combined$Storage.Week, # Adjust column name as needed
    Year = dca_combined$Year,
    Ripening_Stage = dca_combined$Ripening.Stage
  )

  # Plot PC2 vs Storage Week

  pc2_storage_week <- ggplot(combined_pc2_df, aes(x = Storage_Week, y = PC2, color = Year)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = "PC2 vs Storage Week",
      x = "Storage Week",
      y = "PC2 Score (13.59%)"
    ) +
    theme_bw()

  save_plot(
    pc2_storage_week,
    "pc2_vs_storage_week.png",
    PLOTS_DIR
  )


  # Plot PC2 vs Ripening Stage
  pc2_ripening_stage <- ggplot(combined_pc2_df, aes(x = Ripening_Stage, y = PC2, color = Year)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = "PC2 vs Ripening Stage",
      x = "Ripening Stage",
      y = "PC2 Score (13.59%)"
    ) +
    theme_bw()

  save_plot(
    pc2_ripening_stage,
    "pc2_vs_ripening_stage.png",
    PLOTS_DIR
  )
}




# run main function
main(script_dir)
