library(dplyr)
library(tidyr)
library(mixOmics)

load_data <- function(filepath) {
    # Load the data
    data <- read.csv(filepath)
    return(data)
}

infer_mean_for_nan <- function(df) {
    # Infer mean values for NA entries in the dataframe
    # and replace them.

    nan_locations <- which(is.na(df), arr.ind = TRUE)
    for (i in seq_len(nrow(nan_locations))) {
        row <- nan_locations[i, "row"]
        col <- nan_locations[i, "col"]
        mean_value <- mean(df[[col]], na.rm = TRUE)
        df[row, col] <- mean_value
    }
    return(df)
}

extract_sample_names <- function(df) {
    # Extract sample names from the dataframe
    row_names <- df[, 1]
    return(row_names)
}

mutate_df_factors <- function(df) {
    # Convert Storage.Week and ripening.stage to factors
    if ("Storage.Week" %in% colnames(df)) {
        df$Storage.Week <- factor(df$Storage.Week)
    }
    if ("Ripening.Stage" %in% colnames(df)) {
        df$Ripening.Stage <- factor(df$Ripening.Stage)
    }
    return(df)
}

pca_wrapper <- function(df, scale, center) {
    # Perform PCA and return the results

    # Select only numeric columns for PCA
    X <- df[, sapply(df, is.numeric), drop = FALSE]

    # Perform PCA using mixOmics
    pca_result <- pca(X, ncomp = 6, scale = scale, center = center)
    return(pca_result)
}

biplot_wrapper <- function(pca_result, group_var, sample_names) {
    # Create a biplot for the PCA results, colored by the specified group variable

    Var1 <- 100 * (pca_result$prop_expl_var$X[1])
    Var1 <- round(Var1, digits = 2)
    Var2 <- 100 * (pca_result$prop_expl_var$X[2])
    Var2 <- round(Var2, digits = 2)

    biplot(pca_result,
        ind.names = sample_names,
        group = group_var,
        style = "ggplot2"
    ) + ggtitle(paste(deparse(substitute(pca_result)), "grouped by", deparse(substitute(group_var)))) +
        xlab(paste("PC1 (", Var1, "%)", sep = "")) +
        ylab(paste("PC2 (", Var2, "%)", sep = ""))
}

save_plot <- function(plot_obj, filename, output_dir) {
    # Save the plot to the specified filename with a date appended to it

    # Append date to filename
    date_str <- format(Sys.Date(), "%Y-%m-%d")
    filename_with_date <- paste0(date_str, "_", filename)

    # Save the plot
    filepath <- file.path(output_dir, filename_with_date)
    ggplot2::ggsave(filepath, plot = plot_obj, width = 10, height = 8, dpi = 300) # nolint
    cat("Plot saved to:", filepath, "\n")
}

# Function to swap train/test data
switch_training_test_data <- function(train_data, test_data) {
    # Swap training and test data
    swapped_train <- test_data
    swapped_test <- train_data

    # Return as a list
    return(list(
        training_data = swapped_train,
        test_data = swapped_test
    ))
}

# Function to plot confusion matrix as a heatmap
plot_confusion_matrix <- function(conf_matrix, title = "Confusion Matrix") {
    # Convert confusion matrix to data frame
    cm_df <- as.data.frame(conf_matrix$table)
    colnames(cm_df) <- c("Reference", "Predicted", "Count")
    
    # Create and return ggplot
    plot <- ggplot2::ggplot(cm_df, ggplot2::aes(x = Reference, y = Predicted, fill = Count)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = Count), size = 3) +
        ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = title, x = "Predicted", y = "Reference")
    
    return(plot)
}
