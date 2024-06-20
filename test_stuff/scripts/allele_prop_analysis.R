# ---
# title: "proportions_distr"
# author: "Andres Florian"
# date: "2024-05-27"
# output: html_document
# ---

# test
# Read the data from the TSV file, skipping the first two lines
data <- read.table(
  "/Volumes/PiconCossio/biocomp_tools/plodinf/test_stuff/test_data/diploid_sort.bam_all_AfRef_allele_prop_m15_d25.tsv", # nolint
  header = TRUE, skip = 2
)

# Filter data for rows where the third column equals "2"
filtered_data <- data[data[, 4] == 2, ]

# Extract and combine the 4th and 5th column into a unique list item
# unique_item <- list(paste(filtered_data[, 4], filtered_data[, 5], sep = ",")) # nolint

unique_item <- list(filtered_data[, 7])

# Function_def


library(ggplot2)

# Function to read and process data from a single file
read_and_process_file <- function(file_path) {
  # Read the data from the TSV file, skipping the first two lines
  data <- read.table(file_path, header = TRUE, skip = 2)

  # Filter data for rows where the third column equals "2"
  filtered_data <- data[data[, 4] == 2, ]

  # Extract and combine the 4th and 5th column into a unique list item
  # unique_item <- list(paste(filtered_data[, 4], filtered_data[, 5], sep = ","))

  unique_item <- list(filtered_data[, 7])
  return(unique_item)
}


load_and_plot_wrapper <- function(data_path) {
  # List to store unique items from each file
  all_items <- list()
  file_names <- c()
  # Loop through all files ending with ".tsv" in the specified path
  for (file in list.files(data_path, pattern = ".tsv")) {
    # Read and process the current file

    # print(strsplit(file,"_")[[1]][1])
    current_item <- read_and_process_file(paste0(data_path, file))

    file_names <- c(file_names, strsplit(file, "_")[[1]][1])
    # Add the unique item from the current file to the list
    all_items <- c(all_items, current_item)
  }

  names(all_items) <- file_names


  # Loop through each unique item and create a ggplot2 histogram
  for (item_name in names(all_items)) {
    # Split the item back into separate values (assuming comma separation)
    # print(item_name)
    item <- all_items[[item_name]]
    values <- unlist(strsplit(item, ","))

    # Create the ggplot2 histogram
    ploty <- ggplot(aes(x = as.numeric(values)), data = data.frame(values)) +
      geom_histogram(binwidth = 0.03, color = "#e9ecef", alpha = 0.7, position = "identity") +
      labs(title = paste("Histogram of", item_name), x = item_name, y = "Count") +
      xlim(0, 1)
    theme_bw()

    # Print the histogram plot (consider saving to a file if needed)
    #print(ggsave(filename = NULL))  # To display the plot in the console
    print(ploty)
  }
}


# Running


data_path <- "/Volumes/PiconCossio/biocomp_tools/plodinf/test_stuff/test_data/"

load_and_plot_wrapper(data_path)


# Tests QQ plot


# all_items <- list()
# file_names <- c()
# # Loop through all files ending with ".tsv" in the specified path
# for (file in list.files(data_path, pattern = ".tsv")) {
#   # Read and process the current file

#   # print(strsplit(file,"_")[[1]][1])
#   current_item <- read_and_process_file(paste0(data_path, file))
#   # print(file)

#   splitted_file_name <- strsplit(file, "_")

#   if (splitted_file_name[[1]][2] %in% c("1", "2", "3")) {
#     file_name <- paste(splitted_file_name[[1]][1], splitted_file_name[[1]][2], sep = "_")
#   } else {
#     file_name <- splitted_file_name[[1]][1]
#   }


#   file_names <- c(file_names, file_name)
#   # Add the unique item from the current file to the list
#   all_items <- c(all_items, current_item)
# }

# names(all_items) <- file_names

# library(ggplot2)
# ploty <- ggplot()

# # Loop through each unique item and create a ggplot2 histogram
# for (item_name in names(all_items)) {
#   # Split the item back into separate values (assuming comma separation)
#   # print(item_name)
#   item <- all_items[[item_name]]
#   values <- unlist(strsplit(item, ","))

#   # Create the ggplot2 histogram
#   ploty <- ploty + stat_qq(aes(sample = as.numeric(values)), colour = "green")
#   # Print the histogram plot (consider saving to a file if needed)
#   # print(ggsave(filename = NULL))  # To display the plot in the console
#   print(ploty)
# }

# acc_1 <- "VA080"
# acc_2 <- "A001"
# ggplot() +
#   stat_qq(aes(sample = as.numeric(unlist(strsplit(all_items[[acc_1]], ",")))), colour = "green") +
#   stat_qq(aes(sample = as.numeric(unlist(strsplit(all_items[[acc_2]], ",")))), colour = "green")
# # stat_qq(aes(sample = unlist(split_P2_2_all[,c("prop1","prop2")])), colour = "red")
# # geom_abline(aes(slope = 1, intercept = 0), linetype = 2)


# # Function to convert elements of a vector to numeric
# convert_to_numeric <- function(x) {
#   if (is.character(x)) {
#     # Split string by character and convert to numeric
#     as.numeric(unlist(strsplit(x, ",")))
#   } else {
#     # If not character, return as is
#     x
#   }
# }

# # Apply the function to each element of the list
# all_items_mod <- lapply(all_items, convert_to_numeric)


# library(ggplot2)

# # Function to plot qqplot of a vector
# plot_qq <- function(data) {
#   ggplot(aes(sample = qnorm(seq_len(length(data))))) +
#     geom_qq(aes(y = data), color = "blue") +
#     stat_qq_line(color = "red") +
#     labs(title = "QQ Plot", x = "Quantile (Standard Normal)", y = "Quantile (Data)") +
#     theme_bw()
# }

# # Input list of numerical vectors
# data_list <- list(rnorm(10), runif(15), rpois(20)) # Replace with your actual data

# # Create a ggplot with multiple facets for each vector in the list
# ggplot() +
#   facet_wrap(~.names, nrow = 1) + # Wrap plots in a single row
#   plot_qq(data = .x) # Call the plot_qq function for each data in the list

# # Print the plot
# print(ggsave(filename = NULL)) # To display the plot in the console



# library(ggplot2)

# # Function to create a qqplot for a vector
# create_qqplot <- function(data) {
#   # Generate quantiles of a standard normal distribution
#   qnorm_values <- qnorm(seq(0, 1, length = length(data)))

#   # Create the ggplot with qqplot geom
#   ggplot() +
#     geom_qq(aes(sample = data), color = "blue") +
#     labs(title = "QQ Plot", x = "Standard Normal Quantiles", y = "Data Values") +
#     theme_void()
# }

# # Loop through each vector in the list and create qqplot
# for (item in all_items_mod) {
#   # Create the qqplot for the current vector
#   qqplot <- create_qqplot(item)

#   # Add the qqplot to a new layer if not the first iteration
#   if (length(all_items_mod) > 1) {
#     qqplot <- qqplot + geom_qq(aes(sample = item), color = "red", linetype = "dashed")
#   }

#   # Print the combined qqplot (consider saving to a file if needed)
#   print(qqplot)
# }




# library(ggplot2)

# # Function to create a qqplot for a vector
# create_qqplot <- function(data, color = "blue", linetype = solid) {
#   # Generate quantiles of a standard normal distribution
#   qnorm_values <- qnorm(seq(0, 1, length = length(data)))
#   qnorm_data <- sort(data)
#   # Create the qqplot geom with specified color and linetype
#   geom_point(aes(x = qnorm_values, y = qnorm_data), color = color)
# }

# # Empty ggplot object to accumulate qqplots
# all_qqplots <- ggplot()

# # Loop through each vector in the list and add qqplot to the main plot
# for (item in all_items_mod) {
#   # Create the qqplot for the current vector with distinct color/linetype
#   current_qqplot <- create_qqplot(item, color = ifelse(length(all_items_mod) > 1, "red", "blue"), linetype = ifelse(length(all_items_mod) > 1, "dashed", "solid"))

#   # Add the qqplot to the main plot using `+` operator
#   all_qqplots <- all_qqplots + current_qqplot
# }

# # Set common aesthetics and theme for the combined plot
# all_qqplots <- all_qqplots +
#   labs(title = "QQ Plot (Overlaid)", x = "Standard Normal Quantiles", y = "Data Values") +
#   theme_bw()

# # Print the combined qqplot (consider saving to a file if needed)
# print(all_qqplots)




# qnorm_values <- qnorm(seq(0, 1, length = length(all_items_mod[["A001"]])))
# # qnorm_data = sort(ecdf(all_items_mod[["A001"]])(all_items_mod[["A001"]]))
# qnorm_data <- sort(all_items_mod[["A001"]])

# ggplot() +
#   geom_point(aes(x = qnorm_values, y = qnorm_data), color = "blue")




# qnorm_values <- qnorm(seq(0, 1, length = length(all_items_mod[["A001"]])))
# qnorm_data <- sort(ecdf(all_items_mod[["A001"]])(all_items_mod[["A001"]]))
# # qnorm_data = sort(all_items_mod[["A001"]])

# ggplot() +
#   geom_point(aes(x = qnorm_values, y = qnorm_data), color = "blue")



# # Tests pariwise KS_test



# name_file <- "/home/andresflorian/University/Cebolla/ADMIXTURE_test/custom_script_AD3_DP6_filter/name_conversion.tsv"
# names <- read.table(name_file, header = TRUE)
# rownames(names) <- names$Short_name




# data_path <- "/home/andresflorian/University/Cebolla/internship/ploidy/allelic_proportions_libraries/run02/"

# all_items <- list()
# file_names <- c()
# # Loop through all files ending with ".tsv" in the specified path
# for (file in list.files(data_path, pattern = ".tsv")) {
#   # Read and process the current file

#   # print(strsplit(file,"_")[[1]][1])
#   current_item <- read_and_process_file(paste0(data_path, file))
#   # print(file)

#   splitted_file_name <- strsplit(file, "_")

#   if (splitted_file_name[[1]][2] %in% c("1", "2", "3")) {
#     file_name <- paste(splitted_file_name[[1]][1], splitted_file_name[[1]][2], sep = "_")
#   } else {
#     file_name <- splitted_file_name[[1]][1]
#   }


#   file_names <- c(file_names, file_name)
#   # Add the unique item from the current file to the list
#   all_items <- c(all_items, current_item)
# }

# names(all_items) <- file_names


# # Function to convert elements of a vector to numeric
# convert_to_numeric <- function(x) {
#   if (is.character(x)) {
#     # Split string by character and convert to numeric
#     as.numeric(unlist(strsplit(x, ",")))
#   } else {
#     # If not character, return as is
#     x
#   }
# }

# # Apply the function to each element of the list
# all_items_mod <- lapply(all_items, convert_to_numeric)

# rm(all_items)

# library(dgof)

# # Function to perform ks test and return p-value
# ks_test_results <- function(x, y) {
#   # Perform Kolmogorov-Smirnov test and extract p-value
#   ks_test <- ks.test(x = x, y = ecdf(y))
#   results <- c(ks_test$p.value, ks_test$statistic[["D"]])

#   return(results)
# }

# # Get all unique vector names from the list
# vector_names <- names(all_items_mod)

# # Initialize empty dataframe to store p-values
# p_value_matrix <- matrix(nrow = length(vector_names), ncol = length(vector_names), data = NA)
# D_stats_matrix <- matrix(nrow = length(vector_names), ncol = length(vector_names), data = NA)

# # Fill the upper triangle of the matrix with p-values
# for (i in 1:length(vector_names)) {
#   # for (j in (i + 1):length(vector_names)) {
#   for (j in 1:length(vector_names)) {
#     # Get data vectors by name
#     x <- all_items_mod[[vector_names[i]]]
#     y <- all_items_mod[[vector_names[j]]]
#     # print(c(i,j))
#     # print(c(vector_names[i],vector_names[j]))
#     # Calculate p-value using the function and store

#     results <- ks_test_results(x, y)

#     p_value_matrix[i, j] <- results[1]
#     D_stats_matrix[i, j] <- results[2]
#   }
# }

# # Fill the lower triangle by mirroring upper triangle values (assuming symmetry)
# # You can remove this if p-values for both directions are needed
# # diag(p_value_matrix) <- 1
# # p_value_matrix <- t(p_value_matrix) + p_value_matrix

# # Set row and column names of the dataframe
# rownames(p_value_matrix) <- vector_names
# colnames(p_value_matrix) <- vector_names

# rownames(D_stats_matrix) <- vector_names
# colnames(D_stats_matrix) <- vector_names


# # Print the dataframe matrix (consider saving to a file if needed)
# # print(p_value_matrix)


# library(pheatmap)

# pheatmap::pheatmap(D_stats_matrix, annotation_row = names[, c("Population"), drop = FALSE])


# library(pheatmap)

# alpha <- 0.001

# mod_pval_matrix <- p_value_matrix
# mod_pval_matrix[mod_pval_matrix > alpha] <- 1
# mod_pval_matrix[mod_pval_matrix <= alpha] <- 0

# pheatmap::pheatmap(mod_pval_matrix)
