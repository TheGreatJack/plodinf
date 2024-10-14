

# Simulation Functions:

## Read depth generation

pois_read_depth_distr_generator <- function(min_depth, max_depth, average, site_number){
  
  #possible read depths
  read_depth_pos <- seq(min_depth,max_depth)
  #Get probabilities for each read depth from a possion distribution
  pois_distr <- dpois(read_depth_pos,average_pois)
  
  #Sample the final read depth distribution
  read_depth_disrt_pois <- sample(read_depth_pos,number_of_sites,replace = TRUE, prob = pois_distr)
  
  return(read_depth_disrt_pois)
}


## Simulation 2n , 3n , 4n

# Main function to generate simulated read proportions 

snp_simulator_2n_4n <- function(current_ploidy, read_depth_table,central_prop) {
  
  read_depth_table <- as.data.frame(table(read_depth_table))
  read_depth_table$read_depth_table <- as.numeric(as.character(read_depth_table$read_depth_table))
  
  # Define biallelic frequencies for different ploidies
  biallelic_freqs_diploid <- list(c(0.5, 0.5))
  biallelic_freqs_triploid <- list(c(0.33, 0.66))
  biallelic_freqs_tetraploid <- list(c(0.5,0.5),c(0.25, 0.75))
  
  
  # Get corresponding biallelic frequencies
  if (current_ploidy == 2) {
    biallelic_freqs <- biallelic_freqs_diploid
    all_props <- snp_gen_non_tetraploid(biallelic_freqs, read_depth_table)
  } else if (current_ploidy == 3) {
    biallelic_freqs <- biallelic_freqs_triploid
    all_props <- snp_gen_non_tetraploid(biallelic_freqs, read_depth_table)
  } else if (current_ploidy == 4){
    biallelic_freqs <- biallelic_freqs_tetraploid
    all_props <- snp_gen_tetraploid(biallelic_freqs, read_depth_table,central_prop)
  }
  colnames(all_props) <- c("read_depth","reads_a","reads_b","prop_a","prop_b")
  all_props$prop_diff <- all_props$prop_a - all_props$prop_b
  return(all_props)
}

# function for 2n and 3n

snp_gen_non_tetraploid <- function(biallelic_freqs, read_depth_table) {
  # Initialize empty matrix to store allele proportions
  #all_props <- matrix(nrow = n_sites, ncol = 5)
  all_props <- data.frame()
  
  # Simulate for each SNP site
  for (i in 1:nrow(read_depth_table)) {
    read_depth <- read_depth_table[i,1]
    read_depth_freq <- read_depth_table[i,2]
    #print(biallelic_freqs[[1]])
    
    for (site in 1:read_depth_freq) {
      # Sample alleles based on frequencies
      alleles <- sample(c("A", "B"), size = read_depth, prob = biallelic_freqs[[1]], replace = TRUE)
      
      alleles_a <- sum(alleles == "A")
      alleles_b <- sum(alleles == "B")
      
      # Estimate allele proportions
      prop_a <- alleles_a / read_depth
      prop_b <- alleles_b / read_depth  
      
      all_props <- rbind(all_props,c(read_depth,alleles_a,alleles_b,prop_a, prop_b))
    }
  }
  # Return allele proportions
  return(all_props)
}


# function for tetraploids

snp_gen_tetraploid <- function(biallelic_freqs, read_depth_table, central_prop) {
  # Initialize empty matrix to store allele proportions
  #all_props <- matrix(nrow = n_sites, ncol = 5)
  
  out_ward_prop <- (1-central_prop)
  all_props <- data.frame()
  
  # Simulate for each SNP site
  for (i in 1:nrow(read_depth_table)) {
    read_depth <- read_depth_table[i,1]
    read_depth_freq <- read_depth_table[i,2]
    
    biallelic_type_tbl <- table(sample(c("c", "o"), size = read_depth_freq, prob = c(central_prop,out_ward_prop), replace = TRUE))
    
    read_depth_c <- biallelic_type_tbl["c"]
    read_depth_o <- biallelic_type_tbl["o"]
    #print(biallelic_freqs[[1]])
    
    if (!is.na(read_depth_c)) {
      for (site in 1:read_depth_c) {
        # Sample alleles based on frequencies
        alleles <- sample(c("A", "B"), size = read_depth, prob = biallelic_freqs[[1]], replace = TRUE)
        
        alleles_a <- sum(alleles == "A")
        alleles_b <- sum(alleles == "B")
        
        # Estimate allele proportions
        prop_a <- alleles_a / read_depth
        prop_b <- alleles_b / read_depth  
        
        all_props <- rbind(all_props,c(read_depth,alleles_a,alleles_b,prop_a, prop_b))
      }
    }
    #print(!is.na(read_depth_o))
    if (!is.na(read_depth_o)) {
      for (site in 1:read_depth_o) {
        # Sample alleles based on frequencies
        alleles <- sample(c("A", "B"), size = read_depth, prob = biallelic_freqs[[2]], replace = TRUE)
        
        alleles_a <- sum(alleles == "A")
        alleles_b <- sum(alleles == "B")
        
        # Estimate allele proportions
        prop_a <- alleles_a / read_depth
        prop_b <- alleles_b / read_depth  
        
        all_props <- rbind(all_props,c(read_depth,alleles_a,alleles_b,prop_a, prop_b))
      }
    }
  }
  # Return allele proportions
  return(all_props)
}

#Function to parse simulation of read proportions into a table of probabilities per read proportion
# proportions of 0 and 1 (equivalent to -1 and 1) are excluded. as they provide no info

simulation_table_gen <- function(current_ploidy, read_depth_distrib, central_prop){
  
  simulation_table <- snp_simulator_2n_4n(current_ploidy, read_depth_distrib,central_prop)
  
  inversion_vector <- sample(c(-1,1), size = length(simulation_table$prop_diff),replace = TRUE)
  simulation_table$prop_diff <- simulation_table$prop_diff * inversion_vector
  
  simulation_test_table <- as.data.frame(table(simulation_table$prop_diff))
  simulation_test_table$Var1 <- as.numeric(as.character(simulation_test_table$Var1))
  
  simulation_test_table <- simulation_test_table[!(simulation_test_table$Var1 %in% c(-1,1)),]
  
  simulation_test_table$probs <- simulation_test_table$Freq/sum(simulation_test_table$Freq)
  
  return(simulation_test_table)
}


#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

# Model functions


compute_base_probability_table <- function(max_read_depth,peak_probs_list){
  
  probability_table <- data.frame()
  
  for (read_depth in c(2:max_read_depth)){
    for (peak_prob_name in names(peak_probs_list)){
      data_r <-restricted_freq_distr_viable(read_depth, peak_probs_list[[peak_prob_name]])
      data_r$peak_type <- peak_prob_name
      #data_r$read_depth <- read_depth
      probability_table <- rbind(probability_table, data_r)
    }
  }
  return(probability_table)
}



restricted_freq_distr_viable <- function(read_depth, peak_probs){
  reads_a <- seq(1,read_depth-1,1)
  reads_b <- seq(read_depth-1,1,-1)
  
  # Strong future optimization of this double binom estimation 
  read_binom_probs_a <- dbinom(reads_a,size = read_depth, prob = peak_probs[1])
  read_binom_probs_b <- dbinom(reads_a,size = read_depth, prob = peak_probs[2])
  
  final_freqs <- (reads_a - reads_b) / read_depth
  
  final_freqs_probs <- (read_binom_probs_a + read_binom_probs_b ) / (sum(read_binom_probs_a,read_binom_probs_b))
  
  final_df <- data.frame(read_depth = read_depth, reads_a = reads_a , reads_b = reads_b,freqs = final_freqs, probs = final_freqs_probs)
  
  return(final_df)
  
}

compute_pmf_probs <- function(base_prob_table, read_depth_distr, central_prop){
  read_depth_table <- as.data.frame(table(read_depth_distr))
  read_depth_table$depth_prob <- read_depth_table$Freq / sum(read_depth_table$Freq) 
  
  out_ward_prop <- (1-central_prop)
  
  
  filtered_base_table <- base_prob_table[base_prob_table$read_depth %in% read_depth_table$read_depth_distr,]
  
  filtered_base_table[filtered_base_table$peak_type == "0.5_0.5",c("probs")] <- filtered_base_table[filtered_base_table$peak_type == "0.5_0.5",c("probs")] * central_prop
  
  filtered_base_table[filtered_base_table$peak_type == "0.25_0.75",c("probs")] <- filtered_base_table[filtered_base_table$peak_type == "0.25_0.75",c("probs")] * out_ward_prop
  
  
  filtered_base_table <- merge(filtered_base_table, read_depth_table, by.x = "read_depth", by.y = "read_depth_distr", all.x = TRUE)
  
  filtered_base_table$corrected_probs <- filtered_base_table$probs * filtered_base_table$depth_prob
  
  final_weights <- aggregate(corrected_probs ~ freqs, filtered_base_table, FUN = sum)
  
  # I rest 1 at the end because read depths of 1 get ignored.....
  #final_weights$probs <- final_weights$probs / (length(read_depth_table$read_depth_distr))
  
  return(final_weights)
}

simple_pmf_function <- function(x, final_pmf_weights){
  if (x %in% final_pmf_weights$freqs){
    return(final_pmf_weights[final_pmf_weights$freqs == x, c("probs")])
  }else{
    return(0)
  }
}

modified_pmf_function <- 
  
  simple_cdf_func <- function(x,weights) {
    # Check if weights is a data frame with columns 'freqs' and 'probs'
    if (!all(c("freqs", "probs") %in% names(weights))) {
      stop("weights must be a data frame with columns 'freqs' and 'probs'")
    }
    
    # Sort weights by 'freqs' in ascending order
    weights <- weights[order(weights$freqs), ]
    
    
    cdf_probs <-c()
    
    for (val in x){
      #print(val)
      # Initialize cumulative probability
      cdf <- 0
      # Loop through weights and calculate cumulative sum
      for (i in 1:nrow(weights)) {
        if (weights$freqs[i] <= val) {
          cdf <- cdf + weights$probs[i]
        } else {
          cdf_probs <- c(cdf_probs,cdf)
          break  # Stop iterating when freqs exceed x
        }
      }
      #print("Precheck")
      #print(weights$freqs[-1])
      if (val >= tail(weights$freqs,n=1)){
        cdf_probs <- c(cdf_probs,1)       
      }
    }
    return(cdf_probs)
  }


#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

