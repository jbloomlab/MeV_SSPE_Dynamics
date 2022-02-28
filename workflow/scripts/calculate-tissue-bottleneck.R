## ---------------------------
##
## Purpose of script: This script calculates the transmission bottleneck from one tissue to all others.
## Estimation techniques from Sobel Leonard et. al., 2018
##
## Author: Will Hannon
##
## Date Created: 2022-01-31
##
## Copyright (c) Will Hannon, 2022
##
## Email: wwh22@uw.edu
##
## ---------------------------

library(tidyverse)
library(rmutil)

## ---------------------------
## Function to calculate the transmission bottleneck

estimate.transmission.bottleneck = function(donor_recip_df, var_calling_threshold, Nb_min, Nb_max, confidence_interval, method){
  
  # table of SNP frequencies in donor and recipient
  donor_and_recip_freqs_observed = donor_recip_df %>% 
    # filter out variants below thershold in donor
    filter(donor_freqs >= var_calling_threshold)
  
  # the number of variants
  n_variants <- nrow(donor_and_recip_freqs_observed)
  
  if (method == "approx") {
    
    # Approxamate method
    # ----- 
    # Now implement the beta binomial algorithm functions
    
    # this function gives Log Likelihood for every donor recipient SNP frequency pair
    Log_Beta_Binom = function(nu_donor, nu_recipient, NB_SIZE){ 
      # used for recipient frequencies above calling threshold
      LL_val_above = 0
      # used for recipient frequencies below calling threshold
      LL_val_below = 0 
      
      for(k in 0:NB_SIZE){
        LL_val_above = LL_val_above + dbeta(nu_recipient, k, (NB_SIZE - k) ) * dbinom(k, size = NB_SIZE, prob = nu_donor)  
        LL_val_below = LL_val_below + pbeta(var_calling_threshold, k, (NB_SIZE - k)) * dbinom(k,size = NB_SIZE, prob = nu_donor)
      }
      
      # we use LL_val_above above the calling threshold, and LL_val_below below the calling threshold
      LL_val = if_else(nu_recipient >= var_calling_threshold, LL_val_above, LL_val_below)
      # convert likelihood to log likelihood
      LL_val <- log(LL_val) 
      
      return(LL_val)
    }
    
    # this function sums over all SNP frequencies in the donor and recipient
    LL_func_approx = function(Nb_size){ 
      Total_LL = 0
      LL_array = Log_Beta_Binom(donor_and_recip_freqs_observed$donor_freqs, donor_and_recip_freqs_observed$recip_freqs, Nb_size)  
      Total_LL = sum(LL_array)
      return(Total_LL)
    }
    
    # ----- 
    # Calculate the bottlenecks 
    
    # now we define an empty array of Log Likelihoods for all possible bottleneck sizes
    LL_tibble = tibble(bottleneck_size = c(Nb_min:Nb_max), Log_Likelihood = 0*c(Nb_min:Nb_max)) 
    
    # apply LL_func_approx to each row
    for(i in 1:nrow(LL_tibble)){
      LL_tibble$Log_Likelihood[i] = LL_func_approx(LL_tibble$bottleneck_size[i])
    }
    
    # ----- 
    # Find the maximum likelihood estimate and the associated confidence interval
    
    # maximum value of log likelihood
    Max_LL = max(LL_tibble$Log_Likelihood) 
    
    # bottleneck size at which max likelihood occurs and that bottleneck
    Max_LL_bottleneck_index = which(LL_tibble$Log_Likelihood == max(LL_tibble$Log_Likelihood)) 
    Max_LL_bottleneck = Max_LL_bottleneck_index + Nb_min - 1
    
    # necessary ratio of likelihoods set by confidence level
    likelihood_ratio = qchisq(confidence_interval, df=1) 
    # find the condifence intervals
    ci_tibble = filter(LL_tibble, 2*(Max_LL - Log_Likelihood) <= likelihood_ratio) 
    lower_CI_bottleneck = min(ci_tibble$bottleneck_size) # lower bound of confidence interval
    upper_CI_bottleneck = max(ci_tibble$bottleneck_size) # upper bound of confidence interval
    
    # ----- 
    # Check how the analysis when and print messages
    
    # if the condfidence interval table is empty
    if (length(ci_tibble$Log_Likelihood) == 0){
      lower_CI_bottleneck = min(Max_LL_bottleneck) 
      upper_CI_bottleneck = max(Max_LL_bottleneck)
    }
    # if the bottleneck is likely larger then the max set by user
    if(max(Max_LL_bottleneck) == Nb_max){
      upper_CI_bottleneck <- Nb_max
      print("Peak bottleneck value for MLE is at Nb_max! Try raising Nb_max for better bottleneck estimate.")
    }
    # if the bottleneck is likely smaller then the min set by user
    if(min(Max_LL_bottleneck) == Nb_min){
      lower_CI_bottleneck <- Nb_min
      if(Nb_min > 1){
        print("Peak bottleneck value for MLE is at Nb_min! Try lowering Nb_min for better bottleneck estimate.")}
    }
    
    # ----- 
    
  } else if (method == "exact") {
    
    # Exact method
    # ----- 
    # Now implement the beta binomial algorithm functions
    
    # this function gives Log Likelihood for every donor recipient SNP frequency pair
    Log_Beta_Binom = function(nu_donor, recip_total_reads, recip_var_reads, NB_SIZE){ 
      # used for recipient reads above calling threshold
      LL_val_above = 0 
      # used for recipient reads below calling threshold
      LL_val_below = 0 
      
      for(k in 0:NB_SIZE){
        alpha = k + 10^-9
        beta = (NB_SIZE - k) + 10^-9
        m = alpha/(alpha + beta)
        s = (alpha + beta)
        LL_val_above = LL_val_above + dbetabinom(recip_var_reads, recip_total_reads, m, s) * dbinom(k, size = NB_SIZE, prob = nu_donor)
        LL_val_below = LL_val_below + pbetabinom(floor(var_calling_threshold*recip_total_reads), recip_total_reads, m, s) * dbinom(k, size = NB_SIZE, prob = nu_donor)
      }
      
      # we use LL_val_above above the calling threshold, and LL_val_below below the calling threshold
      LL_val = if_else(recip_var_reads >= var_calling_threshold*recip_total_reads, LL_val_above, LL_val_below )
      
      # convert likelihood to log likelihood
      LL_val = log(LL_val) 
      
      return(LL_val)
    }
    
    # this function sums over all SNP frequencies in the donor and recipient
    LL_func_approx = function(Nb_size){  
      Total_LL = 0
      LL_array = Log_Beta_Binom(donor_and_recip_freqs_observed$donor_freqs, donor_and_recip_freqs_observed$recip_total_reads, donor_and_recip_freqs_observed$recip_var_reads, Nb_size)  
      Total_LL = sum(LL_array)
      return(Total_LL)
    }
    
    # ----- 
    # Calculate the bottlenecks 
    
    # now we define array of Log Likelihoods for all possible bottleneck sizes
    LL_tibble = tibble(bottleneck_size = c(Nb_min:Nb_max), Log_Likelihood = 0 * c(Nb_min:Nb_max)) 
    
    # apply LL_func_approx to each row
    for(i in 1:nrow(LL_tibble)){
      LL_tibble$Log_Likelihood[i] = LL_func_approx(LL_tibble$bottleneck_size[i])
    }
    
    # ----- 
    # Find the maximum likelihood estimate and the associated confidence interval
    
    # maximum value of log likelihood
    Max_LL = max(LL_tibble$Log_Likelihood) 
    
    # bottleneck size at which max likelihood occurs
    Max_LL_bottleneck_index = which(LL_tibble$Log_Likelihood == max(LL_tibble$Log_Likelihood))
    Max_LL_bottleneck = Max_LL_bottleneck_index + Nb_min - 1
    
    # necessary ratio of likelihoods set by confidence level
    likelihood_ratio = qchisq(confidence_interval, df=1) 
    # find the condifence intervals
    ci_tibble = filter(LL_tibble, 2*(Max_LL - Log_Likelihood) <= likelihood_ratio)
    lower_CI_bottleneck = min(ci_tibble$bottleneck_size) # lower bound of confidence interval
    upper_CI_bottleneck = max(ci_tibble$bottleneck_size) # upper bound of confidence interval
    
    # ----- 
    # Check how the analysis when and print messages
    
    # if the condfidence interval table is empty
    if (length(ci_tibble$Log_Likelihood) == 0) {
      lower_CI_bottleneck = min(Max_LL_bottleneck) 
      upper_CI_bottleneck = max(Max_LL_bottleneck)
    }
    # if the bottleneck is likely larger then the max set by user
    if(max(Max_LL_bottleneck) == Nb_max) {
      upper_CI_bottleneck = Nb_max
      print("Peak bottleneck value for MLE is at Nb_max! Try raising Nb_max for better bottleneck estimate.")
    }
    # if the bottleneck is likely smaller then the min set by user
    if(min(Max_LL_bottleneck) == Nb_min) {
      lower_CI_bottleneck = Nb_min
      if(Nb_min > 1) {
        print("Peak bottleneck value for MLE is at Nb_min! Try lowering Nb_min for better bottleneck estimate.")
      }
    }
    # ----- 
    
  } else {print("You need to choose the correct method, either 'approx' or 'exact'.")}
  
  print("Bottleneck size")
  if(length(Max_LL_bottleneck) > 1){print("MLE is degenerate. The best bottleneck values are")}
  print(Max_LL_bottleneck)
  print("confidence interval left bound")
  print(lower_CI_bottleneck)
  print("confidence interval right bound")
  print(upper_CI_bottleneck)
  
  return(list(Max_LL_bottleneck, lower_CI_bottleneck, upper_CI_bottleneck))
  
}


## ---------------------------
## Snakemake inputs
variant.path = snakemake@input[[1]] 
output.path = snakemake@output[[1]] 
tissue.name = snakemake@wildcards[["accession"]]
max.bottleneck = snakemake@params[["max_bottleneck"]]

# Read in the variant data
variant.df = read_csv(variant.path)

# Get the tissues 
tissue.target = variant.df %>% filter(Accession == tissue.name) %>% pull(Tissue) %>% unique()

# Get the remaining tissues
tissues = unique(variant.df$Tissue)[!unique(variant.df$Tissue) %in% tissue.target]

tv.df = data.frame()
results.df = data.frame()

# Calculate the bottleneck against all other tissues.
for (tissue in tissues) {

  donor.df = variant.df %>% 
    filter(Tissue == tissue.target) %>% 
    select(SNP, AF, DP)
  
  contact.df = variant.df %>% 
    filter(Tissue == tissue) %>% 
    select(SNP, AF, DP)
  
  bottleneck.df = full_join(donor.df, contact.df,
                            by = 'SNP',
                            suffix = c(".donor", ".contact")) %>% 
    mutate_each(list( ~ replace(., which(is.na(.)), 0))) %>% 
    mutate(Donor = tissue.target,
           Contact = tissue)
  
  tv.df = rbind(tv.df, bottleneck.df)
  
  # -- Calculate the Bottleneck -- #  
  bottleneck.calc.df = bottleneck.df %>% 
    select(donor_freqs = AF.donor,
           recip_freqs = AF.contact,
           recip_total_reads = DP.contact) %>% 
    mutate(recip_var_reads = recip_total_reads*recip_freqs) %>% 
    mutate(recip_var_reads = round(recip_var_reads))
  
  
  print(paste("Calculating bottleneck for", unique(bottleneck.df$Donor), "vs.", unique(bottleneck.df$Contact)))
  
  
  results = estimate.transmission.bottleneck(bottleneck.calc.df,
                                             var_calling_threshold = 0.02,
                                             Nb_min = 1,
                                             Nb_max = max.bottleneck,
                                             confidence_interval = 0.95,
                                             method = "exact")
  
  row = data.frame("donor" = tissue.target,
                   "recipient" = tissue,
                   "lower_ci" = results[[2]],
                   "upper_ci" = results[[3]],
                   "bottleneck" = results[[1]])
  
  results.df = rbind(results.df, row)
      
}

# Write out the results
write_csv(results.df, output.path)
