# Dunn's test for UCSF Network Metrics

install.packages("hash") # version 2.2.6.2
install.packages("dunn.test") # version 1.3.5


# install libraries
library(dplyr)
library(hash)
library(dunn.test)

setwd(paste("\\\\working/directory/to/Tables/R_Dunn/",
            sep=""))

# UCSF
setwd(paste("\\\\working/directory/to/Tables/R_Dunn/",
            "UCSF/",
            sep=""))


# obtain network metrics
metrics <- dir()

# find out datatype of elements in the metrics variable
typeof(metrics)

# find out structure of the metrics variables
str(metrics)

# determine datatype of metrics variable (should be vector)
is.vector(metrics)

# %%%%%%%% AD Comparison %%%%%%%% 
cat(paste("**************",
          "AD Comparisons",
          "**************",
          sep="\n"))


cat("\n")

for (metric in metrics) {
  print(sprintf("Setting working directory for %s...", metric))
  setwd(paste("\\\\working/directory/to/Tables/R_Dunn/",
              "UCSF/",
              metric,
              sep=""))
  
  AD_metric_files <- dir(pattern='AD')
  
  all_AD_metric_values <- list() 
  
  for (AD_metric_file in AD_metric_files) {
    re_category <- substring(AD_metric_file, 1, 4)
    metric_name <- metric
    print(sprintf("Extracting %s for %s...", metric, re_category))
    AD_metric_values <- read.csv(AD_metric_file)
    AD_metric_values <- AD_metric_values %>% select(all_of(metric)) # from dplyr
    
    all_AD_metric_values[[re_category]] <- AD_metric_values
    
    cat("\n")
  }
  
  print(sprintf("Creating numeric vectors for %s...", metric))
  print("Checking order; should be A_AD, B_AD, L_AD, W_AD")
  cat(paste(names(all_AD_metric_values[1]),
            names(all_AD_metric_values[2]),
            names(all_AD_metric_values[3]),
            names(all_AD_metric_values[4]),
            sep="\n"))
  A_AD_vec <- as.numeric(unlist(all_AD_metric_values[1][1][1]))
  B_AD_vec <- as.numeric(unlist(all_AD_metric_values[2][1][1]))
  L_AD_vec <- as.numeric(unlist(all_AD_metric_values[3][1][1]))
  W_AD_vec <- as.numeric(unlist(all_AD_metric_values[4][1][1]))
  
  # Perform Dunn's test
  cat("\n")
  print(sprintf("Performing Dunn's test for %s...", metric))
  cat(paste("1 corresponds to Asian-identifed patients with AD",
            "2 corresponds to Black-identifed patients with AD",
            "3 corresponds to Latine-identifed patients with AD",
            "4 corresponds to White-identifed patients with AD",
            sep="\n"))
  dunn.test(list(A_AD_vec, B_AD_vec, L_AD_vec, W_AD_vec), method="bonferroni")
  
  cat("\n\n")
}

# %%%%%%%% Control Comparison %%%%%%%% 
cat(paste("*******************",
          "Control Comparisons",
          "*******************",
          sep="\n"))


cat("\n")

for (metric in metrics) {
  print(sprintf("Setting working directory for %s...", metric))
  setwd(paste("\\\working/directory/to/Tables/R_Dunn/",
              "UCSF/",
              metric,
              sep=""))
  
  con_metric_files <- dir(pattern='con')
  
  all_con_metric_values <- list() 
  
  for (con_metric_file in con_metric_files) {
    re_category <- substring(con_metric_file, 1, 5)
    metric_name <- metric
    print(sprintf("Extracting %s for %s...", metric, re_category))
    con_metric_values <- read.csv(con_metric_file)
    con_metric_values <- con_metric_values %>% select(all_of(metric)) # from dplyr
    
    all_con_metric_values[[re_category]] <- con_metric_values
    names(all_con_metric_values[1])
    
    cat("\n")
  }
  
  print(sprintf("Creating numeric vectors for %s...", metric))
  print("Checking order; should be A_con, B_con, L_con, W_con")
  cat(paste(names(all_con_metric_values[1]),
            names(all_con_metric_values[2]),
            names(all_con_metric_values[3]),
            names(all_con_metric_values[4]),
            sep="\n"))
  A_con_vec <- as.numeric(unlist(all_con_metric_values[1][1][1]))
  B_con_vec <- as.numeric(unlist(all_con_metric_values[2][1][1]))
  L_con_vec <- as.numeric(unlist(all_con_metric_values[3][1][1]))
  W_con_vec <- as.numeric(unlist(all_con_metric_values[4][1][1]))
  
  # Perform Dunn's test
  cat("\n")
  print(sprintf("Performing Dunn's test for %s...", metric))
  cat(paste("1 corresponds to Asian-identifed control patients",
            "2 corresponds to Black-identifed control patients",
            "3 corresponds to Latine-identifed control patients",
            "4 corresponds to White-identifed control patients",
            sep="\n"))
  dunn.test(list(A_con_vec, B_con_vec, L_con_vec, W_con_vec), method="bonferroni")
  
  cat("\n\n")
}
