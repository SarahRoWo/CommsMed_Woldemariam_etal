# Dunn's test for UCSF UMAP RE

install.packages("hash") # version 2.2.6.2
install.packages("dunn.test") # version 1.3.5


# install libraries
library(dplyr)
library(hash)
library(dunn.test)
library(stringr) # for getting re_category

setwd(paste("\\\\working/directory/to/Tables/UMAP/",
            sep=""))


# obtain re umap dataframes
umap_dfs <- dir()

# find out datatype of elements in the umap_dfs variable
typeof(umap_dfs)

# find out structure of the umap_dfs variables
str(umap_dfs)

# determine datatype of umap_dfs variable (should be vector)
is.vector(umap_dfs)

# %%%%%%%% RE Comparison %%%%%%%% 
cat(paste("************************************",
          "RE Comparison UCSF - first component",
          "************************************",
          sep="\n"))


cat("\n")
  
all_RE_umap_values <- list() 

for (umap_df in umap_dfs) {
  re_category <- str_sub(umap_df, 14, -5)
  RE_umap_name <- re_category
  print(sprintf("Extracting values for %s-identified patients...", re_category))
  RE_metric_values <- read.csv(umap_df)
  
  # X0 is the name of the first umap component column
  RE_metric_values <- RE_metric_values %>% select(all_of("X0")) # from dplyr
  
  all_RE_umap_values[[re_category]] <- RE_metric_values
  
  cat("\n")
}

print(sprintf("Creating numeric vectors for component 1..."))
print("Checking order; should be Asian, Black, Latine, White")
cat(paste(names(all_RE_umap_values[1]),
          names(all_RE_umap_values[2]),
          names(all_RE_umap_values[3]),
          names(all_RE_umap_values[4]),
          sep="\n"))

# not sure why Asian category doesn't have column name
A_RE_vec <- as.numeric(unlist(all_RE_umap_values[1][1]))
B_RE_vec <- as.numeric(unlist(all_RE_umap_values[2][1]))
L_RE_vec <- as.numeric(unlist(all_RE_umap_values[3][1]))
W_RE_vec <- as.numeric(unlist(all_RE_umap_values[4][1]))

# Perform Dunn's test
cat("\n")
print(sprintf("Performing Dunn's test for first component:"))
cat(paste("1 corresponds to Asian-identifed patients",
          "2 corresponds to Black-identifed patients",
          "3 corresponds to Latine-identifed patients",
          "4 corresponds to White-identifed patients",
          sep="\n"))
dunn.test(list(A_RE_vec, B_RE_vec, L_RE_vec, W_RE_vec), method="bonferroni")

cat("\n\n")



cat(paste("*************************************",
          "RE Comparison UCSF - second component",
          "*************************************",
          sep="\n"))


cat("\n")

all_RE_umap_values <- list() 

for (umap_df in umap_dfs) {
  re_category <- str_sub(umap_df, 14, -5)
  RE_umap_name <- re_category
  print(sprintf("Extracting values for %s-identified patients...", re_category))
  RE_metric_values <- read.csv(umap_df)
  
  # X0 is the name of the first umap component column
  RE_metric_values <- RE_metric_values %>% select(all_of("X1")) # from dplyr
  
  all_RE_umap_values[[re_category]] <- RE_metric_values
  print(head(all_RE_umap_values[[re_category]]))
  
  cat("\n")
}

print(sprintf("Creating numeric vectors for component 1..."))
print("Checking order; should be Asian, Black, Latine, White")
cat(paste(names(all_RE_umap_values[1]),
          names(all_RE_umap_values[2]),
          names(all_RE_umap_values[3]),
          names(all_RE_umap_values[4]),
          sep="\n"))

# not sure why Asian category doesn't have column name
A_RE_vec <- as.numeric(unlist(all_RE_umap_values[1][1]))
B_RE_vec <- as.numeric(unlist(all_RE_umap_values[2][1]))
L_RE_vec <- as.numeric(unlist(all_RE_umap_values[3][1]))
W_RE_vec <- as.numeric(unlist(all_RE_umap_values[4][1]))

# Perform Dunn's test
cat("\n")
print(sprintf("Performing Dunn's test for second component"))
cat(paste("1 corresponds to Asian-identifed patients",
          "2 corresponds to Black-identifed patients",
          "3 corresponds to Latine-identifed patients",
          "4 corresponds to White-identifed patients",
          sep="\n"))
dunn.test(list(A_RE_vec, B_RE_vec, L_RE_vec, W_RE_vec), method="bonferroni")

cat("\n\n")