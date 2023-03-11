## Table One, Stratified by AD Status and identified race and ehtnicity

## Install tableone and excel 
install.packages("tableone")
#install.packages("openxlsx") - NOT USED

## tableone package itself
library(tableone)
#library(openxlsx) - NOT USED

## Changed the working directory via GUI:
## Session -> Set Working Directory -> Choose Directory...
## NOTE: This is equivalent to the following command below:
setwd("set/to/working/directory")

## AD MatchIt
AD <- read.csv('Demographics/RE_MI_ad_demo.csv',
               header=TRUE)
## Control MatchIt
con <- read.csv('Demographics/RE_MI_con_demo.csv',
                header=TRUE) 


## Add ADStatus column for both
AD$ADStatus = 'AD'
con$ADStatus = 'control'

## consolidate dataframes
## https://stackoverflow.com/questions/8169323/r-concatenate-two-dataframes
all <- rbind(AD, con)

## Check all has the following dimensions: 5064 x 11
dim(all)

## Create sex column: only two sexes in dataset, no missing values (see above)
## Reference for which ID corresponds to which sex can be found in link below: 
## https://www.ohdsi.org/web/wiki/doku.php?id=documentation:vocabulary:gender
all$sex <- apply(all[,"gender_concept_id", drop=F], 1, function(x) {ifelse(test = x == 8532,
                                                                      yes = 'Female',
                                                                      no = 'Male')})

## Change death_status column: Alive (0) or Deceased (1)
## Reference for which ID corresponds to which sex can be found in link below: 
## https://www.ohdsi.org/web/wiki/doku.php?id=documentation:vocabulary:gender
all$death_status <- apply(all[,"death_status", drop=F], 1, function(x) {ifelse(test = x == 0,
                                                                           yes = 'Alive',
                                                                           no = 'Deceased')})
head(all)

## Check number of Females (n=3588) and Males (n=1476) is correct
sex_value_counts <- as.data.frame(table(all$sex))

sex_value_counts[order(-sex_value_counts$Freq),]

## Only keep demographic data for table one
## https://stackoverflow.com/questions/10085806/extracting-specific-columns-from-a-data-frame
all_demo <- all[c('estimated_age', 
                  'death_status', 
                  'UCSFDerivedRaceEthnicity_Clean', 
                  'sex', 
                  'ADStatus')]

## Change identified race and ethnicity names so that it's consistent
## with what will be submitted
## https://stackoverflow.com/questions/13871614/replacing-values-from-a-column-using-a-condition-in-r
all_demo$UCSFDerivedRaceEthnicity_Clean[all_demo$UCSFDerivedRaceEthnicity_Clean == 'Black or African American'] <- 'Black'
all_demo$UCSFDerivedRaceEthnicity_Clean[all_demo$UCSFDerivedRaceEthnicity_Clean == 'Latinx'] <- 'Latine'
all_demo$UCSFDerivedRaceEthnicity_Clean[all_demo$UCSFDerivedRaceEthnicity_Clean == 'White or Caucasian'] <- 'White'

## Check that identified race and ethnicity names have been modified (n=1266 each) 
re_value_counts <- as.data.frame(table(all_demo$UCSFDerivedRaceEthnicity_Clean))

re_value_counts[order(-re_value_counts$Freq),]

## Create Table One
## Main reference:
## http://rstudio-pubs-static.s3.amazonaws.com/13321_da314633db924dc78986a850813a50d5.html


vars <- c('estimated_age',
          'death_status', 
          'UCSFDerivedRaceEthnicity_Clean', 
          'sex', 
          'ADStatus')

strata <- c('ADStatus', 'UCSFDerivedRaceEthnicity_Clean')

## Create Table 1 stratified by trt (omit strata argument for overall table)
tableOne <- CreateTableOne(vars=vars, strata=strata, data=all_demo, smd=TRUE)
## Just typing the object name will invoke the print.TableOne method
## Tests are by oneway.test/t.test for continuous, chisq.test for categorical
tableOne


## To get SMD printed out:
## https://cran.r-project.org/web/packages/tableone/vignettes/smd.html
## And to show all levels (i.e., all levels for all variables, including binary):
## https://cran.r-project.org/web/packages/tableone/tableone.pdf
print(tableOne, smd=TRUE, showAllLevels=TRUE)

## Save
## https://charliemarks.com/r-tutorials/createtableone
tableOne_csv <- print(tableOne, smd=TRUE, showAllLevels=TRUE)

write.csv(tableOne_csv, "Tables/TableOne.csv")
