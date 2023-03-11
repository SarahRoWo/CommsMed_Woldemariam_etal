
# install.packages("DBI")
library(DBI)
library(dplyr)
library(dbplyr)
library(odbc)
library(MatchIt)
library(ggplot2)
con <- dbConnect(odbc::odbc(),
                 Driver = "Driver",
                 Server = "Server",
                 Trusted_Connection = "yes")

library(odbc)

library(DBI)

# Set working directory to ensure it ends up in the directory's Demographics folder

# AD cohort (OMOP)
createTag = paste("SELECT * FROM personal_database")

alz_cohort <- dbGetQuery(con, createTag)
alz_cohort[,'isalz'] <- rep(1,nrow(alz_cohort))
sprintf("alz_cohort dimensions: %s %s", dim(alz_cohort)[1],dim(alz_cohort)[2])
head(alz_cohort)

createTag = paste("SELECT * FROM personal_database")

background_cohort <- dbGetQuery(con, createTag)
background_cohort[,'isalz'] <- rep(0, nrow(background_cohort))
sprintf("background_cohort dimensions: %s %s", dim(background_cohort)[1],dim(background_cohort)[2])

# put them in one giant dataframe. isalz column labels 1 as alzheimer cohort, and 0 as background cohort.
alz_background_pts <- rbind(alz_cohort,background_cohort)
alz_background_pts <- alz_background_pts[complete.cases(alz_background_pts),]
print(sum(is.na(alz_background_pts))) # Make sure no missing values
print(dim(alz_background_pts))

# run matchit. 
start_time <- Sys.time()

m.out = matchit(isalz ~ UCSFDerivedRaceEthnicity_Clean + estimated_age + gender_concept_id + death_status, 
                data = alz_background_pts, method = "nearest", ratio = 2)
m.data <- match.data(m.out,distance='pscore')

print('Saving data...')
save(m.data, m.out, file = paste0('AD_control_demographics.RData'))
write.csv(m.data,paste0('AD_control_demographics.csv'))

end_time <- Sys.time()
print(end_time - start_time)

summary(m.data)
summary(m.out)
