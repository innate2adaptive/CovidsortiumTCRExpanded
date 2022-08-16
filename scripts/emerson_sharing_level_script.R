#file from Benny 25-07-2022
covid_tcrs <- read.csv("~/projects/benny_chain/csv/CDR3_up_invitro.csv", as.is = T)


adaptive_sharing <- readRDS( "~/projects/adaptive786/rds/sharing_aa.rds")
adaptive_sharing <- as.data.table(adaptive_sharing)
covid_tcrs1 <- left_join(covid_tcrs, adaptive_sharing, by  = c("CDR3" = "aminoAcid"))

#deal with missing sequences and missing lines
covid_tcrs1[is.na(covid_tcrs1$sharing_level),]$sharing_level <- 0
covid_tcrs1[covid_tcrs1 == "",]$sharing_level <- NA

#rename for consistency
covid_tcrs1 <- dplyr::rename(covid_tcrs1,  emerson_sharing_level = sharing_level)
write.csv(covid_tcrs1, "~/projects/benny_chain/csv/CDR3_up_invitro_with_sharing_level.csv")
covid_tcrs1 <- read.csv("~/projects/benny_chain/csv/CDR3_up_invitro_with_sharing_level.csv", as.is = T)

#remove large objects from the working environment
rm(adaptive_sharing)
rm(covid_tcrs)
rm(covid_tcrs1)