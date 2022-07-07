url <- "https://raw.githubusercontent.com/AMLab-Amsterdam/CEVAE/master/datasets/IHDP/csv/ihdp_npci_1.csv"

ihdp <- read.csv(url, header = FALSE)
colnames(ihdp) <- c("treatment", "y_factual", "y_cfactual", "mu0", "mu1",
                    paste0("X", seq_len(25)))

usethis::use_data(ihdp, overwrite = TRUE)
