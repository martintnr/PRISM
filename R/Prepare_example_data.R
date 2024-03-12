#' Title
#'
#' @return
#' @export
#'
#' @examples
Prepare_example_data <- function(){
  #setwd("/home/martin/Script/PleioVar")
  #load(file="data/simulated_example_data.rda")
  #load(file="data/VAR2.rda")
  #setwd("/home/martin/Script/PleioVar/test")

  if(!file.exists("Data/")){system("mkdir Data")}


  if(length(list.files("Data/")) != 0){
    print("Error: the Data folder must be empty")
    return()
  }
  # We set up the folder architecture


  ParametersTable <- simulated_example_data[[2]]
  write.table(ParametersTable, "Data/ParametersTable.csv", row.names = F, sep = ",")
  #We write the LHC-MR-like table parameters
  Traits <- unique(c(ParametersTable$X, ParametersTable$Y))


  createfiles <- function(A){

    Dat <- simulated_example_data[[1]][[A]]
    Dat$variant <- VAR2$variant[c(1:nrow(Dat))]
    Dat$rsid <- VAR2$rsid[c(1:nrow(Dat))]
    write.table(Dat, paste0("Data/", Traits[A],".csv"), row.names = F, sep = ",")

  }
  lapply(X = c(1:length(Traits)), FUN = createfiles)
  #We write all simulated GWAS data


  # We gzip files (if possible) to gain disk space
  tryCatch(
    { system("gzip Data/*.csv")
      message("Files were gzipped")
    },
    error = function(cond) {
      message("Files were not gzipped but that's okay")
    })


  print("Files created")


}
