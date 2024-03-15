#' Prepare example data
#'
#' `Prepare_example_data()` creates a Data/ folder in the current working directory
#' and writes the dataset necessary to run the PleioVar example.
#'
#' @param gzip TRUE if you want the example data to be gzipped, FALSE otherwise
#'
#' @return
#' @export
#'
#' @examples
Prepare_example_data <- function(gzip = F){
  #setwd("/home/martin/Script/PleioVar")
  #load(file="data/simulated_example_data.rda")
  #load(file="data/Index.rda")
  #setwd("/home/martin/Script/PleioVar/test")

  if(!file.exists("Data/")){system("mkdir Data")}


  if(length(list.files("Data/")) != 0){
    print("Error: the Data folder must be empty")
    return()
  }
  # We set up the folder architecture


  ParametersTable <- simulated_example_data[[2]]
  #We write the LHC-MR-like table parameters
  Traits <- unique(c(ParametersTable$X, ParametersTable$Y))


  createfiles <- function(A){

    Dat <- simulated_example_data[[1]][[A]]

    write.table(Dat, paste0("Data/", Traits[A],".csv"), row.names = F, sep = ",")

  }
  lapply(X = c(1:length(Traits)), FUN = createfiles)
  #We write all simulated GWAS data


  # We gzip files (if possible) to gain disk space
    if(gzip == T){
      tryCatch(
        { system("gzip Data/*.csv")
          message("Summary statistics files were gzipped")
        },
        error = function(cond) {
          message("Summary statistics files were not gzipped but that's okay")
        })
    }

  write.table(ParametersTable, "Data/ParametersTable.csv", row.names = F, sep = ",")
  write.table(Index, "Data/Index.csv", row.names = F, sep = ",")


  gc()
  message("Files created")


}
