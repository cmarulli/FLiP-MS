initialize_folders <- function(delete_existing_files = FALSE){
  
  dir.create("final", showWarnings = FALSE)
  dir.create("intermediate", showWarnings = FALSE)
  dir.create("plots", showWarnings = FALSE)
  
  if(delete_existing_files){
    unlink("final/*")
    unlink("intermediate/*")
    unlink("plots/*")
  }
  
  
  
}