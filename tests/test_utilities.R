getTestFile <- function(filename, package="rqubic") {
  res <- system.file("extdata", filename,
                     package=package)
  if(file.exists(res))
    return(res)
  res <- sprintf("../inst/extdata/%s", filename)
  if(file.exists(res))
    return(res)
  stop(sprintf("%s cannot be found either in installed pack directory or in the source directory",
               filename))
}
