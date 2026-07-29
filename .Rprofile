## Rtools45 ucrt64 runtime on PATH (needed for binary pkg loading, e.g. stringfish).
## R 4.5 is a ucrt build; its compiled-package runtime lives in rtools45/ucrt64/bin.
## Without this, binary packages depending on that runtime fail to load with
## "LoadLibrary failure: The specified module could not be found".
local({
  .rt <- c("C:/rtools45/ucrt64/bin", "C:/rtools45/usr/bin")
  .rt <- .rt[dir.exists(.rt)]
  if (length(.rt)) {
    .cur <- strsplit(Sys.getenv("PATH"), ";", fixed = TRUE)[[1]]
    .add <- .rt[!(.rt %in% .cur)]
    if (length(.add)) {
      Sys.setenv(PATH = paste(paste(.add, collapse = ";"), Sys.getenv("PATH"), sep = ";"))
    }
  }
})

## Report the active nlmixr2est build (now installed in the default system library;
## the former R-nlmixr2-dev scratch library has been retired).
local({
  ver <- tryCatch(
    as.character(utils::packageVersion("nlmixr2est")),
    error = function(e) "NOT INSTALLED"
  )
  message(sprintf("[.Rprofile] nlmixr2est %s", ver))
})
