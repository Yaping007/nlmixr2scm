## Project startup: prefer the dev nlmixr2 stack in the scratch library.
## The scratch lib holds nlmixr2est 6.2.0 (analytic ODE gradient, ifoceif alias),
## which supersedes the 6.1.0 build in the default system library.
local({
  scratch <- file.path(path.expand("~"), "R-nlmixr2-dev")
  if (dir.exists(scratch)) {
    .libPaths(c(scratch, .libPaths()))
    ver <- tryCatch(
      as.character(utils::packageVersion("nlmixr2est", lib.loc = scratch)),
      error = function(e) "NOT INSTALLED"
    )
    message(sprintf("[.Rprofile] scratch lib on path; nlmixr2est %s", ver))
  } else {
    message("[.Rprofile] scratch lib not found: ", scratch)
  }
})
