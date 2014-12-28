setup <- function(need_logging = TRUE) {
  setwd("~/Yandex.Disk/Курсовая 4 курс/BioProject")
  if (need_logging) {
    import_logging()
  } 
}

import_logging <- function() {
  if (!require("futile.logger")) {
    install.packages("futile.logger")
    library(futile.logger)
  }
}