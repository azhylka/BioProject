setup <- function(need_logging) {
  setwd("~/Yandex.Disk/Курсовая 4 курс/")
  
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