#####################
#Additional functions
#####################
roundup <- function(x, to = 10){
      to*(x%/%to + as.logical(x%%to))
}

DCount <- function(StartDate, EndDate){
      return(as.Date(EndDate, format = "%d-%m-%Y") - as.Date(StartDate, format = "%d-%m-%Y"))
}

DAddSub <- function(StartDate, NoOfDays){
      return(as.Date(StartDate, format = "%d-%m-%Y") - NoOfDays)
}