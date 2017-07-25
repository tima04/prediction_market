library(data.table)
library(lubridate)

## isthere("ab", c("Aab", "daf")) == c(TRUE, FALSE)
isthere <- function(pat, xs) {
    match <- function(x)
        length(grep(pat, x, ignore.case = T, value = T)) > 0
    rslt <- unlist(Map(match, xs))
    names(rslt) <- NULL
    return (rslt)
}

## lowercase("Ab_CdE") ==  "ab_cde"
lowercase <- function(chars) {
    big <- "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    small <- "abcdefghijklmnopqrstuvwxyz"
    rslt <- ""
    for (char in strsplit(chars, "")[[1]]) {
            pos <- regexpr(char, big)
        if (pos %in% seq(1,26)) 
            rslt <- paste(rslt, (strsplit(small, "")[[1]][pos]), sep="")
        else 
            rslt <- paste(rslt, char, sep="")
    }
    rslt
}

to.seconds <- function(times) {
   as.numeric(dmy_hms(times)) 
}
