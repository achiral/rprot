# functions.R

detach_all <- function() {
  basic.pkg <- c("package:stats", "package:graphics", "package:grDevices", 
                 "package:utils", "package:datasets", "package:methods", "package:base")
  pkg.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1 ,TRUE, FALSE)]
  pkg.list <- setdiff(pkg.list, basic.pkg)
  lapply(pkg.list, detach, character.only = TRUE)
}

# dir.choose <- function() {
#   system("osascript -e 'tell app \"RStudio\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
#          intern = FALSE, ignore.stderr = TRUE)
#   p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
#   return(ifelse(length(p), p, NA))
# }