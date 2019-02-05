setwd("/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code")

start <- as.Date("2018-07-30")
today <- as.Date("2018-08-16")
all_days <- seq(start, today, by="day")

year <- as.POSIXlt(all_days)$year + 1900
urls <- paste0("http://cran-logs.rstudio.com/", year, "/", all_days, ".csv.gz")

# only download the files you don't have:
missing_days <- setdiff(as.character(all_days), 
                        tools::file_path_sans_ext(dir("CRANlogs"), TRUE))

for (i in 1:length(missing_days)) {
  print(paste0(i, "/", length(missing_days)))
  download.file(urls[i], paste0("CRANlogs/", missing_days[i], ".csv.gz"))
}


file_list <- list.files("CRANlogs", full.names=TRUE)

logs <- list()
for (file in file_list) {
  print(paste("Reading", file, "..."))
  clog <- read.table(file, header=TRUE, sep=",", quote='"', dec=".", 
                     fill=TRUE, comment.char="", as.is=TRUE)
  logs[[file]] <- clog[clog$package=="gren" & !is.na(clog$package), ]
}
rm(clog)

# rbind together all files
library(data.table)
dat <- rbindlist(logs)

# add some keys and define variable types
dat[, date:=as.Date(date)]
dat[, package:=factor(package)]
dat[, country:=factor(country)]
dat[, weekday:=weekdays(date)]
dat[, week:=strftime(as.POSIXlt(date),format="%Y-%W")]
setkey(dat, package, date, week, country)

# save data
save(dat, file="CRANlogs/CRANlogs.RData")


plot(unique(dat$date), sapply(unique(dat$date), function(d) {sum(dat$dat==d)}))


# for later analyses: load the saved data.table
load("CRANlogs/CRANlogs.RData")