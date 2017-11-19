library(readr)
library(lubridate)

san_juan_dengue <- read_csv("data-raw/San_Juan_Training_Data.csv")

save(san_juan_dengue, file = "data/san_juan_dengue.rdata")
