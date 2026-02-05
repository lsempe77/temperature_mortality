# Quick check of data availability
library(arrow)
library(data.table)

mort <- as.data.table(read_parquet("phase0_data_prep/results/mortality_intermediate_daily_elderly.parquet"))
weather <- as.data.table(read_parquet("phase0_data_prep/results/era5_intermediate_daily.parquet"))

cat("Mort date range:", as.character(range(mort$date)), "\n")
cat("Weather date range:", as.character(range(weather$date)), "\n")

cat("\nMort date class:", class(mort$date), "\n")
cat("Weather date class:", class(weather$date), "\n")

cat("\nMort date sample:", head(mort$date, 3), "\n")
cat("Weather date sample:", head(weather$date, 3), "\n")

# Filter to 2022-2024 and check
mort_post <- mort[year(date) >= 2022 & year(date) <= 2024]
weather_post <- weather[year(date) >= 2022 & year(date) <= 2024]

cat("\nMort 2022-2024 rows:", nrow(mort_post), "\n")
cat("Weather 2022-2024 rows:", nrow(weather_post), "\n")

# Try merge with explicit date conversion
mort$date <- as.Date(mort$date)
weather$date <- as.Date(weather$date)

data <- merge(mort, weather, by = c("region_code", "date"), all.x = TRUE)
data_post <- data[year(date) >= 2022 & year(date) <= 2024]
cat("\nAfter date conversion - NA in temp_mean:", sum(is.na(data_post$temp_mean)), "\n")
cat("Non-NA temps:", sum(!is.na(data_post$temp_mean)), "\n")
