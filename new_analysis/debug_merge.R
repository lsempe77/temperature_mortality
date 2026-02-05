library(arrow)
library(data.table)

era5 <- as.data.table(read_parquet("c:/Users/LucasSempe/OneDrive - International Initiative for Impact Evaluation/Desktop/sim_data/new_analysis/phase0_data_prep/results/era5_intermediate_daily.parquet"))
mort <- as.data.table(read_parquet("c:/Users/LucasSempe/OneDrive - International Initiative for Impact Evaluation/Desktop/sim_data/new_analysis/phase0_data_prep/results/mortality_intermediate_daily_age_60_69.parquet"))

cat("ERA5:\n")
cat("  date class:", class(era5$date), "\n")
cat("  date sample:", head(as.character(era5$date), 3), "\n")
cat("  rows:", nrow(era5), "\n")

cat("\nMort:\n")
cat("  date class:", class(mort$date), "\n")
cat("  date sample:", head(as.character(mort$date), 3), "\n")
cat("  rows:", nrow(mort), "\n")

# Convert to character - IMPORTANT: mort has POSIXct with time, convert to Date first
era5[, date := as.character(date)]
mort[, date := as.character(as.Date(date))]

cat("\nAfter conversion:\n")
cat("  ERA5 date sample:", head(era5$date, 3), "\n")
cat("  Mort date sample:", head(mort$date, 3), "\n")

# Check overlap
cat("\nERA5 regions:", length(unique(era5$region_code)), "\n")
cat("Mort regions:", length(unique(mort$region_code)), "\n")
common_regions <- intersect(era5$region_code, mort$region_code)
cat("Common regions:", length(unique(common_regions)), "\n")

cat("\nERA5 date range:", min(era5$date), "-", max(era5$date), "\n")
cat("Mort date range:", min(mort$date), "-", max(mort$date), "\n")

# Try merge
setkey(mort, region_code, date)
setkey(era5, region_code, date)
merged <- merge(mort, era5, by = c("region_code", "date"))
cat("\nMerged rows:", nrow(merged), "\n")
