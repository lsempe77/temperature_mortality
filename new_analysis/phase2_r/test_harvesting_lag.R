# Test crosspred lag parameter behavior
library(dlnm)
library(data.table)
library(splines)
library(arrow)

cat("Loading data...\n")

# Use same loading logic as 02b_harvesting_v2.R
BASE_DIR <- "c:/Users/LucasSempe/OneDrive - International Initiative for Impact Evaluation/Desktop/sim_data/new_analysis"
DATA_DIR <- file.path(BASE_DIR, "phase0_data_prep", "results")
EXPOSURE_TYPE <- "intermediate"

mort <- as.data.table(read_parquet(file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_elderly.parquet"))))
era5 <- as.data.table(read_parquet(file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))))

setnames(mort, "intermediate_code", "region_code", skip_absent = TRUE)
setnames(era5, "intermediate_code", "region_code", skip_absent = TRUE)

mort[, date := as.character(as.Date(date))]
era5[, date := as.character(date)]

setkey(mort, region_code, date)
setkey(era5, region_code, date)
data <- merge(mort, era5, by = c("region_code", "date"))
data[, date := as.Date(date)]
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]
data <- data[!is.na(deaths) & !is.na(tmean)]
setDT(data)

# Take one region with lots of data
regions <- unique(data$region_code)
cat(sprintf("Regions: %d\n", length(regions)))
daily <- data[region_code == regions[1]][order(date)]
cat(sprintf("Region %s data: %d rows\n", regions[1], nrow(daily)))

# Fit model with lag=35
temp_pcts <- quantile(data$tmean, c(0.10, 0.75, 0.90), na.rm = TRUE)
temp_boundary <- c(0, 40)

cat("\nFitting model with lag=35...\n")
cb35 <- crossbasis(
  daily$tmean,
  lag = 35,
  argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
  arglag = list(fun = "ns", df = 4)
)

model35 <- glm(
  deaths ~ cb35 + ns(as.numeric(date), df = 7 * length(unique(format(daily$date, "%Y")))) +
    factor(format(date, "%u")),
  data = daily,
  family = quasipoisson(link = "log"),
  na.action = na.exclude
)

# Find MMT
mmt <- median(daily$tmean)
cat(sprintf("MMT: %.2f\n", mmt))

# Test different approaches for cumulative RR at different horizons
eval_temps <- c(11.5, 30)

cat("\n=== Full cumulative (default lag 0-35) ===\n")
cp_full <- crosspred(cb35, model35, at = eval_temps, cumul = TRUE, cen = mmt)
cat(sprintf("  Cold (11.5C): RR = %.4f\n", cp_full$allRRfit[1]))
cat(sprintf("  Heat (30C): RR = %.4f\n", cp_full$allRRfit[2]))

cat("\n=== Testing lag parameter with different specs ===\n")

# Test lag parameter
for (h in c(7, 14, 21, 28, 35)) {
  cat(sprintf("\nLag horizon %d:\n", h))
  
  # Try lag = c(0, h)
  tryCatch({
    cp <- crosspred(cb35, model35, at = eval_temps, cumul = TRUE, cen = mmt, lag = c(0, h))
    cat(sprintf("  lag=c(0,%d): Heat=%.4f Cold=%.4f\n", h, cp$allRRfit[2], cp$allRRfit[1]))
  }, error = function(e) {
    cat(sprintf("  Error with lag=c(0,%d): %s\n", h, e$message))
  })
  
  # Try just lag = h 
  tryCatch({
    cp <- crosspred(cb35, model35, at = eval_temps, cumul = TRUE, cen = mmt, lag = h)
    cat(sprintf("  lag=%d: Heat=%.4f Cold=%.4f\n", h, cp$allRRfit[2], cp$allRRfit[1]))
  }, error = function(e) {
    cat(sprintf("  Error with lag=%d: %s\n", h, e$message))
  })
}

cat("\n=== Method 2: Manual summation from matRRfit ===\n")
cat("(Sum coefficients directly from matRRfit which has lag-specific effects)\n")

cp_full_mat <- crosspred(cb35, model35, at = eval_temps, cumul = FALSE, cen = mmt)

# matRRfit has dimensions [n_temps x n_lags]
cat(sprintf("matRRfit dimensions: %d temps x %d lags\n", nrow(cp_full_mat$matRRfit), ncol(cp_full_mat$matRRfit)))

for (h in c(7, 14, 21, 28, 35)) {
  # Manual cumulative: sum log-RRs from lag 0 to lag h, then exponentiate
  # Actually for RR: multiply RRs or sum log-RRs
  heat_logRR_cumul <- sum(log(cp_full_mat$matRRfit[2, 1:(h+1)]))
  cold_logRR_cumul <- sum(log(cp_full_mat$matRRfit[1, 1:(h+1)]))
  
  cat(sprintf("  Lag 0-%d: Heat=%.4f Cold=%.4f\n", 
              h, exp(heat_logRR_cumul), exp(cold_logRR_cumul)))
}

cat("\nDone!\n")
