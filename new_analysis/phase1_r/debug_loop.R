# Debug script to test first 5 regions
suppressPackageStartupMessages({
  library(dlnm)
  library(data.table)
  library(arrow)
  library(splines)
})

cat("Loading data...\n")
mort <- as.data.table(read_parquet('../phase0_data_prep/results/mortality_immediate_daily_elderly.parquet'))
era5 <- as.data.table(read_parquet('../phase0_data_prep/results/era5_immediate_daily.parquet'))
setnames(mort, 'immediate_code', 'region_code')
setkey(mort, region_code, date)
setkey(era5, region_code, date)
data <- merge(mort, era5, by = c('region_code', 'date'))
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]
setorder(data, region_code, date)

MAX_LAG <- 21
TEMP_DF <- 4
LAG_DF <- 4
temp_boundary <- c(0.7, 35.1)

regions <- unique(data$region_code)[1:5]
cat("Testing first 5 regions...\n")

success_count <- 0

for (i in seq_along(regions)) {
  reg <- regions[i]
  daily <- data[region_code == reg][order(date)]
  n_days <- nrow(daily)
  
  cat("Region", i, ":", reg, "- n_days =", n_days, "\n")
  
  if (n_days < 365) {
    cat("  SKIP: insufficient data\n")
    next
  }
  
  result <- tryCatch({
    reg_temp_pcts <- quantile(daily$tmean, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
    cb <- crossbasis(daily$tmean, lag = MAX_LAG,
      argvar = list(fun = "ns", knots = reg_temp_pcts, Boundary.knots = temp_boundary),
      arglag = list(fun = "ns", df = LAG_DF))
    
    n_years <- length(unique(format(daily$date, "%Y")))
    model <- glm(deaths ~ cb + ns(as.numeric(date), df = 7 * n_years) + factor(format(date, "%u")),
      data = daily, family = quasipoisson(link = "log"), na.action = na.exclude)
    
    cb_idx <- grep("cb", names(coef(model)))
    cb_coef <- coef(model)[cb_idx]
    list(status = "success", n_coef = length(cb_coef))
  }, error = function(e) {
    list(status = "error", msg = conditionMessage(e))
  })
  
  cat("  Result status:", result$status, "\n")
  if (result$status == "success") {
    cat("  N coefficients:", result$n_coef, "\n")
    success_count <- success_count + 1
  } else if (result$status == "error") {
    cat("  Error:", result$msg, "\n")
  }
}

cat("\nTotal successful:", success_count, "\n")
