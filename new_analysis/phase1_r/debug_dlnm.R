# Debug script for DLNM analysis
suppressPackageStartupMessages({
  library(dlnm)
  library(mvmeta)
  library(data.table)
  library(arrow)
  library(splines)
})

# Load data
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

# Test on one region
reg <- unique(data$region_code)[1]
cat("Testing region:", reg, "\n")

daily <- data[region_code == reg][order(date)]
cat("  N days:", nrow(daily), "\n")
cat("  Deaths range:", min(daily$deaths), "-", max(daily$deaths), "\n")
cat("  Temp range:", round(min(daily$tmean), 1), "-", round(max(daily$tmean), 1), "\n")

# Parameters
MAX_LAG <- 21
TEMP_DF <- 4
LAG_DF <- 4
temp_boundary <- c(0.7, 35.1)

reg_temp_pcts <- quantile(daily$tmean, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
cat("  Temp knots:", paste(round(reg_temp_pcts, 1), collapse = ", "), "\n")

# Try creating crossbasis
cat("Creating crossbasis...\n")
cb <- tryCatch({
  crossbasis(
    daily$tmean,
    lag = MAX_LAG,
    argvar = list(fun = "ns", knots = reg_temp_pcts, Boundary.knots = temp_boundary),
    arglag = list(fun = "ns", df = LAG_DF)
  )
}, error = function(e) {
  cat("ERROR in crossbasis:", conditionMessage(e), "\n")
  return(NULL)
})

if (!is.null(cb)) {
  cat("  Crossbasis dim:", dim(cb), "\n")
  
  # Try fitting model
  cat("Fitting GLM...\n")
  daily[, date_num := as.numeric(date)]
  n_years <- length(unique(format(daily$date, "%Y")))
  cat("  N years:", n_years, "\n")
  
  model <- tryCatch({
    glm(
      deaths ~ cb + ns(date_num, df = 7 * n_years) + factor(format(date, "%u")),
      data = daily,
      family = quasipoisson(link = "log"),
      na.action = na.exclude
    )
  }, error = function(e) {
    cat("ERROR in GLM:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (!is.null(model)) {
    cat("  Model fitted!\n")
    cat("  N coefficients:", length(coef(model)), "\n")
    
    cb_idx <- grep("cb", names(coef(model)))
    cat("  Crossbasis indices:", length(cb_idx), "\n")
    
    cb_coef <- coef(model)[cb_idx]
    cat("  Coef range:", round(min(cb_coef), 4), "to", round(max(cb_coef), 4), "\n")
    
    # Try crosspred
    cat("Creating crosspred...\n")
    temp_seq <- seq(temp_boundary[1], temp_boundary[2], length.out = 100)
    
    cp <- tryCatch({
      crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(daily$tmean))
    }, error = function(e) {
      cat("ERROR in crosspred:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (!is.null(cp)) {
      cat("  Crosspred created!\n")
      mmt_idx <- which.min(cp$allRRfit)
      mmt <- temp_seq[mmt_idx]
      cat("  MMT:", round(mmt, 1), "\n")
      cat("  RR range:", round(min(cp$allRRfit), 3), "-", round(max(cp$allRRfit), 3), "\n")
    }
  }
}

cat("\nDone!\n")
