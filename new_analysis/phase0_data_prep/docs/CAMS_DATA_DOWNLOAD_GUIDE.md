# CAMS Global Reanalysis (EAC4) Data Download Guide

## Complete Step-by-Step Instructions for PM2.5 and Ozone Data

**Last Updated:** Phase 0 - Data Preparation  
**Purpose:** Air pollution data for confounding adjustment in heat-mortality analysis

---

## Table of Contents
1. [Overview](#1-overview)
2. [Account Registration](#2-account-registration)
3. [API Configuration](#3-api-configuration)
4. [Variable Reference](#4-variable-reference)
5. [Python Script Usage](#5-python-script-usage)
6. [Data Processing](#6-data-processing)
7. [Troubleshooting](#7-troubleshooting)
8. [Citation](#8-citation)

---

## 1. Overview

### Dataset Information
| Property | Value |
|----------|-------|
| **Dataset Name** | CAMS Global Reanalysis (EAC4) |
| **Dataset ID** | `cams-global-reanalysis-eac4` |
| **Spatial Resolution** | 0.75° × 0.75° (~83 km at equator) |
| **Temporal Resolution** | 3-hourly (analyses); hourly (surface forecasts) |
| **Temporal Coverage** | January 2003 – December 2024 |
| **Data Source** | Atmosphere Data Store (ADS) |
| **DOI** | [10.24381/d58bbf47](https://doi.org/10.24381/d58bbf47) |

### Key References
- Inness et al. (2019): "The CAMS reanalysis of atmospheric composition" - Atmos. Chem. Phys., 19, 3515–3556
- [Full documentation](https://confluence.ecmwf.int/x/OIX4B)

### Why CAMS EAC4?
- **Consistency:** Reanalysis data is consistent over time (unlike operational products)
- **Completeness:** No missing data (unlike ground-based monitoring)
- **Variables:** PM2.5, PM10, Ozone, NO2 and many others
- **Validation:** Extensively validated against surface observations

---

## 2. Account Registration

### Step 1: Create ADS Account

1. Go to: **https://ads.atmosphere.copernicus.eu/**
2. Click "Login/Register" (top right)
3. Select "Register" if you don't have an account
4. Fill in:
   - Email address
   - First name, Last name
   - Organization
   - Purpose of use
5. Verify your email
6. Complete profile setup

### Step 2: Accept License Terms

**CRITICAL:** You must accept the license for each dataset before downloading!

1. Go to dataset page: https://ads.atmosphere.copernicus.eu/datasets/cams-global-reanalysis-eac4
2. Click "Download" tab
3. Scroll down to "Terms of use"
4. Check the box to accept the CC-BY licence
5. Click "Accept terms"

---

## 3. API Configuration

### Step 1: Get Your Personal Access Token

1. Log into ADS: https://ads.atmosphere.copernicus.eu/
2. Click your username (top right)
3. Click "Your profile" or go to: https://ads.atmosphere.copernicus.eu/profile
4. Find your **Personal Access Token** in the "API Key" section
5. Copy the token (long alphanumeric string)

### Step 2: Create .cdsapirc File

**IMPORTANT:** The ADS API uses a DIFFERENT URL than the CDS API!

#### Windows Location
Create file: `C:\Users\<YourUsername>\.cdsapirc`

#### File Contents
```
url: https://ads.atmosphere.copernicus.eu/api
key: YOUR_PERSONAL_ACCESS_TOKEN_HERE
```

**Example:**
```
url: https://ads.atmosphere.copernicus.eu/api
key: 12345678-abcd-1234-efgh-123456789012
```

#### Creating the File on Windows

**Option A: Using PowerShell**
```powershell
# Navigate to home directory
cd $HOME

# Create the file
@"
url: https://ads.atmosphere.copernicus.eu/api
key: YOUR_TOKEN_HERE
"@ | Out-File -FilePath ".cdsapirc" -Encoding utf8
```

**Option B: Using Notepad**
1. Open Notepad
2. Paste the url and key content
3. Save as: `C:\Users\YourUsername\.cdsapirc`
4. **Important:** Save as "All Files" type, not .txt

### Step 3: Install Python Packages

```bash
pip install "cdsapi>=0.7.7"
pip install xarray netcdf4 numpy pandas
```

### Step 4: Verify Configuration

```python
import cdsapi
client = cdsapi.Client()
print("API URL:", client.url)
print("Connection successful!")
```

---

## 4. Variable Reference

### Variables for Heat-Mortality Analysis

| Variable | API Name | Units | Description | Access |
|----------|----------|-------|-------------|--------|
| **PM2.5** | `particulate_matter_2.5um` | kg m⁻³ | Fine particulate matter | Fast |
| **PM10** | `particulate_matter_10um` | kg m⁻³ | Coarse particulate matter | Fast |
| **PM1** | `particulate_matter_1um` | kg m⁻³ | Ultrafine particulate matter | Fast |
| **O3 (column)** | `total_column_ozone` | kg m⁻² | Total column ozone | Fast |

### Unit Conversions

**PM2.5 (kg m⁻³ → µg m⁻³):**
```python
pm25_ugm3 = pm25_kgm3 * 1e9  # Multiply by 10^9
```

**Surface Ozone:** For surface-level ozone, use multi-level `ozone` variable at model level 60 (surface). Units are kg kg⁻¹ (mass mixing ratio).

### Brazil Bounding Box

```python
# Brazil coverage with buffer
area = [5, -75, -35, -30]  # [North, West, South, East]

# Explanation:
# North:  5°N  (covers Roraima)
# West:  75°W  (covers Acre, Amazon)  
# South: 35°S  (covers Rio Grande do Sul)
# East:  30°W  (covers Fernando de Noronha)
```

---

## 5. Python Script Usage

### Basic API Request Structure

```python
import cdsapi

client = cdsapi.Client()

# CAMS EAC4 request for PM2.5
request = {
    'variable': ['particulate_matter_2.5um'],
    'date': '2022-01-01/2022-12-31',
    'time': ['00:00', '06:00', '12:00', '18:00'],
    'area': [5, -75, -35, -30],  # [N, W, S, E]
    'data_format': 'netcdf_zip',
}

client.retrieve(
    'cams-global-reanalysis-eac4',
    request,
    'cams_pm25_2022.zip'
)
```

### Full Download Script

See: `00b_download_cams_pollution.py`

This script:
- Downloads PM2.5, PM10, and total column ozone
- Uses quarterly downloads to avoid timeout
- Processes to daily means
- Saves as NetCDF and Parquet

### Data Request Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `variable` | List of variable names | `['particulate_matter_2.5um']` |
| `date` | Date range (YYYY-MM-DD/YYYY-MM-DD) | `'2022-01-01/2022-03-31'` |
| `time` | Analysis times (3-hourly) | `['00:00', '03:00', '06:00', ...]` |
| `area` | Bounding box [N, W, S, E] | `[5, -75, -35, -30]` |
| `data_format` | Output format | `'netcdf_zip'` or `'grib'` |

### Available Times (3-hourly analyses)

```python
times = ['00:00', '03:00', '06:00', '09:00', 
         '12:00', '15:00', '18:00', '21:00']
```

---

## 6. Data Processing

### Loading Downloaded Data

```python
import xarray as xr
import zipfile
import os

# Extract zip file
with zipfile.ZipFile('cams_pm25_2022.zip', 'r') as z:
    z.extractall('cams_data/')

# Load NetCDF
ds = xr.open_dataset('cams_data/data.nc')

# Check variables
print(ds.data_vars)
# Expected: pm2p5 (or similar short name)
```

### Computing Daily Means

```python
# Convert from kg/m³ to µg/m³
ds['pm25_ugm3'] = ds['pm2p5'] * 1e9

# Compute daily mean
daily = ds.resample(time='1D').mean()

# Save
daily.to_netcdf('cams_pm25_daily_2022.nc')
```

### Spatial Aggregation to Brazilian States

```python
import geopandas as gpd
import regionmask
import numpy as np

# Load Brazil shapefile
brazil = gpd.read_file('brazil_states.gpkg')

# Create region mask
mask = regionmask.mask_geopandas(brazil, ds['longitude'], ds['latitude'])

# Aggregate by state (population-weighted or area-weighted)
for state_id in brazil.index:
    state_mask = (mask == state_id)
    ds_state = ds.where(state_mask, drop=True)
    state_mean = float(ds_state.mean())
```

---

## 7. Troubleshooting

### Error: 404 Not Found

**Problem:** API returns 404 error for dataset.

**Solution:** Check your `.cdsapirc` file:
```
# WRONG - This is the CDS URL (Climate Data Store)
url: https://cds.climate.copernicus.eu/api/v2

# CORRECT - This is the ADS URL (Atmosphere Data Store)
url: https://ads.atmosphere.copernicus.eu/api
```

### Error: "Please agree to terms"

**Problem:** License terms not accepted.

**Solution:**
1. Go to https://ads.atmosphere.copernicus.eu/datasets/cams-global-reanalysis-eac4
2. Click "Download" tab
3. Accept the CC-BY licence at the bottom

### Error: Connection/Authentication Failed

**Solutions:**
1. Check token is correct (no spaces, complete string)
2. Verify file is named exactly `.cdsapirc` (not `.cdsapirc.txt`)
3. Check file location: `C:\Users\<YourUsername>\.cdsapirc`
4. Try regenerating your access token in profile

### Error: Request Too Large

**Problem:** Request times out for large date ranges.

**Solution:** Split into smaller chunks:
```python
# Instead of full year, download by quarter
quarters = [
    ('2022-01-01', '2022-03-31'),
    ('2022-04-01', '2022-06-30'),
    ('2022-07-01', '2022-09-30'),
    ('2022-10-01', '2022-12-31'),
]
```

### Slow Downloads

**Note:** CAMS data downloads can take 10-30 minutes per quarter depending on:
- Server load
- Request size
- Network speed

Consider:
- Downloading during off-peak hours (European nighttime)
- Using background downloads
- Breaking into smaller spatial regions

---

## 8. Citation

### Recommended Citation

```
Inness, A., Ades, M., Agustí-Panareda, A., et al. (2019). 
The CAMS reanalysis of atmospheric composition. 
Atmospheric Chemistry and Physics, 19, 3515–3556.
https://doi.org/10.5194/acp-19-3515-2019
```

### Data Citation

```
CAMS Reanalysis data (2003-2024). 
Copernicus Atmosphere Monitoring Service (CAMS).
https://doi.org/10.24381/d58bbf47
```

### Acknowledgement Text

```
This study used CAMS Reanalysis data from the Copernicus Atmosphere 
Monitoring Service (CAMS), implemented by ECMWF. CAMS data is published 
under a CC-BY licence. Neither ECMWF nor Copernicus endorse any work 
created using this data.
```

---

## Quick Reference Card

```
┌────────────────────────────────────────────────────────────────────┐
│ CAMS EAC4 Quick Reference                                         │
├────────────────────────────────────────────────────────────────────┤
│ API URL:    https://ads.atmosphere.copernicus.eu/api              │
│ Dataset:    cams-global-reanalysis-eac4                           │
│ Config:     C:\Users\<USER>\.cdsapirc                             │
├────────────────────────────────────────────────────────────────────┤
│ PM2.5:      particulate_matter_2.5um  (kg/m³ × 1e9 = µg/m³)      │
│ PM10:       particulate_matter_10um   (kg/m³ × 1e9 = µg/m³)      │
│ O3 column:  total_column_ozone        (kg/m²)                     │
├────────────────────────────────────────────────────────────────────┤
│ Brazil:     area = [5, -75, -35, -30]  # [N, W, S, E]            │
│ Times:      00:00, 03:00, 06:00, ... (every 3 hours)             │
│ Coverage:   2003-01-01 to 2024-12-31                              │
└────────────────────────────────────────────────────────────────────┘
```

---

*Document created for Heat-Mortality Brazil Analysis - Phase 0 Data Preparation*
