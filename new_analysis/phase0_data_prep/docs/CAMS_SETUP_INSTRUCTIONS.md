# CAMS Air Pollution Data - Setup Instructions

## Current Status

**Issue Found:** Your `.cdsapirc` file is configured for the Climate Data Store (CDS), 
but CAMS pollution data is on the Atmosphere Data Store (ADS).

These are two separate systems with different accounts and API keys!

**Your current configuration:**
```
url: https://cds.climate.copernicus.eu/api   ← This is CDS (for ERA5)
key: 2a0468ca-...                             ← This is your CDS key
```

**What you need for CAMS:**
```
url: https://ads.atmosphere.copernicus.eu/api   ← This is ADS (for CAMS)
key: YOUR_ADS_TOKEN_HERE                        ← You need a NEW ADS key
```

---

## Quick Solution (10 minutes)

### Step 1: Create ADS Account

1. Go to: **https://ads.atmosphere.copernicus.eu/**
2. Click "Login/Register" (top right)
3. If you already have an ECMWF account, try logging in with same credentials
4. If not, create a new account

### Step 2: Get Your ADS API Key

1. Log into ADS
2. Click your username → "Your profile"
3. Find the "API Key" section
4. Copy the Personal Access Token

### Step 3: Update Your Configuration

Create a file at `C:\Users\LucasSempe\.cdsapirc_ads` with:

```
url: https://ads.atmosphere.copernicus.eu/api
key: YOUR_NEW_ADS_TOKEN
```

Then when running CAMS downloads, temporarily swap the files:
```powershell
# Before CAMS download:
Copy-Item "$HOME\.cdsapirc" "$HOME\.cdsapirc_cds"
Copy-Item "$HOME\.cdsapirc_ads" "$HOME\.cdsapirc"

# After CAMS download (restore CDS for ERA5):
Copy-Item "$HOME\.cdsapirc_cds" "$HOME\.cdsapirc"
```

### Step 4: Accept License Terms

**CRITICAL!** Before any download:
1. Go to: https://ads.atmosphere.copernicus.eu/datasets/cams-global-reanalysis-eac4
2. Click "Download" tab
3. Scroll to bottom
4. Accept the CC-BY licence

### Step 5: Run the Download Script

```powershell
cd "c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase0_data_prep"
python 00b_download_cams_pollution_v2.py
```

---

## Alternative: Manual Download

If API setup is problematic, download manually:

1. Go to: https://ads.atmosphere.copernicus.eu/datasets/cams-global-reanalysis-eac4
2. Click "Download" tab
3. Select:
   - Variable: `Particulate matter d < 2.5 µm`
   - Date: One month at a time (e.g., 2022-01-01 to 2022-01-31)
   - Time: 00:00, 06:00, 12:00, 18:00
   - Area: Sub-region → North: 5, West: -75, South: -35, East: -30
   - Format: Zipped netCDF
4. Submit and download

---

## Files Created

- `docs/CAMS_DATA_DOWNLOAD_GUIDE.md` - Complete documentation
- `00b_download_cams_pollution_v2.py` - Corrected download script

---

## Backup

Your original CDS config has been backed up to:
`C:\Users\LucasSempe\.cdsapirc_cds_backup`

This preserves your ERA5 access.

---

## Why Two Separate Systems?

| Feature | CDS (Climate) | ADS (Atmosphere) |
|---------|---------------|------------------|
| URL | cds.climate.copernicus.eu | ads.atmosphere.copernicus.eu |
| Data | ERA5, CMIP6, seasonal forecasts | CAMS pollution, composition |
| Account | Separate | Separate |
| API Key | Different | Different |

Both are part of Copernicus, but historically managed by different teams.
The good news: they're migrating to a unified system, but for now, you need both.

---

*Created: Phase 0 - Data Preparation*
