"""
Download ERA5 temperature and related variables for Brazil (2010-2020)
For studying the relationship between temperature and mortality in older people

Variables downloaded:
- 2m_temperature: Main exposure variable (hourly mean temperature)
- 2m_dewpoint_temperature: For calculating humidity/heat index
- maximum_2m_temperature_since_previous_post_processing: Daily maximum temperature
- minimum_2m_temperature_since_previous_post_processing: Daily minimum temperature
- total_precipitation: Weather patterns
- 10m_u_component_of_wind: Wind component (for wind chill/heat dissipation)
- 10m_v_component_of_wind: Wind component
- surface_pressure: Atmospheric pressure
- surface_solar_radiation_downwards: Solar radiation

Brazil bounding box: North: 5.3, South: -33.8, West: -73.9, East: -34.8
"""

import cdsapi


def download_era5_brazil():
    """Download ERA5 data for Brazil from 2010 to 2020."""
    
    client = cdsapi.Client()
    
    # Brazil bounding box (North, West, South, East)
    area = [6, -74, -34, -34]
    
    # Years to download
    # Modified to 2010-2020 since 2021-2024 are already downloaded
    years = [str(year) for year in range(2010, 2021)]
    
    # All months
    months = [f'{m:02d}' for m in range(1, 13)]
    
    # All days
    days = [f'{d:02d}' for d in range(1, 32)]
    
    # Selected hours (every 3 hours to reduce file size, adjust as needed)
    hours = ['00:00', '03:00', '06:00', '09:00', '12:00', '15:00', '18:00', '21:00']
    
    # =========================================================================
    # Download 1: Instantaneous variables (temperature, dewpoint, wind, pressure)
    # These are the core variables for temperature-mortality studies
    # =========================================================================
    
    instantaneous_variables = [
        '2m_temperature',                 # Main exposure variable
        '2m_dewpoint_temperature',        # For humidity/heat index calculation
        '10m_u_component_of_wind',        # Wind - eastward component
        '10m_v_component_of_wind',        # Wind - northward component
        'surface_pressure',               # Atmospheric pressure
    ]
    
    print("=" * 60)
    print("Downloading ERA5 data for Brazil (2010-2020)")
    print("For temperature-mortality research in older populations")
    print("=" * 60)
    
    print("\n[1/3] Downloading instantaneous variables...")
    print(f"      Variables: {', '.join(instantaneous_variables)}")
    
    for year in years:
        output_file = f'era5_brazil_hourly_{year}.nc'
        print(f"      Downloading {year}...")
        
        client.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': instantaneous_variables,
                'year': year,
                'month': months,
                'day': days,
                'time': hours,
                'area': area,
                'data_format': 'netcdf',
            },
            output_file
        )
        print(f"      Saved: {output_file}")
    
    # =========================================================================
    # Download 2: Min/Max temperature (critical for mortality studies)
    # Daily extremes are key predictors of heat/cold-related mortality
    # =========================================================================
    
    minmax_variables = [
        'maximum_2m_temperature_since_previous_post_processing',
        'minimum_2m_temperature_since_previous_post_processing',
    ]
    
    print("\n[2/3] Downloading min/max temperature variables...")
    print(f"      Variables: {', '.join(minmax_variables)}")
    
    for year in years:
        output_file = f'era5_brazil_minmax_{year}.nc'
        print(f"      Downloading {year}...")
        
        client.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': minmax_variables,
                'year': year,
                'month': months,
                'day': days,
                'time': hours,
                'area': area,
                'data_format': 'netcdf',
            },
            output_file
        )
        print(f"      Saved: {output_file}")
    
    # =========================================================================
    # Download 3: Accumulated variables (precipitation, radiation)
    # Additional environmental factors that may confound or modify effects
    # =========================================================================
    
    accumulated_variables = [
        'total_precipitation',                 # Rainfall
        'surface_solar_radiation_downwards',   # Solar radiation
    ]
    
    print("\n[3/3] Downloading accumulated variables...")
    print(f"      Variables: {', '.join(accumulated_variables)}")
    
    for year in years:
        output_file = f'era5_brazil_accumulated_{year}.nc'
        print(f"      Downloading {year}...")
        
        client.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': accumulated_variables,
                'year': year,
                'month': months,
                'day': days,
                'time': hours,
                'area': area,
                'data_format': 'netcdf',
            },
            output_file
        )
        print(f"      Saved: {output_file}")
    
    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 60)
    print("DOWNLOAD COMPLETE!")
    print("=" * 60)
    print("\nFiles created for each year (2010-2020):")
    print("  - era5_brazil_hourly_YYYY.nc")
    print("      Contains: 2m temperature, dewpoint, wind components, pressure")
    print("  - era5_brazil_minmax_YYYY.nc")
    print("      Contains: Daily max/min 2m temperature")
    print("  - era5_brazil_accumulated_YYYY.nc")
    print("      Contains: Total precipitation, solar radiation")
    
    print("\n" + "-" * 60)
    print("NOTES FOR ANALYSIS:")
    print("-" * 60)
    print("• Temperature is in Kelvin. Convert to Celsius: T(°C) = T(K) - 273.15")
    print("• Use 2m_temperature + 2m_dewpoint_temperature to calculate:")
    print("  - Relative Humidity")
    print("  - Heat Index / Apparent Temperature")
    print("  - Humidex")
    print("• Wind speed = sqrt(u² + v²)")
    print("• For daily aggregation, consider using daily mean, max, min")
    print("• Precipitation is in meters (multiply by 1000 for mm)")
    print("• Solar radiation is in J/m² (divide by 3600 for W/m²)")


if __name__ == '__main__':
    download_era5_brazil()
