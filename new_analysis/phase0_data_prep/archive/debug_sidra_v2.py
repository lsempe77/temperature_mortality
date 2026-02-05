import sidrapy
import pandas as pd

def debug_sidra():
    print("Debugging SIDRA API calls...")
    
    # Test 1: Total Population
    print("\n--- Test 1: Total Population ---")
    try:
        data = sidrapy.get_table(
            table_code="9514",
            territorial_level="6",
            ibge_territorial_code="1100015", # Alta Floresta D'Oeste
            variable="93",
            period="2022",
            classificacao="2/6794|287/652" # Total Sex, Total Age
        )
        print(data.iloc[:5] if data is not None else "No data")
    except Exception as e:
        print(f"Error: {e}")

    # Test 2: Elderly Population (Specific Age Group)
    print("\n--- Test 2: Elderly Population (60-64) ---")
    try:
        data = sidrapy.get_table(
            table_code="9514",
            territorial_level="6",
            ibge_territorial_code="1100015",
            variable="93",
            period="2022",
            classificacao="2/6794|287/3244" # Total Sex, Age 60-64
        )
        print(data.iloc[:5] if data is not None else "No data")
    except Exception as e:
        print(f"Error: {e}")

    # Test 3: Urban Population
    print("\n--- Test 3: Urban Population ---")
    try:
        data = sidrapy.get_table(
            table_code="9514",
            territorial_level="6",
            ibge_territorial_code="1100015",
            variable="93",
            period="2022",
            classificacao="2/6794|287/652|1/1" # Total Sex, Total Age, Urban
        )
        print(data.iloc[:5] if data is not None else "No data")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    debug_sidra()
