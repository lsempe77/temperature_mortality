import sidrapy
import pandas as pd

def check_metadata():
    print("Checking metadata for Table 9514 (Census 2022)...")
    
    # We can't easily "list variables" with sidrapy directly without knowing the API structure well,
    # but we can try to fetch a small sample with the standard variable '93' (People) 
    # instead of '1000093' to see if that fixes it.
    
    # Variable 93 is usually "População residente" in absolute numbers.
    # Variable 1000093 might be a percentage or derived indicator.
    
    try:
        print("Attempting fetch with Variable 93 (Absolute Count)...")
        data = sidrapy.get_table(
            table_code="9514",
            territorial_level="6",
            ibge_territorial_code="1100015", # Test with one municipality (Alta Floresta D'Oeste)
            variable="93",
            period="2022",
            classificacao="2/6794|287/652" # Total Sex, Total Age
        )
        
        if data is not None and not data.empty:
            print("Success! Sample data:")
            print(data.iloc[:2])
        else:
            print("Variable 93 returned no data.")
            
    except Exception as e:
        print(f"Error with Var 93: {e}")

if __name__ == "__main__":
    check_metadata()
