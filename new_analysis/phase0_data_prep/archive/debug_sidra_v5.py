import requests
import pandas as pd

def debug_requests_v5():
    print("Debugging SIDRA API v5...")
    
    # Test 1: List all Age Groups in Table 9514
    # We use 'all' for c287
    # We must specify c10013 (Forma de declaração) = Total (113635)
    url_age_list = "https://apisidra.ibge.gov.br/values/t/9514/n6/1100015/v/93/p/2022/c2/6794/c287/all/c10013/113635"
    print(f"\nAge List URL: {url_age_list}")
    try:
        r = requests.get(url_age_list)
        data = r.json()
        print(f"Returned {len(data)} rows.")
        # Print first 5 rows to see codes
        for row in data[:5]:
            print(f"Age Code: {row.get('D5C')}, Age Name: {row.get('D5N')}, Value: {row.get('V')}")
    except Exception as e:
        print(f"Failed: {e}")

if __name__ == "__main__":
    debug_requests_v5()
