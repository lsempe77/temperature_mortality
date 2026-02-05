import requests
import pandas as pd

def debug_requests():
    print("Debugging SIDRA API with direct requests...")
    
    # Base URL pattern
    # /t/{table}/n6/{muni}/v/{var}/p/{period}/c{class_id}/{cat_id}/...
    
    # Test 1: Total Population
    # Table 9514, Muni 1100015, Var 93, Period 2022, Sex(2)=Total(6794), Age(287)=Total(652)
    url1 = "https://apisidra.ibge.gov.br/values/t/9514/n6/1100015/v/93/p/2022/c2/6794/c287/652"
    print(f"\nURL 1: {url1}")
    try:
        r = requests.get(url1)
        print(r.json()[0]) # Header
        print(r.json()[1]) # Data
    except Exception as e:
        print(e)

    # Test 2: Elderly Population (60-64)
    # Age(287)=3244
    url2 = "https://apisidra.ibge.gov.br/values/t/9514/n6/1100015/v/93/p/2022/c2/6794/c287/3244"
    print(f"\nURL 2: {url2}")
    try:
        r = requests.get(url2)
        print(r.json()[1]) # Data
    except Exception as e:
        print(e)

    # Test 3: Urban Population
    # Situation(1)=Urban(1)
    # Note: Table 9514 might not have Situation(1). Let's check Table 9514 metadata or try.
    # Actually, Table 9514 is "Population resident by Age, Sex and Form of Age Declaration".
    # It DOES NOT seem to have Urban/Rural (Situation).
    # We need to find the correct table for Urban/Rural.
    # Table 4714 (Census 2022) - Population by Situation?
    # Or Table 9513?
    
    # Let's try adding c1/1 to Table 9514 just in case, but likely it will fail or ignore.
    url3 = "https://apisidra.ibge.gov.br/values/t/9514/n6/1100015/v/93/p/2022/c2/6794/c287/652/c1/1"
    print(f"\nURL 3: {url3}")
    try:
        r = requests.get(url3)
        print(r.json()[1])
    except Exception as e:
        print(e)

if __name__ == "__main__":
    debug_requests()
