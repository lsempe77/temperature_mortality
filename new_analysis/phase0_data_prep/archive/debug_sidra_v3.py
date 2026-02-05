import requests
import pandas as pd

def debug_requests_v3():
    print("Debugging SIDRA API v3...")
    
    # Test 1: GDP (Table 5938) - Should work
    # Var 37 = PIB
    url_gdp = "https://apisidra.ibge.gov.br/values/t/5938/n6/1100015/v/37/p/2021"
    print(f"\nGDP URL: {url_gdp}")
    try:
        r = requests.get(url_gdp)
        print(r.json()[1])
    except Exception as e:
        print(e)

    # Test 2: Urban/Rural (Table 4714 - Census 2022)
    # Var 93 = População residente
    # C1 = Situação do domicílio (1=Urbana, 2=Rural, 0=Total?)
    # Let's try c1/1 (Urban) and c1/2 (Rural)
    url_urban = "https://apisidra.ibge.gov.br/values/t/4714/n6/1100015/v/93/p/2022/c1/1"
    print(f"\nUrban URL: {url_urban}")
    try:
        r = requests.get(url_urban)
        print(r.json()[1])
    except Exception as e:
        print(f"Failed: {e}")

    # Test 3: Age (Table 9514) - Try to get ANY value
    # Maybe we need to specify 'Forma de declaração' explicitly?
    # c10013/113635 (Total)
    url_age = "https://apisidra.ibge.gov.br/values/t/9514/n6/1100015/v/93/p/2022/c2/6794/c287/652/c10013/113635"
    print(f"\nAge URL: {url_age}")
    try:
        r = requests.get(url_age)
        print(r.json()[1])
    except Exception as e:
        print(f"Failed: {e}")

if __name__ == "__main__":
    debug_requests_v3()
