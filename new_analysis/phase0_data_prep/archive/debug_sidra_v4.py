import requests

def debug_requests_v4():
    print("Debugging SIDRA API v4...")
    
    # Test 1: Table 9513 (Census 2022 - Age/Sex)
    # Try fetching Total
    url_9513 = "https://apisidra.ibge.gov.br/values/t/9513/n6/1100015/v/93/p/2022/c2/6794/c287/652"
    print(f"\nTable 9513 URL: {url_9513}")
    try:
        r = requests.get(url_9513)
        print(r.json()[1])
    except Exception as e:
        print(f"Failed: {e}")

    # Test 2: Table 4714 (Census 2022 - Situation)
    # Try fetching Total (no classification)
    url_4714 = "https://apisidra.ibge.gov.br/values/t/4714/n6/1100015/v/93/p/2022"
    print(f"\nTable 4714 URL: {url_4714}")
    try:
        r = requests.get(url_4714)
        print(r.json()[1])
    except Exception as e:
        print(f"Failed: {e}")

    # Test 3: Table 200 (Census 2010 - Situation) - Fallback
    url_200 = "https://apisidra.ibge.gov.br/values/t/200/n6/1100015/v/93/p/2010/c1/1"
    print(f"\nTable 200 (2010) URL: {url_200}")
    try:
        r = requests.get(url_200)
        print(r.json()[1])
    except Exception as e:
        print(f"Failed: {e}")

if __name__ == "__main__":
    debug_requests_v4()
