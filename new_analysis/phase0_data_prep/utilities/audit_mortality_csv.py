"""
Audit all mortality CSV files to document format variations.
"""
import pandas as pd
from pathlib import Path

INPUT_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\Input_data")

print("="*80)
print("MORTALITY CSV FORMAT AUDIT")
print("="*80)

for year in range(2010, 2025):
    yy = str(year)[-2:]
    file = INPUT_DIR / f"DO{yy}OPEN.csv"
    
    if not file.exists():
        print(f"{year}: FILE NOT FOUND")
        continue
    
    # Read first line to detect separator
    with open(file, 'r', encoding='latin1') as f:
        header_line = f.readline()
    
    # Count separators
    comma_count = header_line.count(',')
    semicolon_count = header_line.count(';')
    sep = ',' if comma_count > semicolon_count else ';'
    
    # Read sample
    try:
        df = pd.read_csv(file, sep=sep, encoding='latin1', nrows=3, dtype=str)
        df.columns = df.columns.str.strip().str.strip('"')
        
        # Get key columns info
        dtobito_sample = df['DTOBITO'].iloc[0] if 'DTOBITO' in df.columns else 'N/A'
        idade_sample = df['IDADE'].iloc[0] if 'IDADE' in df.columns else 'N/A'
        sexo_sample = df['SEXO'].iloc[0] if 'SEXO' in df.columns else 'N/A'
        
        # Detect date format
        if '-' in str(dtobito_sample):
            date_fmt = 'YYYY-MM-DD'
        else:
            date_fmt = 'DDMMYYYY'
        
        # Detect age format
        try:
            age_val = int(str(idade_sample).replace(',', ''))
            if age_val >= 400:
                age_fmt = '4xx coded'
            else:
                age_fmt = 'plain years'
        except:
            age_fmt = 'unknown'
        
        print(f"{year}: sep='{sep}'  date={date_fmt:12}  age={age_fmt:12}  DTOBITO={dtobito_sample}, IDADE={idade_sample}, SEXO={sexo_sample}")
    except Exception as e:
        print(f"{year}: ERROR - {e}")

print("="*80)
