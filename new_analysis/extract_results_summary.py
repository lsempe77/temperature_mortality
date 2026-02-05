"""
Extract and compare results across all phases for intermediate and immediate levels.
"""
import json
import os
from pathlib import Path

BASE_DIR = Path(__file__).parent

def load_json(path):
    try:
        with open(path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading {path}: {e}")
        return {}

# Load all results
results = {
    'dlnm_int': load_json(BASE_DIR / 'phase1_r/results/dlnm_r_intermediate_results_v2.json'),
    'dlnm_imm': load_json(BASE_DIR / 'phase1_r/results/dlnm_r_immediate_results_v2.json'),
    'burden_int': load_json(BASE_DIR / 'phase1_r/results/attributable_burden_r_intermediate.json'),
    'burden_imm': load_json(BASE_DIR / 'phase1_r/results/attributable_burden_r_immediate.json'),
    'sens_int': load_json(BASE_DIR / 'phase2_r/results/sensitivity_r_intermediate.json'),
    'sens_imm': load_json(BASE_DIR / 'phase2_r/results/sensitivity_r_immediate.json'),
    'harv_int': load_json(BASE_DIR / 'phase2_r/results/harvesting_r_intermediate.json'),
    'harv_imm': load_json(BASE_DIR / 'phase2_r/results/harvesting_r_immediate.json'),
    'heatwave_int': load_json(BASE_DIR / 'phase2_r/results/heatwave_r_intermediate.json'),
    'heatwave_imm': load_json(BASE_DIR / 'phase2_r/results/heatwave_r_immediate.json'),
    'supp_int': load_json(BASE_DIR / 'phase3_r/results/supplementary_r_intermediate.json'),
    'supp_imm': load_json(BASE_DIR / 'phase3_r/results/supplementary_r_immediate.json'),
    'meta_int': load_json(BASE_DIR / 'phase4_r/results/meta_regression_intermediate.json'),
    'meta_imm': load_json(BASE_DIR / 'phase4_r/results/meta_regression_immediate.json'),
    'age_int': load_json(BASE_DIR / 'phase4_r/results/age_stratification_intermediate.json'),
    'age_imm': load_json(BASE_DIR / 'phase4_r/results/age_stratification_immediate.json'),
    'sex_int': load_json(BASE_DIR / 'phase4_r/results/sex_stratification_intermediate.json'),
    'sex_imm': load_json(BASE_DIR / 'phase4_r/results/sex_stratification_immediate.json'),
    'cause_int': load_json(BASE_DIR / 'phase4_r/results/cause_stratification_intermediate.json'),
    'cause_imm': load_json(BASE_DIR / 'phase4_r/results/cause_stratification_immediate.json'),
}

output_lines = []

def p(text=""):
    print(text)
    output_lines.append(text)

p("=" * 80)
p("COMPARATIVE ANALYSIS: TEMPERATURE-MORTALITY IN BRAZIL (2010-2024)")
p("Immediate (510 regions) vs Intermediate (133 regions) Geographic Levels")
p("=" * 80)
p()

# Phase 1: DLNM
p("=" * 80)
p("PHASE 1: DISTRIBUTED LAG NON-LINEAR MODEL (DLNM) - CORE RESULTS")
p("=" * 80)
p()

for level in ['int', 'imm']:
    name = 'INTERMEDIATE (133 regions)' if level == 'int' else 'IMMEDIATE (510 regions)'
    d = results[f'dlnm_{level}'].get('pooled', {})
    heat = d.get('rr_heat_p99', {})
    cold = d.get('rr_cold_p1', {})
    heat95 = d.get('rr_heat_p95', {})
    cold5 = d.get('rr_cold_p5', {})
    
    p(f"{name}:")
    p(f"  Pooled regions: {d.get('n_regions', 'N/A')}")
    p(f"  MMT (Minimum Mortality Temperature): {d.get('mmt', 0):.1f}°C")
    p(f"  Extreme Heat (P99): RR = {heat.get('rr', 0):.3f} (95% CI: {heat.get('rr_lo', 0):.3f}-{heat.get('rr_hi', 0):.3f})")
    p(f"  Moderate Heat (P95): RR = {heat95.get('rr', 0):.3f} (95% CI: {heat95.get('rr_lo', 0):.3f}-{heat95.get('rr_hi', 0):.3f})")
    p(f"  Extreme Cold (P1):  RR = {cold.get('rr', 0):.3f} (95% CI: {cold.get('rr_lo', 0):.3f}-{cold.get('rr_hi', 0):.3f})")
    p(f"  Moderate Cold (P5):  RR = {cold5.get('rr', 0):.3f} (95% CI: {cold5.get('rr_lo', 0):.3f}-{cold5.get('rr_hi', 0):.3f})")
    p()

# Phase 1: Attributable Burden
p("=" * 80)
p("PHASE 1: ATTRIBUTABLE MORTALITY BURDEN")
p("=" * 80)
p()

for level in ['int', 'imm']:
    name = 'INTERMEDIATE' if level == 'int' else 'IMMEDIATE'
    d = results[f'burden_{level}']
    
    p(f"{name}:")
    p(f"  Study period: 2010-2024 (15 years)")
    p(f"  Total elderly deaths: {d.get('total_deaths', 0):,}")
    p(f"  Heat-attributable deaths: {d.get('heat_attributable_deaths', 0):,.0f} ({d.get('heat_af_pct', 0):.2f}% AF)")
    p(f"  Cold-attributable deaths: {d.get('cold_attributable_deaths', 0):,.0f} ({d.get('cold_af_pct', 0):.2f}% AF)")
    p(f"  Total temperature-attributable: {d.get('total_attributable_deaths', 0):,.0f} ({d.get('total_af_pct', 0):.2f}% AF)")
    annual_heat = d.get('heat_attributable_deaths', 0) / 15
    annual_cold = d.get('cold_attributable_deaths', 0) / 15
    p(f"  Annual average: Heat={annual_heat:,.0f}, Cold={annual_cold:,.0f}")
    p()

# Phase 2: Sensitivity
p("=" * 80)
p("PHASE 2: SENSITIVITY ANALYSES - LAG STRUCTURE")
p("=" * 80)
p()

for level in ['int', 'imm']:
    name = 'INTERMEDIATE' if level == 'int' else 'IMMEDIATE'
    d = results[f'sens_{level}']
    
    p(f"{name}:")
    lag_sens = d.get('lag_sensitivity', {})
    for lag_key in ['lag_7', 'lag_14', 'lag_21', 'lag_28']:
        lag_data = lag_sens.get(lag_key, {})
        heat = lag_data.get('heat', {})
        cold = lag_data.get('cold', {})
        lag_days = lag_key.replace('lag_', '')
        if heat.get('rr'):
            p(f"  Lag {lag_days} days: Heat RR={heat.get('rr', 0):.3f} ({heat.get('rr_lo', 0):.3f}-{heat.get('rr_hi', 0):.3f}), Cold RR={cold.get('rr', 0):.3f} ({cold.get('rr_lo', 0):.3f}-{cold.get('rr_hi', 0):.3f})")
    p()

# Phase 2: Heatwave
p("=" * 80)
p("PHASE 2: HEATWAVE EFFECT MODIFICATION")
p("=" * 80)
p()

for level in ['int', 'imm']:
    name = 'INTERMEDIATE' if level == 'int' else 'IMMEDIATE'
    d = results[f'heatwave_{level}']
    
    p(f"{name}:")
    hw = d.get('heatwave_effect', {})
    p(f"  Heatwave definition: ≥3 consecutive days above P95")
    p(f"  Heatwave days: {d.get('heatwave_days', 'N/A')} ({d.get('heatwave_pct', 0):.2f}% of total)")
    p(f"  Regions with heatwaves: {d.get('regions_with_heatwaves', 'N/A')}")
    p(f"  Additive heatwave effect: RR = {hw.get('rr', 0):.3f} (95% CI: {hw.get('rr_lo', 0):.3f}-{hw.get('rr_hi', 0):.3f})")
    sig = "SIGNIFICANT" if hw.get('rr_lo', 0) > 1 else "NOT SIGNIFICANT"
    p(f"  Interpretation: {sig} - Multi-day heat events have {'additional' if hw.get('rr_lo', 0) > 1 else 'no additional'} mortality risk")
    p()

# Phase 3: Confounding
p("=" * 80)
p("PHASE 3: CONFOUNDING CONTROL - MODEL SPECIFICATIONS")
p("=" * 80)
p()

for level in ['int', 'imm']:
    name = 'INTERMEDIATE' if level == 'int' else 'IMMEDIATE'
    d = results[f'supp_{level}']
    
    p(f"{name}:")
    models = d.get('models', {})
    for model_name, model_data in models.items():
        heat = model_data.get('heat', {})
        cold = model_data.get('cold', {})
        if heat.get('rr'):
            p(f"  {model_name}: Heat RR={heat.get('rr', 0):.3f}, Cold RR={cold.get('rr', 0):.3f}")
    p()

# Phase 4: Age Stratification
p("=" * 80)
p("PHASE 4: HETEROGENEITY - AGE STRATIFICATION")
p("=" * 80)
p()

for level in ['int', 'imm']:
    name = 'INTERMEDIATE' if level == 'int' else 'IMMEDIATE'
    d = results[f'age_{level}']
    
    p(f"{name}:")
    age_groups = d.get('age_groups', {})
    for age_name, age_data in age_groups.items():
        heat = age_data.get('heat', {})
        cold = age_data.get('cold', {})
        if heat.get('rr'):
            p(f"  {age_name}: Heat RR={heat.get('rr', 0):.3f} ({heat.get('rr_lo', 0):.3f}-{heat.get('rr_hi', 0):.3f}), Cold RR={cold.get('rr', 0):.3f} ({cold.get('rr_lo', 0):.3f}-{cold.get('rr_hi', 0):.3f})")
    p()

# Phase 4: Sex Stratification
p("=" * 80)
p("PHASE 4: HETEROGENEITY - SEX STRATIFICATION")
p("=" * 80)
p()

for level in ['int', 'imm']:
    name = 'INTERMEDIATE' if level == 'int' else 'IMMEDIATE'
    d = results[f'sex_{level}']
    
    p(f"{name}:")
    sex_groups = d.get('sex_groups', {})
    for sex_name, sex_data in sex_groups.items():
        heat = sex_data.get('heat', {})
        cold = sex_data.get('cold', {})
        if heat.get('rr'):
            p(f"  {sex_name}: Heat RR={heat.get('rr', 0):.3f} ({heat.get('rr_lo', 0):.3f}-{heat.get('rr_hi', 0):.3f}), Cold RR={cold.get('rr', 0):.3f} ({cold.get('rr_lo', 0):.3f}-{cold.get('rr_hi', 0):.3f})")
    p()

# Phase 4: Cause Stratification
p("=" * 80)
p("PHASE 4: HETEROGENEITY - CAUSE OF DEATH STRATIFICATION")
p("=" * 80)
p()

for level in ['int', 'imm']:
    name = 'INTERMEDIATE' if level == 'int' else 'IMMEDIATE'
    d = results[f'cause_{level}']
    
    p(f"{name}:")
    cause_groups = d.get('cause_groups', {})
    for cause_name, cause_data in cause_groups.items():
        heat = cause_data.get('heat', {})
        cold = cause_data.get('cold', {})
        if heat.get('rr'):
            p(f"  {cause_name}: Heat RR={heat.get('rr', 0):.3f} ({heat.get('rr_lo', 0):.3f}-{heat.get('rr_hi', 0):.3f}), Cold RR={cold.get('rr', 0):.3f} ({cold.get('rr_lo', 0):.3f}-{cold.get('rr_hi', 0):.3f})")
    p()

# Meta-regression
p("=" * 80)
p("PHASE 4: META-REGRESSION - REGIONAL HETEROGENEITY")
p("=" * 80)
p()

for level in ['int', 'imm']:
    name = 'INTERMEDIATE' if level == 'int' else 'IMMEDIATE'
    d = results[f'meta_{level}']
    
    p(f"{name}:")
    p(f"  Regions fitted: {d.get('n_regions_fitted', 'N/A')}")
    p(f"  Regions with valid estimates: {d.get('n_regions_valid', 'N/A')}")
    het = d.get('heterogeneity', {})
    p(f"  I² (heterogeneity): {het.get('i2', 'N/A')}")
    p()

# Save to file
output_path = BASE_DIR / 'phase5_outputs' / 'results_comparison_analysis.md'
output_path.parent.mkdir(exist_ok=True)

with open(output_path, 'w') as f:
    f.write('\n'.join(output_lines))

print(f"\nResults saved to: {output_path}")
