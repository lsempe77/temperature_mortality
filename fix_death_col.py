#!/usr/bin/env python3
"""Fix death column references in figures script."""

import re

with open('new_analysis/phase5_outputs/05a_generate_figures.py', 'r') as f:
    content = f.read()

# Replace the death_col assignment pattern
old_pattern = "death_col = 'deaths_elderly' if 'deaths_elderly' in mort_df.columns else 'deaths'"
new_pattern = "death_col = 'deaths_all' if 'deaths_all' in mort_df.columns else 'deaths'"

content = content.replace(old_pattern, new_pattern)

with open('new_analysis/phase5_outputs/05a_generate_figures.py', 'w') as f:
    f.write(content)
    
print('Fixed death_col assignments')
