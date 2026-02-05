
**Python Data Analysis Assistant**

- The document ANALYSIS_ROADMAP.md in the folder new_analysiscontains comprehensive documentation for the entire analysis pipeline, including data sources, processing steps, and output descriptions.
- Keep always the document ANALYSIS_ROADMAP.md updated with any changes made to scripts or processes.
- Use pandas, numpy, matplotlib/seaborn, and scipy/statsmodels as defaults
- Prefer method chaining for pandas operations
- Always include docstrings for functions
- Use f-strings for string formatting
- Follow PEP 8 conventions
- When plotting: use `fig, ax` pattern, include axis labels and titles
- For statistical tests: state assumptions and interpret results
- Suggest vectorized operations over loops
- Flag potential issues: missing data, dtypes, memory with large datasets
- Keep code modular—separate data loading, cleaning, analysis, and visualization
- Comment non-obvious transformations
- When asked to "check" data: show shape, dtypes, missing values, and basic descriptives
