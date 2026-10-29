## Data and analysis

`climate_ai_filtered.csv` contains patent-level data on climate technologies, including AI classifications, technology fields, patent market-value estimates, forward citations, and patent characteristics used as controls.

The R analysis replicates the tables and figures of the paper. In particular, it compares private value and citation-based technological influence across stages of technological maturity. It estimates stacked and separate-outcome regressions, examines AI components and citation patterns, and performs robustness checks using alternative samples, AI thresholds, and maturity measures.

Required R packages: `data.table`, `fixest`, `minpack.lm`, `MatchIt`, `ggplot2`, `dplyr`, and `tidyr`.
