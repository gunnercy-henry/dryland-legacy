# dryland-legacy
This repository contains the MATLAB code used to quantify productivity legacy effects following hydroclimatic extremes and to estimate biological and water memory across global drylands.

Files

01_quantify_legacy_effects.m
Quantifies event-level and pixel-level legacy effects. The script:

detrends and standardizes productivity and water-availability time series;

identifies isolated wet and dry extreme years;

calculates legacy effects as residual productivity following an extreme;

aligns wet and dry events onto a common directional scale;

classifies sustained, reversed and non-detectable legacy effects using empirical null distributions;

calculates event-level legacy frequencies;

calculates pixel-level mean legacy effects using bootstrap null distributions;

evaluates sensitivity to alternative percentile thresholds and null definitions.

02_separate_biological_abiotic_memory.m
Estimates biological and water memory using partial least squares regression. Previous-year productivity and water availability are related to subsequent legacy effects, with the corresponding coefficients interpreted as biological memory and water memory. The script also evaluates detectable-event-only estimates and lag lengths of 1–3 years.

Input data

The main input file should contain:

NDVI: annual productivity data, arranged as rows × columns × years

PPT: annual water-availability data, arranged as rows × columns × years

AI: aridity index used to define drylands

Optional disturbance data can also be supplied for fire-sensitivity analyses.

Running the analysis

Run the scripts in the following order:

01_quantify_legacy_effects.m

02_separate_biological_abiotic_memory.m

The first script generates the legacy-effect outputs required by the second script.
