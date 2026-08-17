# Bee nesting strategy, sociality, and climate

This repository contains the R scripts used to curate data and conduct comparative analyses of bee nesting strategy, sociality, geographic distributions, climatic niches, and trait evolution.

The analyses combine species-level sociality and nesting traits with a bee phylogeny, occurrence records, and climatic data. The workflow includes trait and occurrence-data curation, geographic and climatic niche analyses, phylogenetic comparative analyses, ancestral-state reconstruction, tests of correlated trait evolution, and models of evolutionary climatic optima.

---

## Scripts

### `00_utility_functions.R`

Contains custom functions used by multiple scripts for data processing, climatic summaries, phylogenetic analyses, and other utility operations.

### `01.1_trait_data_curation.R`

Curates the species-level sociality and nesting trait dataset and corresponding phylogeny. Recodes traits for binary and three-state analyses, handles parasitic species and uncertain trait states, and prunes the phylogeny to match the final trait dataset.

### `01.2_check_distributions.R`

Provides a visual quality-control step for spatially thinned occurrence data by plotting occurrence records separately for each species.

### `01.3_occurrence_data_curation.R`

Curates bee occurrence records and generates species-level climatic summaries. The workflow removes occurrences outside the native ranges of introduced species, spatially thins occurrence records, extracts climatic data, and calculates species-level climate summary statistics. Removal of introduced occurrences includes a manual GIS review step.

### `02_kruskal_wallis.R`

Tests for differences in climatic variables among sociality and nesting groups using Kruskal–Wallis tests and post hoc comparisons and generates associated climatic distribution plots.

### `03_heatmaps.R`

Visualizes geographic variation in bee sociality and nesting strategy across the Americas, including regional trait proportions and their relationships with latitude.

### `04_climatic_niche_breadth.R`

Evaluates climatic differentiation among sociality and nesting groups using principal component analysis (PCA), PERMANOVA, and multivariate climatic niche volumes.

### `05_phylANOVA.R`

Uses phylogenetic ANOVA to test associations between sociality, nesting strategy, and climatic variables while accounting for phylogenetic relationships among species.

### `06_corHMM.R`

Models the joint evolution of sociality and nesting strategy using `corHMM`, including model comparison, transition-rate estimation, and ancestral-state reconstruction.

### `07_corHMM_simmaps.R`

Generates stochastic character maps from the selected `corHMM` model and summarizes evolutionary transitions among sociality × nesting states.

### `08_fitCorrelationTest.R`

Tests for correlated evolution between sociality and nesting strategy using `corHMM::fitCorrelationTest` across alternative trait encodings.

### `09_hOUwie_bio1.R`, `09_hOUwie_bio4.R`, `09_hOUwie_bio12.R`, `09_hOUwie_bio15.R`

Fit evolutionary models relating sociality and nesting regimes to four climatic variables:

- **BIO1:** Mean annual temperature
- **BIO4:** Temperature seasonality
- **BIO12:** Annual precipitation
- **BIO15:** Precipitation seasonality

### `10_retrieve_hOUwie_results.R`

Summarizes and compares fitted hOUwie models, retrieves parameter estimates, and generates figures of estimated evolutionary climatic optima.

### `11_ASR_tree.R`

Generates the circular phylogenetic visualization of the ancestral-state reconstruction, combining reconstructed evolutionary states, observed tip states, a temporal axis, and family-level annotations.
