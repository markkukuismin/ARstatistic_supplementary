# ARstatistic_supplementary
Supplementary materials for the article "Using rejection sampling probability of acceptance as a measure of independence"

The R package `ARstatistic` can be found at [https://github.com/markkukuismin/ARstatistic](https://github.com/markkukuismin/ARstatistic). DOI: 10.5281/zenodo.12740094

# Author Contributions Checklist Form

# Data

## Abstract
The data set analyzed in the paper was simulated using R. For example, use R scripts `power_analysis.R`, `power_analysis_monotonic.R`, `power_analysi_loess` and `power_analysis_small_sample.R` to reproduce the simulation results discussed in the paper. The real data are available publicly in R packages `NISTnls`, `GGMnonreg`, and `mlbench`.

# Availability
The real data are available publicly in R packages `NISTnls`, `GGMnonreg`, and `mlbench`.

# Code

## Abstract
The R codes provided with the manuscript reproduce all the simulation, Figures, Tables, and real data results discussed in the paper. R source codes of the new method proposed in the paper are found in the folder "functions" ("rho.R", "rho_rank.R", "rho_test_perm.R")

## Description
The R code (R script files with extension .R) are provided at Github.

* `rho.R`: Returns estimates of the A-R statistic (expected value) or simulated distribution.
* `rho_test_perm.R`: Permutation test.
* `rho_rank.R`: Returns estimate of the A-R statistic calculated using ranks.
* `rho_rank_qs.R`: Used only to demonstrate how the A-R statistic using ranks can be computed in $O(n \log n)$ time. Should not be used in real data-analysis!

The following R script files are needed to reproduce Figures and Tables in the paper.

### Main article

Tables and Figures are numbered according to the revision files.

* `Illustrate_different_dependencies/`
  - `Illustrate_different_dependencies_ggplot.R`: reproduce Figure 1 in the paper.
* `power_analysis/`
  - `power_analysis.R`: reproduce the simulation analysis and Figure 3 discussed in Section 3.1 ("Simulation Results").
  - `relative_power.R`: reproduce estimates shown in Table 2 discussed in Section 3.1 ("Simulation Results").
* `power_analysis_loess/`
  - `power_analysis_loess.R`: reproduce the simulation analysis (a simple random walk) discussed in Section 3.1 ("Simulations").
  -`plot_power.R`: reproduces Figure 4 discussed in Section 3.1 ("Simulation Results").
  - `plot_scatterplots.R`: reproduces Figure 2 in Section 3 ("Simulations").
  - `relative_power.R`: reproduce the simulation analysis results (Table 4 and Supplementary Figure S2) of the simple random walk discussed in Section 3.1 ("Simulation Results").
* `power_analysis_monotonic/power_analysis_monotonic.R`: reproduce the simulation analysis and Figure 5 discussed in Section 3.1 ("Simulation Results").
* `ar_mic_rdc_dcor_xicor_under_null/statistic_dist_under_null_and_alt.R`: reproduce Figure 6 in the paper.
* `power_analysis_small_sample/`
  - `power_analysis_small_sample.R`: reproduce the small sample simulation analysis and Supplementary Figure S3 discussed in Section 3.1 ("Simulation Results").
  - `relative_power_small_sample.R`: reproduce estimates shown in Table 6 and Supplementary Figure S4 discussed in Section 3.1 ("Simulation Results").
* `Eckerle4_data_example/Eckerle4_data_demo.R`: reproduce estimates and p-values discussed in Section 4.1 ("Example 1").
* `Sachs_network_example/`
  - `Sachs_network_demo.R`: reproduce estimates, undirected networks, and Figures 7 and 8 discussed in Section 4.2 ("Example 2").
  - `Sachs_dag_graph.R`: Create an undirected version of the directed acyclic graph (dag) reported in Sachs et al (2005).

### Supplementary Materials

* `bootstrap_example/bootstrap_demo.R`: reproduces Supplementary Examples (Example 1, 2, and Supplementary Figures S5, S6).
* `ar_credible_interval/`
  - `ar_vs_cor.R`: reproduces Supplementary Examples (Example 3, and Supplementary Figures S7 and S8).
  - `interval_estimates_under_null.R`: reproduces Supplementary Figure S9.
* `sonar_rf_example/`
  - `sonar_rf_demo.R`: reproduces the Random Forest (RF) example scenario (i), and Supplementary Figures S10 and S11 discussed in Supplementary Materials ("Example 3".)
  - `sonar_rf_demo_test_accuracy.R`: reproduces the RF example scenario (ii), and Supplementary Figure S12 discussed in Supplementary Materials ("Example 3".)
* `running_time_comparison/runtime_comparison.R`: reproduces run time comparison (Supplementary Figure S13).

Additional R scripts.

* `rdc.R`: The Randomized Dependence Coefficient (see [link](https://proceedings.neurips.cc/paper/2013/file/aab3238922bcc25a6f606eb525ffdc56-Paper.pdf))
* `get_legend.R`: For legend placement in Figures 3 and 5.
* `findCorrelation_length.R`: A help function used in the RF example.
* `qsort.R`: Quick Sort algorithm.
* `clusters.R`: Simulate the 5 cluster example data.

# Reproducibility
Tables and Figures in the paper can be reproduced using the supplementary R scripts.
