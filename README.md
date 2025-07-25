# hespresso: Methods for detecting shifts in homeolog expression ratios in allopolyploids
    
_hespresso_ blends homeolog expression with espresso, capturing
the essence of power, energy, and focus in homeolog expression analysis.
It provides statistical tests for detecting shifts in homeolog expression
ratios among subgenomes of allopolyploid species across diverse conditions,
using RNA-Seq read count data.
Take a sip of _hespresso_, and start your analysis strong!


## Installation

To install the _hespresso_ package,
start R (version 4.5 or later) and run the following commands:

```r
# Bioconductor Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("edgeR"))

# CRAN Packages
install.packages(c("remotes", "foreach", "doParallel", "progressr", "ggplot2"))

# cmdstanr
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
install_cmdstan()

# hespresso
remotes::install_github("bitdessin/hespresso")
```



## Documentation

To list all available vignettes included in the _hespresso_ package,
run as below.

```r
browseVignettes("hespresso")
```


## Citations

If you use HOBIT, please cite:

Sun J, Sese J, and Shimizu KK.
A moderated statistical test for detecting shifts in homeolog expression ratios
in allopolyploids.
bioRxiv 2025;2025.07.01.660977.
[10.1101/2025.07.01.660977](https://doi.org/10.1101/2025.07.01.660977)


If you use HomeoRoq, please cite:

Akama S, Shimizu-Inatsugi R, Shimizu KK, Sese J.
Genome-wide quantification of homeolog expression ratio revealed nonstochastic
gene regulation in synthetic allopolyploid Arabidopsis.
Nucleic Acids Res. 2014;42(6):e46.
[10.1093/nar/gkt1376](https://doi.org/10.1093/nar/gkt1376)

