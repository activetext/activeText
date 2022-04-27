---
title: "activeText_vignette"
output: rmarkdown::pdf_vignette
vignette: >
  %\VignetteIndexEntry{activeText_vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



After installing the package, load it as normal using the `library()` function.


```r
library(activeText)
#> Error in library(activeText): there is no package called 'activeText'
library(tidyverse)
#> ── Attaching packages ──────────────────────────────────────────────── tidyverse 1.3.1 ──
#> ✔ ggplot2 3.3.5     ✔ purrr   0.3.4
#> ✔ tibble  3.1.6     ✔ dplyr   1.0.8
#> ✔ tidyr   1.2.0     ✔ stringr 1.4.0
#> ✔ readr   2.1.2     ✔ forcats 0.5.1
#> ── Conflicts ─────────────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
```

Let's take a look at the dataset that we'll use for our demonstration of the package.

```r
activeText::bbc_data_all
#> Error in loadNamespace(x): there is no package called 'activeText'
```
The BBC news dataset contains 2,215 news articles that are evenly divided between tech, entertainment, sports, politics, and business. All articles that belong to the politics group are labeled as 1, and all other articles are labeled as 0.


```r
results <- activeText::active_label(
                      docs = activeText::bbc_data_all,
                      handlabel = FALSE,
                      labels = c("Not Political", "Political"),
                      labels_name = "label",
                      lambda = 0.01,
                      max_active = 10,
                      init_size = 1,
                      max_query = 1
                    )
#> Error in loadNamespace(x): there is no package called 'activeText'
print(results$docs)
#> Error in print(results$docs): object 'results' not found
```
