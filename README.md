---
author: Mitchell Bosley
title: README
---

R package for using active learning to classify text documents.

Authors:

-   [Mitchell Bosley](https://mbosley.github.io)
-   [Saki Kuzushima](https://ksaki.github.io)
-   [Ted Enamorado](https://www.tedenamorado.com/)
-   [Yuki Shiraito](https://shiraito.github.io)

# Installation Instructions

To install `activeText`, you must first install the `devtools` using the
following code:

``` r
if(!require(devtools)) install.packages("devtools")
```

Then, load \`devtools\` and use the function
\`install<sub>github</sub>()\` to install \`activeText\`:

``` r
library(devtools)
install_github("activetext/activeText", dependencies = TRUE)
```
