![img](./logo.png)


# Table of Contents

1.  [Overview](#org79cefbb)
2.  [Installation](#org4d5011e)
3.  [Basic Usage](#org43ecad4)


<a id="org79cefbb"></a>

# Overview

R package for using active learning to classify text documents.

Authors:

-   [Mitchell Bosley](https://mbosley.github.io)
-   [Saki Kuzushima](https://sakikuzushima.github.io/)
-   [Ted Enamorado](https://www.tedenamorado.com/)
-   [Yuki Shiraito](https://shiraito.github.io)

Recommended citation: Mitchell Bosley, Saki Kuzushima, Ted Enamorado, Yuki Shiraito. 2022. [“Improving Probabilistic Models in Text Classification via Active Learning”](https://arxiv.org/abs/2202.02629). Preprint, arXiv:2202.02629.




<a id="org4d5011e"></a>
# Installation

To install `activeText`, you must first install the `devtools` using the following code:

```R
if(!require(devtools)) install.packages("devtools")
```

Then, use the function `devtools::install_github()` to install `activeText`:

```R
devtools::install_github("activetext/activeText", ref="main")
```


<a id="org43ecad4"></a>

# Basic Usage

For this tutorial, we&rsquo;re going to use the included BBC News dataset, which is comprised of 2225 documents, 417 of which belong to the &ldquo;politics&rdquo; class. Let&rsquo;s take a quick look at the data:

```R
activeText::bbc_data_all
```

    # A tibble: 2,225 × 3
       id    text                                                              label
       <chr> <chr>                                                             <dbl>
     1 1     "Ink helps drive democracy in Asia\n\nThe Kyrgyz Republic, a sma…     0
     2 2     "China net cafe culture crackdown\n\nChinese authorities closed …     0
     3 3     "Microsoft seeking spyware trojan\n\nMicrosoft is investigating …     0
     4 4     "Digital guru floats sub-$100 PC\n\nNicholas Negroponte, chairma…     0
     5 5     "Technology gets the creative bug\n\nThe hi-tech and the arts wo…     0
     6 6     "Wi-fi web reaches farmers in Peru\n\nA network of community com…     0
     7 7     "Microsoft releases bumper patches\n\nMicrosoft has warned PC us…     0
     8 8     "Virus poses as Christmas e-mail\n\nSecurity firms are warning a…     0
     9 9     "Apple laptop is 'greatest gadget'\n\nThe Apple Powerbook 100 ha…     0
    10 10    "Google's toolbar sparks concern\n\nSearch engine firm Google ha…     0
    # … with 2,215 more rows

For validation purposes, for each document the label `0` indicates a non-political news article, and a label `1` indicates a political news article.

We&rsquo;re going to use the package&rsquo;s workhorse function, `active_label` to use an active learning loop to classify each of these documents as political or not. The simplest possible implementation is to supply the document corpus as the `docs` variable, and the discrete set of possible labels as a vector to the `labels` variable.

```R
results <- activeText::active_label(
                      docs = activeText::bbc_data_all,
                      labels = c("Not Political", "Political"),
                      seed = 1234
                    )
```

A window as shown below will appear. If it does not pop up, check your task bar for an icon wthat looks like a blue feather and click on that. By default, the `active_label` function asks the user to label 10 documents per iteration, up to a maximum of 5 iterations, plus a set of documents for initialization of the algorithm. 

The count of the iteration is at the top left of the application. The number of documents is centered right above the document text. There are radio buttons below the document text to label the document. Press save to save the progress so far, next to go to the next document, and back to go to the prior document.


![The first document to label.](/data/FirstImage.PNG)

After reading the document, the user must then enter select a radio button to continue.

The user will continue to be prompted for labels until the number of labels per iteration is reached. The window will briefly disappear while the EM algorithm runs.

Once the maximum iterations for labeling is reached, an Expectation Maximization (EM) algorithm will run until log-likelihood of the posterior distribution of the parameters stops increasing. The progress of the EM algorithm is tracked by a progress bar in the terminal:

    Active Iter: 1/5 EM Runs: [====>--------------------------------] 14/100 (max)

The numerator of the number on the left hand side indicates the current iteration of the active learning process, while the denominator indicates the maximum number of active learning steps allowed. Similarly, the progress bar and the number on the right indicate the number of iterations of the EM algorithm in a given active step.

Once the maximum number of active iterations have run, or a stopping condition is met (such as running out of documents to label), the process ends. We can inspect the classification results as follows.

```R
results$docs
```

    # A tibble: 2,225 × 5
       id    text                                            label Class_1   Class_2
       <chr> <chr>                                           <dbl>   <dbl>     <dbl>
     1 1     "Ink helps drive democracy in Asia\n\nThe Kyrg…     0       0 1   e+  0
     2 2     "China net cafe culture crackdown\n\nChinese a…     0       1 1.36e- 28
     3 3     "Microsoft seeking spyware trojan\n\nMicrosoft…     0       1 1.56e- 61
     4 4     "Digital guru floats sub-$100 PC\n\nNicholas N…     0       1 1.10e- 59
     5 5     "Technology gets the creative bug\n\nThe hi-te…     0       1 1.69e-109
     6 6     "Wi-fi web reaches farmers in Peru\n\nA networ…     0       1 2.14e- 63
     7 7     "Microsoft releases bumper patches\n\nMicrosof…     0       1 5.92e- 51
     8 8     "Virus poses as Christmas e-mail\n\nSecurity f…     0       1 4.12e- 62
     9 9     "Apple laptop is 'greatest gadget'\n\nThe Appl…     0       1 9.90e-205
    10 10    "Google's toolbar sparks concern\n\nSearch eng…     0       1 2.43e- 90
    # … with 2,215 more rows

The variables `Class_1` and `Class_2` represent the estimated probability that a given document belongs to the &ldquo;Political&rdquo; class or not, respectively. To evaluate the classification accuracy, we can simply calculate the proportion of correctly labeled documents, using a simple cutoff of 0.5. We can also filter out the user labeled documents using `results$hand_labeled_index`, which is vector of the ids for documents that have been user-labeled.

```R
results$docs |>
  dplyr::filter(!(id %in% results$hand_labeled_index)) |>
  dplyr::mutate(
           est_label = ifelse(Class_2 >= 0.5, 1, 0),
           est_correct = ifelse(label == est_label, 1, 0)
         ) |>
  dplyr::pull(est_correct) |>
  sum() / (nrow(results$docs) - length(results$hand_labeled_index))
```

    [1] 0.7524249

Thus, after labeling 60 documents using an active learning process, we are correctly classifying 75 percent of the remaining documents.

If the user needs to take a break at any point during the labeling process, the user selects the save button in the window then exits the window. The save file location is printed to the terminal. To return where the user left off, run (with the real name of the saved file):

```R
activeText::active_label_from_saved("NAME_OF_SAVED_FILE.RDS")
```
