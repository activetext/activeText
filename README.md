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
-   [Saki Kuzushima](https://ksaki.github.io)
-   [Ted Enamorado](https://www.tedenamorado.com/)
-   [Yuki Shiraito](https://shiraito.github.io)

Recommended citation: Mitchell Bosley, Saki Kuzushima, Ted Enamorado, Yuki Shiraito. 2022. [“Improving Probabilistic Models in Text Classification via Active Learning”](https://arxiv.org/abs/2202.02629). arXiv preprint arXiv:2202.02629




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

First, start by setting a seed for reproducability.

```R
set.seed(1990)
```

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
                    )
```

By default, the `active_label` function asks the user to label 10 documents per iteration, up to a maximum of 5 iterations, plus a set of documents for initialization of the algorithm. The count of the iteration and the number of documents is indicated in the first two lines of the prompt.

    [ Iteration 0 of max 5 ]
    [ Document 1 of 10 ]
    
    GM in crunch talks on Fiat future
    
    Fiat will meet car giant General Motors (GM) on Tuesday in an attempt to reach agreement over the future of the Italian firm's loss-making auto group.
    
    Fiat claims that GM is legally obliged to buy the 90% of the car unit it does not already own; GM says the contract, signed in 2000, is no longer valid. Press reports have speculated that Fiat may be willing to accept a cash payment in return for dropping its claim. Both companies want to cut costs as the car industry adjusts to waning demand.
    
    The meeting between Fiat boss Sergio Marchionne and GM's Rick Wagoner is due to take place at 1330 GMT in Zurich, according to the Reuters news agency.
    
    Mr Marchionne is confident of his firm's legal position, saying in an interview with the Financial Times that GM's argument "has no legs". The agreement in question dates back to GM's decision to buy 20% of Fiat's auto division in 2000. At the time, it gave the Italian firm the right, via a 'put option', to sell the remaining stake to GM. In recent weeks, Fiat has reiterated its claims that this 'put' is still valid and legally binding. However, GM argues that a Fiat share sale made last year, which cut GM's holding to 10%, together with asset sales made by Fiat have terminated the agreement.
    
    Selling the Fiat's car-making unit may not prove so simple, analysts say, especially as it is a company that is so closely linked to Italy's industrial heritage. Political and public pressure may well push the two firms to reach a compromise. "We are not expecting Fiat to exercise its put of the auto business against an unwilling GM at this point," brokerage Merrill Lynch said in a note to investors, adding that any legal battle would be protracted and damaging to the business. "As far as we are aware, the Agnelli family, which indirectly controls at least 30% of Fiat, has not given a firm public indication that it wants to sell the auto business. "Fiat may be willing to cancel the 'put' in exchange for money."
    
    
    1: Not Political
    2: Political
    
    Selection: Enter an item from the menu, or 0 to exit
    Selection: _

After reading the above, the user must then enter either `1` for Not Political, `2` for Political, or `0` to exit.

The user will continue to be prompted for labels until the number of labels per iteration is reached. Once this occurs, an Expectation Maximization (EM) algorithm will run until log-likelihood of the posterior distribution of the parameters stops increasing. The progress of the EM algorithm is tracked by a progress bar:

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
