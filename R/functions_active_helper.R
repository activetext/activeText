##################################################
## Project: active
## Script purpose: Helper Functions for Active Learning EM
## Date: 2020/5/22
## Author: Mitchell Bosley
##################################################

#' @importFrom dplyr "%>%"

#' @export
clean_data <- function(docs, n_class, doc_name, index_name, labels_name=NULL,
                       filters=NULL, add_index=T, add_filter=T, keep_labels=F) {

    #' @title           Structure Data
    #' @description     Structures data to prepare for Active-EM implementation.
    #'                  Options to filter documents by chosen character strings, as well as to
    #'                  add index value for each document.
    #'
    #' @param docs       [matrix]      Matrix of labeled and/or unlabeled documents.
    #' @param n_class    [numeric]     Number of classes to be considered.
    #' @param doc_name    [string]      Character string indicating the variable in 'docs'
    #'                                 that denotes the text of the documents to be classified.
    #' @param index_name  [character]   Character string indicating the variable in 'docs' that
    #'                                 denotes the index value of the document to be classified.
    #' @param labels_name [character]   Character string indicating the variable in \code{docs}
    #'                                 that denotes the already known labels of the documents.
    #'                                 By default, value is set to \code{NULL}.
    #' @param filters    [character]   A vector of regular expressions used to filter out unwanted documents.
    #' @param add_index   [logical]     Boolean logical value indicating whether or not add an index
    #'                                 in the restructuring process.
    #' @param add_filter  [logical]     Boolean logical value indicating whether or not to filter documents
    #'                                 in the restructuring process.
    #' @param keep_labels [logical]     Boolean logical value indicating whether or not to keep an existing column
    #'                                 of labels in the dataset.
    #'
    #' @return           [matrix]      Structured matrix of labeled and unlabeled documents, updated with
    #'                                 labels for the documents in 'toLabel'.

    # Turn documents into tibble
    docs <- dplyr::as_tibble(docs)

    # Data cleaning
    if (keep_labels == T) {
        if (add_index == T) {
            docs <- docs %>%
                dplyr::select(doc_name, labels_name)
        } else {
            docs <- docs %>%
                dplyr::select(doc_name, labels_name, index_name)
        }
    } else {
        if (add_index == T) {
            docs <- docs %>%
                dplyr::select(doc_name)
        } else {
            docs <- docs %>%
                dplyr::select(doc_name, index_name)
        }
    }

    # Filter documents
    if (add_filter == T) {
        docs <- docs %>%
            dplyr::filter(!stringr::str_detect(!!dplyr::sym(doc_name),
                                               paste(filters, collapse="|")))
    }

    # Add index
    if (add_index == T) {
        docs <- docs %>%
            dplyr::mutate(!!dplyr::sym(index_name) := 1:nrow(.))
    }

    # Add classes dynamically
    for (i in 1:n_class) {
        docs <- docs %>%
            tibble::add_column(!!paste("Class", i, sep="_") := NA)
    }

    return(docs)
}

#' @export
query_label <- function(docs, label_id_vec, n_class, labels, doc_name,
                        index_name, labels_name=NULL,
                        active_iter=NULL, maxIter=NULL,
                        handlabel=TRUE, metadata_vars = NA) {

#' @title           Label Query
#' @description     Queries documents for classification by oracle.
#'
#' @param docs       [matrix]    Matrix of labeled and/or unlabeled documents.
#' @param label_id_vec [vector]    Matrix of documents to be labeled.
#' @param n_class    [numeric]   Number of classes to be considered.
#' @param labels     [vector]    Vector of character strings indicating classification options.
#' @param doc_name    [character] Character string indicating the variable in \code{docs} that
#'                               denotes the text of the documents to be classified.
#' @param index_name  [character] Character string indicating the variable in 'docs' that
#'                               denotes the index value of the document to be classified.
#' @param labels_name [character] Character string indicating the variable in \code{docs}
#'                               that denotes the already known labels of the documents.
#'                               By default, value is set to \code{NULL}.
#' @param active_iter [numeric]   Numeric value denoting which iteration of the active learning cycle
#'                               the algorithm is in. Appears as header information to the user-labeling
#'                               process.
#' @param maxIter    [numeric]   Numeric value denoting the maximum number of active learning iterations.
#' @param handlabel  [logical]   Boolean logical value indicating whether to initiate user-input script.
#'                               If set to \code{FALSE}, and if \code{labels_name} is provided, the script
#'                               queries the document label directly from the column denoted by \code{labels_name}.
#'
#' @return          [matrix]     Structured matrix of labeled and unlabeled documents, updated with
#'                               labels for the documents in `label_id_vec`.

  if (handlabel) {

    for (i in 1:length(label_id_vec)) {

      lab_row <- docs %>%
        dplyr::filter(!!dplyr::sym(index_name) == label_id_vec[i])
      lab <- ifelse(
        is.null(labels_name), NA,
        lab_row %>% pull(!!dplyr::sym(labels_name))
      )

      if (is.na(lab)) {

        ## Document tracker
        header <- paste("[ Document", i, "of", length(label_id_vec), "]")


        ## Active learning iteration tracker
        if (!is.null(active_iter) & !is.null(maxIter)) {
          active_header <- paste("\n[ Iteration", active_iter, "of max", maxIter, "]")
          header <- paste(active_header, header, sep = "\n")
        }

        to_label_text <- docs %>%
          dplyr::filter(!!dplyr::sym(index_name) == label_id_vec[i]) %>%
          dplyr::pull(!!dplyr::sym(doc_name))

        ## Menu-based classification
        if (!is.na(metadata_vars)) {
          docs %>%
            filter(!!dplyr::sym(index_name) == label_id_vec[i]) %>%
            select(metadata_vars) %>%
            glimpse()
        }

        selection <- menu(labels, title=paste(header, to_label_text, sep="\n\n"))
        ident <- which(docs[[index_name]] == label_id_vec[i])

      } else {

        selection <- lab
        ident <- which(docs[[index_name]] == label_id_vec[i])
      }

        ## Update document matrix based on classifications
        for (j in 1:n_class){
          docs[ident, paste("Class", j, sep="_") ] <- 0
        }
        docs[ident, paste("Class", selection, sep="_")] <- 1

    }

  } else {

    id_vec <- which(docs[[index_name]] %in% label_id_vec)

    for (class in get_classes(n_class)) {
       docs[id_vec, class] <- 0
    }

    for (id in id_vec) {
      if (docs[id, labels_name] == 0) {
        docs[id, "Class_1"] <- 1
      } else {
        docs[id, "Class_2"] <- 1
      }
    }

  }

  return(docs)

}


#' @title Query Eta
#' @description Actively query eta parameter.
#' @param eta Eta parameter from EM.
#' @param class Class under consideration.
#' @param n_query Number of words to query.
#' @param handlabel Boolean indicating whether or not to label by hand.
#' @param true_eta If handlabel is false and an eta matrix is provided
#' to the true_eta param, values from the true_eta are used to determine
#' whether or not a word should be flagged.
#' @param keyword_select_scheme Keyword selection scheme. Either "ranked"
#' for top eta keyword selection, "ratio" for top eta ratio keyword selection,
#' or "combined" for log word count + log eta ratio selection.
#' @param val_scheme If handlabel is TRUE and true_eta is provided, sets
#' automated active updating scheme. Can be "ratio", in which case it
#' decides based on the ratio of the true etas, or "ranked", in which
#' case it decides based on the raw values of the true eta.
#' @param verbose Boolean deciding whether to print for debug.
#' @param dfm Quanteda document feature matrix.
#' @param keywords_list Existing list of keywords.
#' @return A vector of words.
query_eta <- function(eta, class, n_query, handlabel = TRUE,
                      true_eta = NA,
                      keyword_select_scheme = "ratio",
                      val_scheme = "ratio",
                      verbose = FALSE, dfm,
                      keywords_list = NA) {

  ## Get the column of eta into class columns, rank in descending order
  ## and subset.
  other_class <- ifelse(class == 1, 2, 1)
  if (keyword_select_scheme == "ranked") {
    prob_vec <- eta[, class]
  } else if (keyword_select_scheme == "ratio") {
    prob_vec <- eta[, class] - eta[, other_class]
  } else if (keyword_select_scheme == "combined") {
    log_term_freq <- log(colSums(as.matrix(dfm)))
    prob_vec <- (eta[, class] - eta[, other_class]) + log_term_freq
  }

  words_vec_ordered <- names(prob_vec[order(-prob_vec)])

  ## if a list of keywords is provided, remove them from consideration
  if (length(keywords_list) > 1 && !is.na(keywords_list[[class]])) {
    words_vec_ordered <- words_vec_ordered[
      which(!(words_vec_ordered %in% keywords_list[[class]]))
    ]
  }

  words_vec_subset <- words_vec_ordered[1:n_query]

  if (handlabel) {
    ##  Query the user re: the appropriateness of the words.
    wrong_words <- select.list(
      choices = c(words_vec_subset, NA),
      multiple = TRUE,
      title = paste0(
        "The model currently believes that these words are highly representative of Class ",
        class, ". Please indicate if any of these words are incorrectly associated with this class.",
        " If all words are correctly associated, choose 'NA'."
      )
    )

    right_ind <- which(!(words_vec_subset %in% wrong_words))
    right_words <- `if`(
      length(right_ind) != 0, words_vec_subset[right_ind], NA
    )

  } else if (length(true_eta) != 1 && !handlabel) {
    ## Alternatively, use the true eta parameter to decide automatically.

    if (val_scheme == "ratio") {
      true_ratio <- true_eta[, class] - true_eta[, other_class]
      neg_ratio <- true_ratio[true_ratio <= 0]
      pos_ratio <- true_ratio[true_ratio > 0]
      top_neg_ratio <- names(neg_ratio[neg_ratio <= quantile(neg_ratio, 0.1)])
      top_pos_ratio <- names(pos_ratio[pos_ratio >= quantile(pos_ratio, 0.9)])
      wrong_ind <- which(words_vec_subset %in% top_neg_ratio)
      right_ind <- which(words_vec_subset %in% top_pos_ratio)
    } else if (val_scheme == "ranked") {

      true_prob_vec <- true_eta[, class]
      true_words_vec_subset <- names(
        true_prob_vec[order(-true_prob_vec)]
      )[1:n_query]

      true_other_vec <- true_eta[, other_class]
      true_other_vec_ordered <- names(
        true_prob_vec[order(-true_other_vec)]
      )[1:n_query]

      wrong_ind <- which(
        (words_vec_subset %in% true_other_vec_ordered) &
        !(words_vec_subset %in% true_words_vec_subset)
      )

      right_ind <- which(words_vec_subset %in% true_words_vec_subset)
    }

    wrong_words <- `if`(
      length(wrong_ind) != 0, words_vec_subset[wrong_ind], NA
    )
    right_words <- `if`(
      length(right_ind) != 0, words_vec_subset[right_ind], NA
    )

    if (verbose) {
      cat("Class: ", class, "\n")
      cat("wrong words: ", wrong_words, "\n")
      cat("right words: ", right_words, "\n\n")
    }

  } else {
    stop("handlabel variable must be TRUE if true_eta is not provided.")
  }

  ## Export the information about the scaled keywords
  return(list(wrong_words = wrong_words, right_words = right_words))
}

active_eta_update <- function(beta_tbl, prev_model_output,
                              n_class, n_query, gamma,
                              handlabel, true_eta, dfm,
                              keywords_list = NA,
                              update_scheme = c("update_list", "inc_gamma"),
                              verbose = TRUE) {
  for (class in 1:n_class) {
    other_class <- ifelse(class == 1, 2, 1)
    new_keywords <- query_eta(
      eta = prev_model_output[[length(prev_model_output)]]$eta,
      class = class, n_query, handlabel, true_eta, dfm = dfm,
      keywords_list = keywords_list
    )

    if (!is.na(new_keywords$right_words)) {
      if (update_scheme == "inc_gamma") {
        right_keyword_idx <- which(
          rownames(beta_tbl) %in% new_keywords$right_words
        )
        beta_tbl[right_keyword_idx, class] <-
          beta_tbl[right_keyword_idx, class] + gamma
      } else if (update_scheme == "update_list") {
        if (verbose) {
          cat("\nAdding", new_keywords$right_words, "to keywords for Class", class)
        }
        if (is.na(keywords_list[[class]][1])) {
          keywords_list[[class]] <- new_keywords$right_words
        } else {
          keywords_list[[class]] <- c(
            keywords_list[[class]], new_keywords$right_words
          )
        }

      }
    }

    if (!is.na(new_keywords$wrong_words)) {
      if (update_scheme == "inc_gamma") {
        wrong_keyword_idx <- which(
          rownames(beta_tbl) %in% new_keywords$wrong_words
        )
        beta_tbl[wrong_keyword_idx, other_class] <-
          beta_tbl[wrong_keyword_idx, other_class] + gamma
      }
    }

  }

  if (update_scheme == "update_list") {
    ## reinitialize beta table with updated keywords
    beta_tbl <- initialize_beta_tbl(list(dfm), n_class, keywords_list, gamma)
  }

  return(list(beta_tbl = beta_tbl, keywords_list = keywords_list))
}

convert_beta_tbl <- function(beta_tbl) {
  names <- beta_tbl[[1]]
  beta_mtx <- as.matrix(beta_tbl[, -1])
  rownames(beta_mtx) <- names
  return(beta_mtx)
}

#' @export
get_index <- function(docs, index_name) {

    #' @title Get Index
    #' @description Gets index from a subset of doucments.
    #'
    #' @param docs       [numeric]   Documents chosen to initialize system.
    #' @param index_name  [character] Character string indicating the variable in 'docs'
    #'                               that denotes the index value of the document to be classified.
    #'
    #' @return          [vector] Vector of classes.

    index <- docs %>%
        dplyr::pull(!!dplyr::sym(index_name))

    return(index)

}

#' @export
get_clusters <- function(n_cluster) {

    #' @title Get Clusters
    #' @description Builds a vector of classification options.
    #'
    #' @param n_class   [numeric] Number of clusters to be considered.
    #'
    #' @return          [vector] Vector of clusters.

    clusters <- c()
    for (j in 1:n_cluster) {
        clusters <- c(clusters, paste("Cluster", j, sep="_"))
    }

    return(clusters)
}

#' @export
get_classes <- function(n_class) {

    #' @title Get Classes
    #' @description Builds a vector of classification options.
    #'
    #' @param n_class   [numeric] Number of classes to be considered.
    #'
    #' @return          [vector] Vector of classes.

    classes <- c()
    for (j in 1:n_class) {
        classes <- c(classes, paste("Class", j, sep="_"))
    }

    return(classes)
}


#' @export
get_dfm <- function(docs, doc_name = "text", index_name = "id", stem=T, ngrams=1, trimPct=0.0001, min_doc_freq=2, idfWeight=F, removeStopWords=T, minChar=4) {

    #' @title Get Document Feature Matrix
    #' @description Builds document feature matrix using quanteda package.
    #'
    #' @param docs       [matrix]    Matrix of labeled and unlabeled documents.
    #' @param doc_name    [character] Character string indicating the variable in 'docs'
    #'                               that denotes the text of the documents to be classified.
    #' @param index_name  [character] Character string indicating the variable in 'docs'
    #'                               that denotes the index value of the document to be classified.
    #' @param stem       [logical]   Switch indicating whether or not to stem terms.
    #' @param ngrams     [integer]   Integer value indicating the size of the ngram to use to build the dfm.
    #' @param trimPct    [numeric]   Numeric value indicating the threshold of percentage of document
    #'                               membership at which to remove terms from the data-term matrix.
    #'                               E.g., if \code{trimPct = .5}, then all words that are in less than
    #'                               50 percent of the documents will be removed.
    #' @param min_doc_freq [integer]   Minimum number of documents a term must be in to stay in the document term matrix.
    #' @param idfWeight  [logical]   Switch indicating whether to weight the document term matrix by the frequency of
    #'                               word counts. Only works if \code{dfmType = "quanteda"}.
    #'
    #' @return          [matrix]    Document term matrix.

    # If ngrams > 1, this if condition handles removal of stopwords properly.
    if (ngrams == 1) {
        dfm <- docs %>%
            quanteda::corpus(docid_field=index_name, text_field=doc_name) %>%
            quanteda::dfm(tolower=T, remove_numbers=T, remove_url=T, remove_twitter=T,
                          remove_punct=T, remove_hyphens=T, stem=stem, ngrams=ngrams) %>%
            {if (removeStopWords) quanteda::dfm_remove(., quanteda::stopwords(source="stopwords-iso")) else .} %>%
            quanteda::dfm_select(min_nchar=minChar) %>%
            quanteda::dfm_trim(min_termfreq=trimPct, min_docfreq=min_doc_freq,
                               termfreq_type="prop") %>%
            {if (idfWeight) quanteda::dfm_tfidf(.) else .}
    } else {
        dfm <- docs %>%
            quanteda::corpus(docid_field=index_name, text_field=doc_name) %>%
            quanteda::tokens(remove_numbers=T, remove_url=T, remove_twitter=T,
                             remove_punct=T, remove_hyphens=T) %>%
            quanteda::tokens_tolower() %>%
            {if (removeStopWords) quanteda::tokens_remove(., quanteda::stopwords(source="stopwords-iso")) else .} %>%
            {if (stem) quanteda::tokens_wordstem(.) else .} %>%
            quanteda::tokens_select(min_nchar=minChar) %>%
            quanteda::tokens_ngrams(n=ngrams) %>%
            quanteda::dfm() %>%
            quanteda::dfm_trim(min_termfreq=trimPct, min_docfreq=min_doc_freq, termfreq_type="prop")  %>%
            {if (idfWeight) quanteda::dfm_tfidf(.) else .}
    }

   return(dfm)

}

#' @export
split_dfm <- function(dfm, splitIndex) {

    #' @title Splits Document Term Matrix into two parts.
    #' @description Splits a document term matrix according to the provided \code{splitIndex}. Outputs a list of
    #'              two sections: the first that belongs to the \code{splitIndex}, and the second that does not.
    #'
    #' @param dfm        [matrix]    Output of \code{get_dfm}. Document-Term Matrix.
    #' @param splitIndex [vector]    Vector of strings denoting the row values to split on.
    #'
    #' @return           [list]      List of two parts of the document term matrix.

    first_section <- dfm[rownames(dfm) %in% splitIndex, ]
    second_section <- dfm[!(rownames(dfm)) %in% splitIndex, ]

    return(list(first_section=first_section, second_section=second_section))

}

get_class_matrix <- function(docs, n_class, hand_labeled_index,
                             doc_name, index_name, NB_init) {

#' @title        Get Document Class Probability Matrix
#' @description  Gets matrix of document classifications.
#'
#' @param docs               [matrix]    Matrix of labeled and unlabeled documents.
#' @param n_class            [numeric]   Number of classes to be considered.
#' @param hand_labeled_index   [vector]    Vector of index values for hand labeled documents in \code{docs}.
#' @param doc_name            [character] Character string indicating the variable in 'docs'
#'                                       that denotes the text of the documents .
#' @param index_name          [character] Character string indicating the variable in 'docs'
#'                                       that denotes the index value of the documents .
#'
#' @return                   [matrix]    Class assignment matrix.

  C_train <- docs %>%
    dplyr::filter(!!dplyr::sym(index_name) %in% hand_labeled_index) %>%
    dplyr::select_at(dplyr::vars(index_name, dplyr::matches("^Class")))
  index <- C_train %>% pull(!!dplyr::sym(index_name))
  col_names <- colnames(C_train)
  C_train <- C_train %>% dplyr::select(-!!dplyr::sym(index_name))
  C_train <- purrr::map(C_train, Matrix::Matrix, sparse = T) %>%
    purrr::reduce(cbind2)
  rownames(C_train) <- index

  return(C_train)

}

get_entropy <- function(data) {

    #' @title         Get Entropy Value
    #' @description   Gets row-wise entropy values from a rectangular data array.
    #'
    #' @param   [matrix] Matrix of data.
    #'
    #' @return  [vector] Vector of entropy values.

    data <- as.matrix(data)
    entropy <- rep(0, nrow(data))

    for (i in 1:nrow(data)) {
        for (j in 1:ncol(data)) {
            entropy[i] <- entropy[i] + -(data[i, j] * log(data[i, j]))
        }
    }

    entropy[is.nan(entropy)] <- 0
    # add small sum to entropy values, take log
    entropy <- log(entropy + 0.000000000000000001)
    return(entropy)

}

region_sample_edge <- function(docs, max_query, edge = T, regions = "both") {


  #' @title Log-Ratio Region Sampling
  #'
  #' @param docs [dataframe] Documents with log ratio and cumulative sum columns.
  #' @param max_query [integer] Number of documents to be queried.
  #' @param edge [logical] Whether we sapmle from the edge of the posterior distributions as well
  #'
  #' @return Sampled documents.
  #'
  #'

  if (edge == T){
    # If we sample from the edge of the posterior distributions,
    # (i.e. sample documents that has almost 0 or 1 posterior to avoid
    # group of documents from being classified wrongly)
    # we sample 1/2 (default 5 docs) from the uncertain region in the log-likelihood ratio,
    # and 1/2 from the edge cases.
    fraction_edge <- 1/2
    fraction <- fraction_edge / 4
  }else{
    fraction <- 1/4
  }

  neg_sample <- which(docs$cum_sum_neg == 0)[1]
  neg_sample_range <- seq(neg_sample - floor(fraction * max_query),
                          neg_sample + floor(fraction * max_query) - 1)
  pos_sample <- which(docs$cum_sum_pos == 0)[length(which(docs$cum_sum_pos == 0))]
  pos_sample_range <- seq(pos_sample - floor(fraction * max_query),
                          pos_sample + floor(fraction * max_query) - 1)
  # combine sampling regions, and pull associated index values
  comb_sample_range <- c(neg_sample_range, pos_sample_range)


  # removing this for the moment

  # if (length(unique(comb_sample_range)) != max_query) {
  #   comb_sample_range <- unique(comb_sample_range)
  #   diff <- max_query - length(comb_sample_range)
  #   if (diff %% 2 == 0) {
  #     append_neg <- seq(min(comb_sample_range ) - diff / 2,
  #                       min(comb_sample_range) - 1)
  #     append_pos <- seq(max(comb_sample_range) + 1,
  #                       max(comb_sample_range) + diff / 2)
  #     comb_sample_range <- c(append_neg,
  #                            comb_sample_range,
  #                            append_pos)
  #   } else {
  #     paired_vals <- comb_sample_range[-1]
  #     append_neg <- seq(min(paired_vals) - diff / 2,
  #                       min(paired_vals) - 1)
  #     append_pos <- seq(max(paired_vals) + 1,
  #                       max(paired_vals) + diff / 2)
  #     comb_sample_range <- c(append_neg,
  #                            comb_sample_range,
  #                            append_pos)
  #
  #     coin <- rbinom(1, 1, 0.5)
  #     if (coin == 1) {
  #       comb_sample_range <- c(comb_sample_range,
  #                              max(comb_sample_range) + 1)
  #     } else {
  #       comb_sample_range <- c(comb_sample_range,
  #                              min(comb_sample_range) - 1)
  #
  #     }
  #   }
  # }

  if (edge == T){
    # add the samples from the edge

    # for now, we sample some documents from the 100 documents that have very high/low log likelihood
    pos_edge <- sample(1:100, floor(fraction_edge/2 * max_query))
    neg_edge <- sample(seq(nrow(docs)-100, nrow(docs)), floor(fraction_edge/2 * max_query))
    comb_sample_range <- c(pos_edge, comb_sample_range, neg_edge)
  }

  # sampling the appropriate documents
  if (regions == "both") {
          docs <- docs[comb_sample_range, ]
  } else if (regions == "pos") {
      if (edge == T) {
          to_sample_pos <- c(pos_edge, pos_sample_range)
          docs <- docs[to_sample_pos, ]
      } else {
          docs <- docs[pos_sample_range, ]
      }
  } else if (regions == "neg") {
      if (edge == T) {
          to_sample_neg <- c(neg_edge, neg_sample_range)
          docs <- docs[to_sample_pos, ]
      } else {
          docs <- docs[neg_sample_range, ]
      }
  }

  return(docs)
}


region_sample <- function(docs, max_query) {

    #' @title Log-Ratio Region Sampling
    #'
    #' @param docs [dataframe] Documents with log ratio and cumulative sum columns.
    #' @param max_query [integer] Number of documents to be queried.
    #'
    #' @return Sampled documents.

    neg_sample <- which(docs$cum_sum_neg == 0)[1]
    neg_sample_range <- seq(neg_sample - floor(1 / 4 * max_query),
                            neg_sample + floor(1 / 4 * max_query) - 1)
    pos_sample <- which(docs$cum_sum_pos == 0)[length(which(docs$cum_sum_pos == 0))]
    pos_sample_range <- seq(pos_sample - floor(1 / 4 * max_query),
                            pos_sample + floor(1 / 4 * max_query) - 1)
    # combine sampling regions, and pull associated index values
    comb_sample_range <- c(neg_sample_range, pos_sample_range)

    if (length(unique(comb_sample_range)) != max_query) {
        comb_sample_range <- unique(comb_sample_range)
        diff <- max_query - length(comb_sample_range)
        if (diff %% 2 == 0) {
            append_neg <- seq(min(comb_sample_range ) - diff / 2,
                              min(comb_sample_range) - 1)
            append_pos <- seq(max(comb_sample_range) + 1,
                              max(comb_sample_range) + diff / 2)
            comb_sample_range <- c(append_neg,
                                   comb_sample_range,
                                   append_pos)
        } else {
            paired_vals <- comb_sample_range[-1]
            append_neg <- seq(min(paired_vals) - diff / 2,
                              min(paired_vals) - 1)
            append_pos <- seq(max(paired_vals) + 1,
                              max(paired_vals) + diff / 2)
            comb_sample_range <- c(append_neg,
                                   comb_sample_range,
                                   append_pos)

            coin <- rbinom(1, 1, 0.5)
            if (coin == 1) {
                comb_sample_range <- c(comb_sample_range,
                                       max(comb_sample_range) + 1)
            } else {
                comb_sample_range <- c(comb_sample_range,
                                       min(comb_sample_range) - 1)

            }
        }
    }

    docs <- docs[comb_sample_range, ]

    return(docs)
}

log_ratio_sample <- function(docs, out, dfm, mu, tau, max_query, edge, regions) {

    #' @title Log-Ratio Sampling
    #'
    #' @param docs [dataframe] Documents with log ratio and cumulative sum columns.
    #' @param out [list] List of output objects from `EM()` function.
    #' @param dfm [quanteda dfm] Quanteda dfm object.
    #' @param mu Error acceptance rate for first region.
    #' @param eta Error acceptance rate for second region.
    #' @param max_query [integer] Number of documents to be queried.
    #' @param edge [logical] Whether we sample from documents that have very high/low posterior
    #'
    #' @return Sampled documents.

    # getting the log ratio per document
    eta_ratio <- out$eta[, "Class_2"] -
        out$eta[, "Class_1"]
    pos_etas <- out$eta[, "Class_2"]
    neg_etas <- out$eta[, "Class_1"]

    docs <- docs %>%
        dplyr::mutate(log_ratio = as.vector(dfm %*% eta_ratio),
                      pr_d_eta_pos = exp(as.vector(dfm %*% pos_etas)),
                      pr_d_eta_neg = exp(as.vector(dfm %*% neg_etas))) %>%
        dplyr::arrange(desc(log_ratio)) %>%
        dplyr::mutate(cum_sum_neg = as.numeric(cumsum(pr_d_eta_pos) /
                                                   sum(pr_d_eta_pos) < mu),
                      cum_sum_pos = as.numeric(rev(cumsum(rev(pr_d_eta_neg)) /
                                                       sum(pr_d_eta_neg)) < tau)) %>%
        region_sample_edge(max_query, edge = edge, regions = regions)

    return(docs)

}

#' @export
get_uncertain_docs <- function(docs, bound, max_query,
                               index_name, hand_labeled_index, force_list=F,
                               query_type="basic_entropy",
                               quantileBreaks=c(75, 20),
                               sampleProps=c(.5, .3, .2),
                               mu=0.001,
                               tau=0.001,
                               regions="both",
                               dfm = NULL,
                               ## EM_out,
                               seed=NULL,
                               n_cluster = NULL) {

#' @title        Get Uncertain Documents
#' @description  Get documents that the previous iteration of the EM algorithm is least sure about.
#'
#' @param docs               [matrix]    Matrix of labeled and unlabeled documents.
#' @param bound              [numeric]   The choice of lower bound for entropy-based uncertainty selection.
#' @param max_query           [numeric]   Maxmium number of uncertain documents that can be queried.
#' @param index_name          [character] Character string indicating the variable in 'docs'
#'                                       that denotes the index value of the documents .
#' @param hand_labeled_index   [vector]    Vector of index values for hand labeled documents in \code{docs}.
#' @param force_list          [logical]   Switch indicating whether to force the filtering of documents with
#'                                       no entropy. Set to \code{FALSE} by default.
#' @param query_type          [string]    String indicating which type of uncertainty sampling to use. Options are \code{"standard_entropy"},
#'                                       \code{"normalized_entropy"}, \code{"tiered_entropy"}, or
#'                                       \code{"tiered_entropy_weighted"}.
#' @param quantileBreaks     [vector]    Vector of break points to distinguish entropy zones. The first value is
#'                                       the break point between the first and second tier, the second is the
#'                                       break point between the second and third tier.
#' @param sampleProps        [vector]    Vector of sampling proportions for each entropy zone. The first value is
#'                                       the proportion of \code{max_query} to be sampled from the high entropy region,
#'                                       the second value is the proportion to be sampled from the middle entropy region,
#'                                       and the third value is the proportion to be sampled from the lowest entropy region.
#' @param n_cluster          [int]       Number of clusters.
#'
#' @return                   [vector]    Vector of id values of documents that the EM algorithm is uncertain about.

  error_sample <- function(docs) {
    warning("Insufficient entropy, sampling randomly.")
    ## warning("Insufficient entropy, breaking")
    ## break
    return(dplyr::sample_n(docs, max_query))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }


  ## Calculates entropy across clusters, rather than classes, and uses
  ## this quantity for entropy sampling.
  if (query_type == "basic_entropy_cluster") {
    uncertainClass <- docs %>%
      dplyr::filter(!(!!dplyr::sym(index_name)) %in% hand_labeled_index) %>%
      {if (force_list == T) {
       {tryCatch(
          dplyr::filter_at(., dplyr::vars(dplyr::matches("^Cluster")),
                           dplyr::all_vars(get_entropy(.) >= bound)),
          error=function(e) error_sample(.)
        )}
       } else {
         dplyr::filter_at(., dplyr::vars(dplyr::matches("^Cluster")),
                          dplyr::all_vars(get_entropy(.) > bound))
       }} %>%
      dplyr::arrange_at(dplyr::vars(dplyr::matches("^Cluster")),
                        dplyr::all_vars(desc(get_entropy(.)))) %>%
      dplyr::slice(1:max_query)
  }

  ## Get the difference between the positive cluster and the most
  ## negative cluster, and actively label the documents for which
  ## the difference is the lowest.
  if (query_type == "margin_cluster") {
    options(dplyr.show_progress = FALSE)
    pos_cluster <- paste0("Cluster_", n_cluster)
    max_neg_cluster_val <- docs %>%
      dplyr::select(-pos_cluster) %>%
      dplyr::select(matches("^Cluster")) %>%
      dplyr::rowwise() %>%
      do( (.) %>%
          as.data.frame() %>%
          dplyr::mutate(max_neg_cluster_val = max(.))) %>%
      dplyr::pull(max_neg_cluster_val)
    pos_cluster_val <- docs[[pos_cluster]]
    uncertainClass <- docs %>%
      dplyr::mutate(marg_diff = abs((max_neg_cluster_val -
                              pos_cluster_val))) %>%
      dplyr::arrange(marg_diff) %>%
      dplyr::slice(1:max_query)
  }

  if (query_type == "basic_entropy") {
    entropy <- docs %>%
      select_at(vars(matches("^Class"))) %>%
      get_entropy()

    uncertainClass <- docs %>%
      mutate(entropy = entropy) %>%
      dplyr::filter(!(!!dplyr::sym(index_name)) %in% hand_labeled_index) %>%
      arrange(desc(entropy)) %>%
      slice(1:max_query) %>%
      select(-entropy)
  }

  if (query_type == "normalized_entropy") {
    uncertainClass <- docs %>%
      dplyr::filter(!(!!dplyr::sym(index_name)) %in% hand_labeled_index) %>%
      dplyr::mutate(entropy = dplyr::select_at(., dplyr::vars(dplyr::matches("^Class"))) %>%
                      get_entropy()) %>%
      {tryCatch(dplyr::sample_n(., size=max_query, weight=entropy),
                error=function(e) error_sample(.))} %>%
      dplyr::select(-entropy)
  }

  if (query_type == "tiered_entropy") {
    uncertainClass <- docs %>%
      dplyr::filter(!(!!dplyr::sym(index_name)) %in% hand_labeled_index) %>%
      dplyr::mutate(entropy = dplyr::select_at(., dplyr::vars(dplyr::matches("^Class"))) %>%
                      get_entropy(),
                    entropy_rank = ntile(entropy, n=100)) %>%
      {tryCatch(bind_rows(
         uncertainClass %>% dplyr::filter(entropy_rank > quantileBreaks[1]) %>%
         dplyr::sample_n(max_query * sampleProps[1]),
         uncertainClass %>% dplyr::filter(entropy_rank > quantileBreaks[2] & entropy_rank < quantileBreaks[1] + 1) %>%
         dplyr::sample_n(max_query * sampleProps[2]),
         uncertainClass %>% dplyr::filter(entropy_rank < quantileBreaks[2] + 1) %>%
         dplyr::sample_n(max_query * sampleProps[3])),
         error=function(e) error_sample(.))} %>%
      dplyr::select(-entropy, -entropy_rank)

  }

  if (query_type == "tiered_entropy_weighted") {
    uncertainClass <- docs %>%
      dplyr::filter(!(!!dplyr::sym(index_name)) %in% hand_labeled_index) %>%
      dplyr::mutate(entropy = dplyr::select_at(., dplyr::vars(dplyr::matches("^Class"))) %>%
                      get_entropy(),
                    entropy_rank = ntile(entropy, n=100)) %>%
      {tryCatch(dplyr::bind_rows(
                         uncertainClass %>% dplyr::filter(entropy_rank > quantileBreaks[1]) %>%
                         dplyr::sample_n(max_query * sampleProps[1], weight=entropy),
                         uncertainClass %>% dplyr::filter(entropy_rank > quantileBreaks[2] & entropy_rank < quantileBreaks[1] + 1) %>%
                         dplyr::sample_n(max_query * sampleProps[2], weight=entropy),
                         uncertainClass %>% dplyr::filter(entropy_rank < quantileBreaks[2] + 1) %>%
                         dplyr::sample_n(max_query * sampleProps[3], weight=entropy)),
                error=function(e) error_sample(.))} %>%
      dplyr::select(-entropy, -entropy_rank)

  }

  if (query_type == "random") {
    uncertainClass <- docs %>%
      dplyr::filter(!(!!dplyr::sym(index_name)) %in% hand_labeled_index) %>%
      dplyr::sample_n(max_query)
  }

  if (query_type == "log_ratio") {
    print("Log ratio currenting not working")
    break
    ## uncertainClass <- log_ratio_sample(docs=docs, out=EM_out,
                                       ## dfm=dfm, mu=mu, tau=tau,
                                       ## max_query=max_query,
                                       ## edge = F,
                                       ## regions = regions)
  }

  if (query_type == "log_ratio_edge"){
    print("Log ratio currenting not working")
    break
    ## uncertainClass <- log_ratio_sample(docs=docs, out=EM_out,
    ##                                    dfm=dfm, mu=mu, tau=tau,
    ##                                    max_query=max_query,
    ##                                    edge = T,
    ##                                    regions = regions)
  }

  return(uncertainClass[[index_name]])

}

#' @export
matchCluster2Class <- function(output, count, n_cluster, n_class) {
  pos_class <- output[, n_cluster]
  neg_class <- 1 - pos_class
  obj <- output
  obj <- obj[, -(1:(n_cluster - n_class))]
  obj[, 1] <- neg_class
  obj[, 2] <- pos_class
  return(obj)
}

#' @export
match_clusters_to_docs <- function(docs, EMoutput, index_name, n_cluster) {
#' @title Match Multicluster EM Output to Document Matrix
#' @description Matches the output of multicluster EM to the document corpus matrix.
#' @return Matrix of documents.

  colnames(EMoutput) <- get_clusters(n_cluster)

  row_names <- rownames(EMoutput)
  EM_out_tbl <- as.matrix(EMoutput) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(!!dplyr::sym(index_name) := row_names)
  docs <- docs %>%
    left_join(EM_out_tbl, by = index_name)
  return(docs)

}

match_EM_to_docs <- function(docs, EMoutput, classes, doc_name, index_name,
                             labels_name=NULL) {

#' @title        Match EM Output to Document Matrix
#' @description  Matches the output of the previous run of the EM algorithm to the matrix of documents.
#'
#' @param docs        [matrix]    Matrix of labeled and unlabeled documents.
#' @param EMoutput    [matrix]    Matrix of expected class assignments produced by EM algorithm..
#' @param classes     [vector]    Vector of character strings indicating the classes being considered.
#' @param doc_name     [character]
#' @param index_name   [character] Character string indicating the variable in 'docs'
#'                                that denotes the index value of the documents .
#' @param labels_name  [character] Character string indicating the variable in \code{docs}
#'                                that denotes the already known labels of the documents.
#'                                By default, value is set to \code{NULL}.
#'
#' @return            [matrix]    Matrix of documents that the EM algorithm is uncertain about.
  colnames(EMoutput) <- classes


  match_type <- class(docs[[paste0(index_name)]])
  to_join <- dplyr::as_tibble(as.matrix(EMoutput),
                              rownames=paste0(index_name)) %>%
    dplyr::mutate(!!dplyr::sym(index_name) :=
                    if (match_type == "numeric" |
                        match_type == "integer") {
                      as.numeric(!!dplyr::sym(index_name))
                    } else {
                      as.character(!!dplyr::sym(index_name))
                    }
                  )

  if (is.null(labels_name) == T) {
    docs <- dplyr::left_join(
                     docs %>%
                     dplyr::select(doc_name, index_name),
                     to_join, by=paste0(index_name)
                   )
  } else {
    docs <- dplyr::left_join(
                     docs %>%
                     dplyr::select(doc_name, index_name, labels_name),
                     to_join, by=paste0(index_name)
                   )
  }

  return(docs)

}

#' @export
get_term_sparsity <- function(dfm) {
    freq_doc <- quanteda::docfreq(dfm)
    n <- quanteda::ndoc(dfm)

    return(freq_doc / n)
}
active_initial_messages <- function(n_cluster, query_type) {
  #' Prints initial messages for active sampling, if needed.
  if (n_cluster == 2 & query_type %in% c("margin_cluster", "basic_entropy_cluster")) {
    query_type <- "basic_entropy"
    message("Cluster sampling only works with greater than two clusters.
Defaulting to basic_entropy sampling scheme.")
  }
}

get_initial_documents <- function(docs, init_index, index_name,
                                  init_size, whichOutTest) {
  #' gets initial documents for active learning algorithm
  if (is.null(init_index)) {
    if (!is.null(whichOutTest)) {
      to_label_ind <- docs %>%
        dplyr::filter(!(!!dplyr::sym(index_name)) %in% whichOutTest) %>%
        dplyr::pull(!!dplyr::sym(index_name)) %>%
        sample(init_size)
    } else {
      to_label_ind <- docs %>%
        dplyr::pull(!!dplyr::sym(index_name)) %>%
        sample(init_size)
    }
  } else {
    to_label_ind <- init_index
  }

  return(to_label_ind)
}

generate_lambda_vec <- function(lambda_decay, lambda, rate, iters,
                                max_active) {
  #' Creates a vector of lambdas depending on decay value
  if (lambda_decay) {
    decay <- function(rate, iters) {
      data <- c()
      for (i in 0:iters) {
        data[i + 1] <- 1 / (1 + rate * i)
      }
      return(data)
    }
    lambda_vec <- decay(rate = ld_rate, iters = max_active)
  } else {
    lambda_vec <- rep(lambda, max_active + 1)
  }
  return(lambda_vec)
}

tune_lambda_in_active <- function(docs, index_name, hand_labeled_index, n_cluster,
                                  tune_lambda_range, tune_lambda_prop_init,
                                  tune_lambda_parallel, tune_lambda_k, seed) {
  #' Tunes lambda value at each iteraction, if enabled
  tuning_docs <- docs %>%
    dplyr::filter(!!dplyr::sym(index_name)
                  %in% hand_labeled_index[[count]])

  tune_out <- tune_lambda(
    docs = tuning_docs,
    n_clusters = n_cluster,
    lambdas = tune_lambda_range,
    prop_init = tune_lambda_prop_init,
    parallel = tune_lambda_parallel,
    k = tune_lambda_k,
    seed = seed
  )
  return(tune_out)
}

## get_prev_model_params <- function(count, NB_init) {
##   #' Grab previous model parameters if necessary
##   if (count == 1) {
##     choose_NB_init <<- TRUE
##     prev_word_prob <<- prev_class_prob <<- NULL
##   } else if (count > 1 & NB_init == FALSE) {
##     choose_NB_init <<- FALSE
##     prev_word_prob <<- output[[count - 1]]$eta
##     prev_class_prob <<- output[[count - 1]]$pi
##   }
## }

## get_oos_pred <- function(dfm, output, count, n_cluster, n_class,
##                          out_docs_static, doc_name, index_name, labels_name) {
##   #' gets out of sample prediction
##     out_prediction <- E_step(
##       .D_test=dfm,
##       .class_prob=output[[length(output)]]$pi,
##       .word_prob=output[[length(output)]]$eta
##     )
##     if (n_cluster > 2) {
##       EM_out_classlik <- matchCluster2Class(
##         exp(out_prediction),
##         count, n_cluster, n_class
##       )
##     } else {
##       EM_out_classlik <- exp(out_prediction)
##     }
##     out_docs <- match_EM_to_docs(
##       out_docs_static,
##       EMoutput = EM_out_classlik,
##       classes, doc_name, index_name,
##       labels_name
##     )
##   return(out_docs)
## }

check_lr_convergence <- function(output, count, log_ratio_threshold, log_ratio_conv_type) {
  #' check for convergence if using log-ratio sampling
  current_maximand <- output[[count]]$maximands[length(output[[count]]$maximands)]
  last_maximand <- output[[count - 1]]$maximands[length(output[[count - 1]]$maximands)]
  total_diff <- abs(current_maximand - last_maximand)
  if (total_diff < log_ratio_threshold & log_ratio_conv_type == "maximand") {
    message(paste0("\n", "Total maximand change: ", total_diff))
    message(paste0("\n", "Log Ratio Sampling: Convergence Reached"))
    break
  } else {
    message(paste0("\n", "Total maximand change: ", total_diff))
  }
}

#' @export
agg_helper_convert <- function(model_preds,
                               n_cluster_collapse_type = "simple") {
#' @title Aggregation Helper
#' @description helps aggregation function by collapsing clusters to classes (binary only).

  cluster_names <- colnames(select(model_preds, -dfm_id))
  pos_cluster <- cluster_names[length(cluster_names)]
  neg_clusters <- cluster_names[1:length(cluster_names) - 1]

  if (n_cluster_collapse_type == "simple") {
    model_preds <- model_preds %>%
      mutate(Class_1 = 1 - model_preds[[pos_cluster]],
             Class_2 = model_preds[[pos_cluster]])
  } else if (n_cluster_collapse_type == "max_neg") {
    ## Get maximum value of negative clusters by row, then normalize
    model_preds <- model_preds %>%
      mutate(Class_1 = do.call(pmax, model_preds[neg_clusters]),
             Class_2 = model_preds[[pos_cluster]],
             Class_1 = Class_1 / (Class_1 + Class_2),
             Class_2 = Class_2 / (Class_1 + Class_2))
  }

  return(model_preds)
}

#' @export
get_mean_mpe <- function(mod, dfm, val_data, labels_name = "label", index_name = "id",
                    n_cluster_collapse_type = "simple", n_cluster) {
#' @title Get Mean Prediction Error Singular
#' @description gets mean prediction error
  cluster_names <- get_clusters(n_cluster)
  out_prediction <- E_step(
    .D_test = dfm,
    .class_prob = mod$pi,
    .word_prob = mod$eta
  ) %>% as.matrix %>%
  `colnames<-`(cluster_names) %>%
  as_tibble(rownames = "id")

  class_preds <- agg_helper_convert(out_prediction, n_cluster_collapse_type)

  mean_mpe <- val_data %>%
    dplyr::select(!!dplyr::sym(labels_name), !!dplyr::sym(index_name)) %>%
    dplyr::left_join(class_preds, by = index_name) %>%
    dplyr::mutate(mpe = abs(!!dplyr::sym(labels_name) - exp(Class_2))) %>%
    summarize(mean_mpe = mean(mpe)) %>%
    pull(mean_mpe)

  return(mean_mpe)
}

#' @export
get_mean_mpes <- function(dfms, models, val_data, n_cluster) {
#' @title Get Mean Prediction Error
#' @description gets mean mpes across a list of dfms and list of models of equal length
  N <- length(dfms)
  mean_mpes <- c()
  for (i in 1:N) {
    mean_mpes[i] <- get_mean_mpe(models[[i]], dfms[[i]], val_data,
                                 n_cluster = n_cluster)
  }
  return(mean_mpes)
}

#' @export
get_alpha_m <- function(mean_mpe) {
#' @title Get alpha_m
#' @description gets normalized model weight
  alpha_m <- log((1 - mean_mpe) / mean_mpe)
  return(alpha_m)
}

#' @export
get_model_weights <- function(dfms, models, val_data, n_cluster) {
#' @title Get Model Weights
#' @description calculates the weights that each model recieves and normalize
  mean_mpes <- get_mean_mpes(dfms, models, val_data, n_cluster)
  model_weights <- sapply(mean_mpes, get_alpha_m)
  model_weights <- model_weights / sum(model_weights)
  return(model_weights)
}

#' @export
get_weighted_prediction <- function(model_preds, model_weights,
                                    index_name = "id",
                                    labels_name = "label") {
#' @title Get Weighted Prediction
#' @description gets weighted predictions for unlabled documents
  model_preds_w <- tibble(dfm_id = 1:length(model_weights),
                          model_weights = model_weights) %>%
    right_join(model_preds, by = "dfm_id") %>%
    group_by(!!sym(index_name), !!sym(labels_name)) %>%
    dplyr::summarize_at(dplyr::vars(dplyr::matches("^Cluster|^Class")),
                        weighted.mean, w = model_weights) %>%
    ungroup()
  target <- model_preds %>% filter(dfm_id == 1) %>% pull(!!sym(index_name))
  model_preds_w <- model_preds_w[match(target, model_preds_w[[index_name]]),]
  return(model_preds_w)
}

#' @export
choose_best_model <- function(model_preds, model_weights, index_name = "id") {
#' @title Choose best model.
#' @description gets predictions from single best model
  best_model_weight <- max(model_weights)[1]
  best_model_preds <- tibble(dfm_id = 1:length(model_weights),
                             model_weights = model_weights) %>%
    right_join(model_preds, by = "dfm_id")
  best_dfm_id <- best_model_preds %>%
    filter(model_weights == best_model_weight) %>%
    pull(dfm_id) %>%
    unique %>%
    sample(1)
  best_model_preds <- best_model_preds %>%
    filter(dfm_id == best_dfm_id)
  return(best_model_preds)
}


#' @export
aggregate_model_predictions <- function(pred_lst,
                                        dfms = NULL, models = NULL,
                                        val_data = NULL, n_cluster,
                                        agg_type = "random",
                                        n_cluster_collapse_type = "simple") {

#' @title Aggregate Model Predictions
#' @description Processes model predictions according to cluster structure,
#' then chooses between model predictions for each dfm being
#' used to fit a model.

  pred_tbl_in <- agg_helper_convert(
    pred_lst$model_output_in, n_cluster_collapse_type
  )

  if (!is.null(pred_lst$model_output_out)) {
    pred_tbl_out <- agg_helper_convert(
      pred_lst$model_output_out, n_cluster_collapse_type
    )
  }

  if (max(pred_tbl_in$dfm_id) == 1) {

    in_agg <- select(pred_tbl_in, -dfm_id)
    if (!is.null(pred_lst$model_output_out)) {
      out_agg <- select(pred_tbl_out, -dfm_id)
    }

  } else if (agg_type == "random") {

    rdm_dfm_id <- sample(1:max(unique(pred_tbl_in$dfm_id)), 1)
    in_agg <- pred_tbl_in %>%
      filter(dfm_id == rdm_dfm_id)
    if (!is.null(pred_lst$model_output_out)) {
      out_agg <- pred_tbl_out %>%
        filter(dfm_id == rdm_dfm_id)
    }

  } else if (agg_type == "best") {

    model_weights <- get_model_weights(dfms, models, val_data, n_cluster)
    in_agg <- choose_best_model(pred_tbl_in, model_weights)
    if (!is.null(pred_lst$model_output_out)) {
      out_agg <- choose_best_model(pred_tbl_out, model_weights)
    }

  } else if (agg_type == "weighted_avg") {

    model_weights <- get_model_weights(dfms, models, val_data, n_cluster)
    in_agg <- get_weighted_prediction(pred_tbl_in, model_weights)
    if (!is.null(pred_lst$model_output_out)) {
      out_agg <- get_weighted_prediction(pred_tbl_out, model_weights)
    }

  }

  out_lst <- list(
    in_agg = in_agg, out_agg = NULL,
    agg_type = agg_type,
    n_cluster_collapse_type = n_cluster_collapse_type
  )

  if (!is.null(pred_lst$model_output_out)) {
    out_lst[["out_agg"]] <- out_agg
  }

  return(out_lst)
}

#' @export
gen_results_tbl <- function(include_out_stats, metadata, max_iters, model_name) {
#' @title Generate Results Table
#' @description generates an object for storing model results

  res_obj <- tibble(
    model_name = rep(model_name, max_iters),
    iter = 0,
    docs_labeled = 0,
    accuracy_in = 0,
    precision_in = 0,
    recall_in = 0,
    F1_in = 0
  )
  if (include_out_stats) {
    res_obj <- bind_cols(
      res_obj,
      tibble(
        accuracy_out = rep(0, max_iters),
        precision_out = 0,
        recall_out = 0,
        F1_out = 0
      )
    )
  }
  res_obj <- bind_cols(
    res_obj,
    as_tibble(metadata, .rows = max_iters)
  )


  return(res_obj)
}

#' @export
update_results <- function(include_out_stats, res_obj,
                           agg_output_in, agg_output_out,
                           hl_index, i, time_sec) {
#' @title Update model results.
#' @description updates model results

  cf_in <- get_conf_matrix(
    docs = agg_output_in,
    labeledIndex = hl_index,
    useLabeled = FALSE
  )

  res_obj[i, ]$iter <- i - 1
  res_obj[i, ]$docs_labeled <- length(hl_index)
  res_obj[i, ]$accuracy_in <- get_classification_accuracy(cf_in)
  res_obj[i, ]$precision_in <- get_precision_binary(cf_in)
  res_obj[i, ]$recall_in <- get_recall_binary(cf_in)
  res_obj[i, ]$F1_in <- get_F1_binary(cf_in)

  if (include_out_stats) {
    cf_out <- get_conf_matrix(
      docs = agg_output_out,
      labeledIndex = hl_index,
      useLabeled = FALSE
    )
    res_obj[i, ]$accuracy_out <- get_classification_accuracy(cf_out)
    res_obj[i, ]$precision_out <- get_precision_binary(cf_out)
    res_obj[i, ]$recall_out <- get_recall_binary(cf_out)
    res_obj[i, ]$F1_out <- get_F1_binary(cf_out)
  }

  # make a time column if it does not exist
  if (all(colnames(res_obj) != "time_sec")){
    res_obj$time_sec <- NULL
  }
  res_obj$time_sec[i] <- time_sec
  return(res_obj)
}

initialize_beta_tbl <- function(dfms, n_class, keywords_list = NA, gamma) {
  #' initializes a table of prior values for eta

  if (length(dfms) != 1) {
    stop("Active eta tuning feature only works with singular DFM.")
  }

  beta_tbl <- matrix(nrow = ncol(dfms[[1]]), ncol = n_class, data = 2)
  rownames(beta_tbl) <- colnames(dfms[[1]])

  for (i in 1:n_class) {
    if (!is.na(keywords_list[[i]])) {
      key_class_idx <- which(rownames(beta_tbl) %in% keywords_list[[i]])
      beta_tbl[key_class_idx, i] <- beta_tbl[key_class_idx, i] + gamma
    }
  }

  return(beta_tbl)
}

update_em_param_tbl <- function(em_param_tbl, model_output, base_index, id) {
  model_predictions <- tibble() %>%
    mutate(dfm_id = id) %>%
    left_join(
      as_tibble(
        as.matrix(model_outputs[[i]]$classLik),
        rownames = paste0(index_name)
      ),
      by = paste0(index_name)
    ) %>%
    bind_rows(model_predictions_in)
}

model_out_to_tbl <- function(model_outputs) {
  model_output_in_lst <- list()

  if (length(model_outputs[[1]]$out_prediction) != 0) {
    model_output_out_lst <- list()
  }

  for (i in 1:length(model_outputs)) {
    model_output_in_lst[[i]] <- model_outputs[[i]]$classLik %>%
      as.matrix() %>%
      as_tibble(rownames = "id") %>%
      `colnames<-`(c("id", get_clusters(ncol(.)))) %>%
      mutate(dfm_id = i)

    if (length(model_outputs[[i]]$out_prediction) != 0) {
      model_output_out_lst[[i]] <- model_outputs[[i]]$out_prediction %>%
        as.matrix() %>%
        as_tibble(rownames = "id") %>%
        `colnames<-`(c("id", get_clusters(ncol(.)))) %>%
        mutate(dfm_id = i)
    }
  }

  out_lst <- list(model_output_in = bind_rows(model_output_in_lst))
  if (length(model_outputs[[1]]$out_prediction != 0)) {
    out_lst[["model_output_out"]] <- bind_rows(model_output_out_lst)
  }

  return(out_lst)
}

update_docs <- function(docs_old, new_data, classes = get_classes(n_class)) {
  if (!is.null(docs_old)) {
    cols_to_add <- new_data %>%
      dplyr::select(id, dplyr::matches("^Class"))
    base <- docs_old %>%
      dplyr::select(!dplyr::matches("^Class"))

    return(
      left_join(base, cols_to_add, by = "id")
    )
  } else {
    return(NULL)
  }
}

#' @export
#' @title Get Keywords
#' @description Gets keywords to feed to `active_EM()`,
#' depending on on scheme type.
#' @param docs Documents table, same as for `active_EM()`.
#' @param dfm Quanteda document-feature matrix matching `docs`.
#' @param num_keywords Number of keywords selected for each class.
#' @param scheme Keyword selection scheme. "max_eta_raw" finds selects based
#' on maximum eta values for positive and negative classes. "max_eta_ratio"
#' selects based on ratio between eta values for positive and negative classes.
#' @param verbose If 'TRUE', prints out keywords to console.
#' @return List of length 2. First element is vector of keywords for negative class,
#' second element is vector of keywords for positive class.
get_keywords <- function(docs, dfm, num_keywords = 10,
                         scheme = c("max_eta_raw", "max_eta_ratio"),
                         verbose = TRUE) {

  ## use naive bayes step to get wrod probabilities
  word_prob_mtx <- get_true_eta(docs, dfm, n_class)
  neg_word_prob <- word_prob_mtx[, 1]
  pos_word_prob <- word_prob_mtx[, 2]

  if (scheme == "max_eta_raw") {
    ## get word strings from top pos and neg words
    ordered_pos_word_prob <- pos_word_prob[order(-pos_word_prob)]
    ordered_neg_word_prob <- pos_word_prob[order(-neg_word_prob)]
    top_pos_words <- names(ordered_pos_word_prob)[1:num_keywords]
    top_neg_words <- names(ordered_neg_word_prob)[1:num_keywords]
    output_list <- list(top_neg_words, top_pos_words)
  } else if (scheme == "max_eta_ratio") {
    neg_word_ratio <- neg_word_prob - pos_word_prob
    neg_word_ratio_ordered <- neg_word_ratio[order(-neg_word_ratio)]
    top_neg_ratio <- names(neg_word_ratio_ordered)[1:num_keywords]
    pos_word_ratio <- pos_word_prob - neg_word_prob
    pos_word_ratio_ordered <- pos_word_ratio[order(-pos_word_ratio)]
    top_pos_ratio <- names(pos_word_ratio_ordered)[1:num_keywords]
    output_list <- list(top_neg_ratio, top_pos_ratio)
  }

  if (verbose) {
    print(output_list)
  }

  return(output_list)
}

get_true_eta <- function(docs, dfm, n_class, n_cluster) {
  ## create a full class matrix
  docs$Class_1 <- 0
  docs$Class_2 <- 0

  ## for (i in 1:nrow(docs)) {
  ##   if (docs[i, "label"] == 1) {
  ##     docs[i, "Class_2"] <- 1
  ##   } else {
  ##     docs[i, "Class_1"] <- 1
  ##   }
  ## }

  docs <- docs %>%
    mutate(
      Class_1 = ifelse(label == 1, 0, 1),
      Class_2 = ifelse(label == 1, 1, 0)
    )

  class_matrix <- get_class_matrix(
    docs, n_class = 2, hand_labeled_index = docs$id,
    doc_name = "text", index_name = "id"
  )

  ## use naive bayes step to get wrod probabilities
  return(get_word_prob_NB(dfm, class_matrix))
}
