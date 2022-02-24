##################################################
## Project: active
## Script purpose: Functions for Output Analysis
## Date: 2019/7/11
## Author: Mitchell Bosley
##################################################

#' @importFrom foreach "%dopar%"
#' @importFrom foreach "%do%"

#' @export
get_conf_matrix <- function(docs, labelsName = "label",
                            index_name = "id", labeledIndex = NULL,
                            levels = c(0, 1), n_class = 2,
                            useLabeled = TRUE) {

    if (useLabeled == F) {
        docs <- docs %>%
            dplyr::filter(!(!!dplyr::sym(index_name)) %in% labeledIndex)
    }
    margin <- 1 / n_class
    objects <- c()
    for (i in 1:length(levels)) {
        for (j in 1:length(levels)) {
            objects <- c(
                objects, docs %>%
                    dplyr::filter(!!dplyr::sym(paste0("Class_", j)) > margin &
                               !!dplyr::sym(labelsName) == levels[i]) %>%
                    nrow()
            )
        }
    }

    conf_matrix <- matrix(objects, ncol=n_class, byrow=T, dimnames=list(levels, levels))

    return(conf_matrix)

}

#' @export
get_classification_accuracy <- function(conf_matrix) {

    total <- 0
    for (i in 1:ncol(conf_matrix)) {
        total <- total + conf_matrix[i, i]
    }

    return(total / sum(conf_matrix))
}

#' @export
get_precision_binary <- function(conf_matrix) {

  if (sum(conf_matrix[ , 2]) == 0) {
    precision <- 0
  } else {
    precision <- conf_matrix[2, 2] / sum(conf_matrix[ , 2])
  }

  return(precision)

}

#' @export
get_recall_binary <- function(conf_matrix) {

  if (sum(conf_matrix[ , 2]) == 0) {
    recall <- 0
  } else {
    recall <- conf_matrix[2, 2] / sum(conf_matrix[2 , ])
  }


  return(recall)
}

#' @export
get_F1_binary <- function(conf_matrix, precision=NULL, recall=NULL) {

  if (is.null(precision) == T){
    precision <- get_precision_binary(conf_matrix)
  }

  if (is.null(recall) == T) {
    recall <- get_recall_binary(conf_matrix)
  }

  if (precision != 0 & recall != 0) {
    F1 <- 2 * ((precision * recall) / (precision + recall))
  } else {
    F1 <- 0
  }

  return(F1)

}

#' @export
get_precision_binary_weighted <- function(conf_matrix, FP_w) {
  
  precision <- conf_matrix[2,2] / (conf_matrix[2,2] + FP_w)
  
  return(precision)
  
}

#' @export
get_recall_binary_weighted <- function(conf_matrix, FN_w) {
  
  if (sum(conf_matrix[ , 2]) == 0) {
    recall <- 0
  } else {
    
    
    recall <- conf_matrix[2,2] / (conf_matrix[2,2] + FN_w)
  }
  
  
  return(recall)
}





#' @export
get_conf_matrix_weight <- function(docs, labelsName, index_name=NULL,
                                   labeledIndex=NULL, levels=c(0, 1),
                                   n_class=2, n_cluster=2, useLabeled=T) {
  #' @title Obtain weighted False negative and false positive in the Confusio nmatrix
  if (useLabeled == F) {
    docs <- docs %>%
      dplyr::filter(!(!!dplyr::sym(index_name)) %in% labeledIndex)
  }

  pos_class <- paste0("Class_", n_cluster)
  neg_class <- c()
  for (i in 1:(n_cluster - 1)) {
    neg_class <- c(neg_class, paste0("Class_", i))
  }

  TN <- sum(docs[, paste0(pos_class)] < 1 / n_cluster & docs[, labelsName] == 0)
  FP <- sum(docs[, paste0(pos_class)] >= 1 / n_cluster & docs[, labelsName] == 0)
  FN <- sum(docs[, paste0(pos_class)] < 1 / n_cluster & docs[, labelsName] == 1)
  TP <- sum(docs[, paste0(pos_class)] >= 1 / n_cluster & docs[, labelsName] == 1)
  objects <- objects_weighed <- c()

  neg_class_sum <- numeric(nrow(docs))
  objects <- c(TN, FP, FN, TP)
  for (class in neg_class) {
    neg_class_sum <- neg_class_sum + docs[, paste0(class)]
  }
  
  # weighted FP and FN
  FP_weighted <- sum((docs[, paste0(pos_class)] >= 1 / n_cluster) &
                     (docs[,labelsName] == 0) * (docs[, paste0(pos_class)]))
  FN_weighted <- sum((docs[, paste0(pos_class)] < 1 / n_cluster) &
                     (docs[,labelsName] == 1) * (neg_class_sum))
  
  
  conf_matrix <- matrix(objects, ncol=n_class, byrow=T, dimnames=list(levels, levels))
  
  return(list(conf_matrix=conf_matrix, FP_w=FP_weighted, FN_w=FN_weighted))
  
}
#' @export
get_empirical_dist <- function(docs, docsType, docName, index_name, labelsName, n_extreme=50,
     dfm=NULL, extreme_type="abs_diff", ...) {

    #' @title Get True Empirical Word Distributions.
    #' @description Get data showing the true empirical distribution of word likelihoods from labeled data.

    # get index of docs with positive classification
    pos_index <- docs %>%
        dplyr::filter(!!sym(labelsName)==1) %>%
        dplyr::pull(id)

    # get dfm of all sampled documents

    if (is.null(dfm)) {
        dfm <- get_dfm(docs=docs, docName=docName, index_name=index_name,
                       stem=stem, trimPct=trimPct, minDocFreq=minDocFreq,
                       idfWeight=idfWeight, ngrams=ngrams)
    }


    # split dfm into positive and negative documents
    pos_dfm <- dfm[which(rownames(dfm) %in% pos_index), ]
    neg_dfm <- dfm[-which(rownames(dfm) %in% pos_index), ]

    # get empirical distribution of word probabilities by class
    Class_1 <- log((1 + quanteda::colSums(neg_dfm)) /
                       (quanteda::nfeat(neg_dfm) + sum(quanteda::colSums(neg_dfm))))
    Class_2 <- log((1 + quanteda::colSums(pos_dfm)) /
                       (quanteda::nfeat(pos_dfm) + sum(quanteda::colSums(pos_dfm))))

    # get class probabilities
    num_pos <- docs %>%
        filter(!!sym(labelsName)==1) %>%
        nrow()
    class_probs <- log(c(num_pos, nrow(docs) - num_pos)) - log(2 + nrow(docs))

    # get term sparsity statistic
    term_sparsity <- get_term_sparsity(dfm)

    # turn data into tibble for ggplots
    dat <- cbind(Class_1, Class_2, term_sparsity)
    word_data <- tibble::add_column(as_tibble(dat, rownames="term"), docsType)

    # find words with most extreme likelihoods
    if (extreme_type == "abs_diff") {
        extreme_words_pos <- word_data %>%
            dplyr::mutate(diff = abs(Class_2 - Class_1)) %>%
            dplyr::filter(Class_2 > Class_1) %>%
            dplyr::arrange(desc(diff)) %>%
            dplyr::slice(1:n_extreme) %>%
            dplyr::pull(term)

        extreme_words_neg <- word_data %>%
            dplyr::mutate(diff = abs(Class_2 - Class_1)) %>%
            dplyr::filter(Class_2 < Class_1) %>%
            dplyr::arrange(desc(diff)) %>%
            dplyr::slice(1:n_extreme) %>%
            dplyr::pull(term)
    } else if (extreme_type == "log_ratio") {
        # find extreme words wth log likelihood ratio
        extreme_words_pos <- word_data %>%
            dplyr::mutate(ratio = Class_2 / Class_1) %>%
            dplyr::filter(ratio > 0) %>%
            dplyr::arrange(ratio) %>%
            dplyr::slice(1:n_extreme) %>%
            dplyr::pull(term)

        extreme_words_neg <- word_data %>%
            dplyr::mutate(ratio = Class_1 / Class_2) %>%
            dplyr::filter(ratio > 0) %>%
            dplyr::arrange(ratio) %>%
            dplyr::slice(1:n_extreme) %>%
            dplyr::pull(term)
    }
    # trim non-extreme words from dfms, and find the number of documents for
    # each class that contain at least one extreme word from the vector of extreme words
    vec <- which(colnames(pos_dfm) %in% extreme_words_pos)
    pos_classification_pct_in <- pos_dfm %>%
        quanteda::dfm(select=extreme_words_pos) %>%
        quanteda::dfm_subset(quanteda::rowSums(.) > 0) %>%
        quanteda::ndoc(.) / quanteda::ndoc(pos_dfm)

    pos_classification_pct_out <-  neg_dfm %>%
        quanteda::dfm(select=extreme_words_pos) %>%
        quanteda::dfm_subset(quanteda::rowSums(.) > 0) %>%
        quanteda::ndoc(.) / quanteda::ndoc(pos_dfm)

    neg_classification_pct_in <- neg_dfm %>%
        quanteda::dfm(select=extreme_words_neg) %>%
        quanteda::dfm_subset(quanteda::rowSums(.) > 0) %>%
        quanteda::ndoc(.) / quanteda::ndoc(pos_dfm)

    neg_classification_pct_out <- pos_dfm %>%
        quanteda::dfm(select=extreme_words_neg) %>%
        quanteda::dfm_subset(quanteda::rowSums(.) > 0) %>%
        quanteda::ndoc(.) / quanteda::ndoc(pos_dfm)

    # create a dummy variable indicating whether a term is extreme or not
    word_data <- word_data %>%
        dplyr::mutate(extreme = dplyr::case_when(
            term %in% extreme_words_pos ~ "Extreme Positive",
            !(term %in% c(extreme_words_pos, extreme_words_neg)) ~ "Not Extreme",
            term %in% extreme_words_neg ~ "Extreme Negative"
        ))

    return(list(word_data=word_data, class_probs=class_probs,
                extreme_words_pos=extreme_words_pos,
                extreme_words_neg=extreme_words_neg,
                pos_classification_pct_in=pos_classification_pct_in,
                pos_classification_pct_out=pos_classification_pct_out,
                neg_classification_pct_in=neg_classification_pct_in,
                neg_classification_pct_out=neg_classification_pct_out))
}

#' @export
get_fixed_words_mtx <- function(data) {

    #' @title Get a Matrix of Fixed Words and Class Probabilities.
    #' @description Gets a matrix of extreme words and respective word probabilities. Used to
    #' supply to EM algorithm to fix word probabilities.

    out <- data %>%
        filter(extreme != "Not Extreme") %>%
        select(term, neg_vec, pos_vec)

    out2 <- as.matrix(out[, 2:3])
    rownames(out2) <- out[["term"]]
    colnames(out2) <- c("Class_1", "Class_2")

    return(out2)
}

#' @export
run_unsupervised_EM <- function(docs, docName, docsType, index_name, labelsName, dfm=NULL,
                                extreme_words_pos=NULL, extreme_words_neg=NULL,
                                stem=T, trimPct=0, minDocFreq=1, idfWeight=F,
                                ngrams=1, n_extreme=50, n_class=2, export_plot=T,
                                export_all=T, class_prob=NULL, word_prob=NULL, ...) {

    #' @title Run Unsupervised EM
    #' @description Runs supervised EM on a corpus of documents cleaned by \code{clean_data()}.
    #' @return  Document level probabilities, word probabilities, scatterplot of word probablities.

    if (is.null(dfm)) {
        dfm <- get_dfm(docs=docs, docName=docName, index_name=index_name,
                       stem=stem, trimPct=trimPct, minDocFreq=minDocFreq,
                       idfWeight=idfWeight, ngrams=ngrams)
    }

    term_sparsity <- get_term_sparsity(dfm)

    # Run EM Algorithm
    output <- EM(.D_test=dfm, .supervise=F, .class_prob=class_prob,
                 .word_prob=word_prob, .export_all=export_all)

    # Assign colnames if needed
    if (is.null(colnames(output$classLik[[1]])) | is.null(colnames(output$eta[[1]]))) {
        class_names <- get_classes(n_class)
        if (length(output$classLik) > 1) {
            for (i in 1:length(output$classLik)) {
                colnames(output$classLik[[i]]) <- class_names
            }
        } else {
            colnames(output$classLik) <- class_names
        }

        if (length(output$eta) > 1) {
            for (i in 1:length(output$eta)) {
                colnames(output$eta[[i]]) <- class_names
            }
        } else {
            colnames(output$eta) <- class_names
        }

    }

    # match EM results to documents
    # docs <- match_EM_to_docs(docs=docs, EMoutput=output$classLik[[length(output$classLik)]], classes=get_classes(n_class),
    #                          docName=docName, index_name=index_name, labelsName="label")
    #
    # Restructure data
    if (export_all == T) {
        word_data <- tibble()
        for (i in seq_along(output$eta)) {
            word_data <- output$eta[[i]] %>%
                cbind(., term_sparsity, EM_iteration=i) %>%
                as.matrix() %>%
                as_tibble(rownames="term") %>%
                bind_rows(word_data, .)
        }

        docs_data <- tibble()
        for (i in seq_along(output$classLik)) {
            docs_data <- output$classLik[[i]] %>%
                as.matrix() %>%
                as_tibble(rownames="id") %>%
                bind_cols(., label=docs %>% pull(!!dplyr::sym(labelsName))) %>%
                mutate(EM_iteration=i) %>%
                bind_rows(docs_data, .)
        }

    } else {
        word_data <- output$eta[[1]] %>%
            cbind(., term_sparsity, EM_iteration=i) %>%
            as.matrix() %>%
            as_tibble(rownames="term")
    }

    # Add extreme category
    if (!is.null(extreme_words_pos) & !is.null(extreme_words_neg)) {
        word_data <- word_data %>%
            mutate(extreme = case_when(
                term %in% extreme_words_pos ~ "Extreme Positive",
                !(term %in% c(extreme_words_pos, extreme_words_neg)) ~ "Not Extreme",
                term %in% extreme_words_neg ~ "Extreme Negative"
            ))
    }

    return(list(docs_data=docs_data, word_data=word_data))

}

#' @export
# KL divergence function
KL <- function(p, q){

    #' @title KL Divergence Function
    #' @description Gets KL divergence between two distributions.
    #'
    #' @param p Empirical distribution.
    #' @param q Posterior distribution.
    #'

    return(sum(exp(p)*(p - q)))
}

#' @export
tuning_algorithm <- function(docs, docName, index_name, labelsName, measure="accuracy", max_runs=30, seed=123,
                             stem=c(F, T), ngrams=c(1:3), trimPct=c(0, 0.0001, 0.001, 0.01),
                             minDocFreq=c(1, 2, 3), idfWeight=c(F, T), removeStopWords=c(F, T), minChar=1:4,
                             weight=seq(0, 1, .1), bound=seq(0, 1, .1), initSize=5, maxActive=20,
                             maxQuery=1, queryType=c("basic_entropy", "random"), forceSweep=NA, dfmOut=F) {

    #' @title Active EM Parameter Tuning Algorithm
    #' @description Tries to find optimal combination of pre-processing and active algorithm parameters.

    params <- expand.grid(ngrams=ngrams, stem=stem, trimPct=trimPct, minDocFreq=minDocFreq,
                          idfWeight=idfWeight, removeStopWords=removeStopWords, minChar=minChar,
                          weight=weight, bound=bound, initSize=initSize, maxActive=maxActive,
                          maxQuery=maxQuery, queryType=queryType, accuracy=NA, F1=NA, KL=NA)

    params_record <- data.frame()
    last <- 0
    # set conditional counter
    i <- 1
    # set explicit counter
    total_runs <- 0
    max_runs <- 100

    if (dfmOut) {
        dfm <- get_dfm(docs=docs, docName=docName, index_name=index_name, stem=params[i, "stem"], ngrams=params[i, "ngrams"],
                       trimPct=params[i, "trimPct"], minDocFreq=params[i, "minDocFreq"], idfWeight=params[i, "idfWeight"],
                       removeStopWords=params[i, "removeStopWords"], minChar=params[i, "minChar"])
    }

    repeat {

        if (!dfmOut) {
            dfm <- get_dfm(docs=docs, docName=docName, index_name=index_name, stem=params[i, "stem"], ngrams=params[i, "ngrams"],
                           trimPct=params[i, "trimPct"], minDocFreq=params[i, "minDocFreq"], idfWeight=params[i, "idfWeight"],
                           removeStopWords=params[i, "removeStopWords"], minChar=params[i, "minChar"])
        }

        active_output <- active_EM(docs=docs, weight=params[i, "weight"], docName=docName, maxActive=params[i, "maxActive"],
                                   maxQuery=params[i, "maxQuery"], initSize=params[i, "initSize"], index_name=index_name, labelsName=labelsName,
                                   bound=params[i, "bound"], handLabel=F, dfm=dfm, forceList=T, queryType=params[i, "queryType"], seed=seed)

        cf <- get_conf_matrix(docs=active_output$docs, labelsName=labelsName, index_name=index_name,
                              labeledIndex=active_output$handLabeledIndex, n_class=2, levels=c(0,1), useLabeled=F)

        params[i, "accuracy"] <- get_classification_accuracy(cf)
        params[i, "F1"] <- get_F1_binary(cf)
        params[i, "KL"] <-  get_empirical_dist(docs=docs, docsType="", docName=docName,
                                               index_name=index_name, labelsName=labelsName,
                                               dfm=dfm)$word_data %>%
            select(Class_1, Class_2) %>%
            as.matrix() %>%
            KL(p=., q=active_output$EMoutput[[length(active_output$EMoutput)]]$eta)

        # params[i, "measure"] <- dplyr::case_when(
        #     measure == "accuracy" ~ get_classification_accuracy(cf),
        #     measure == "F1" ~ get_F1_binary(cf),
        #     measure == "KL" ~ get_empirical_dist(docs=docs, docsType="", docName=docName,
        #                                          index_name=index_name, labelsName=labelsName,
        #                                          dfm=dfm)$word_data %>%
        #         select(Class_1, Class_2) %>%
        #         as.matrix() %>%
        #         KL(p=., q=active_output$EMoutput[[length(active_output$EMoutput)]]$eta)
        # )
        params_record <- bind_rows(params_record, params[i, ])

        cond <- function(measureVal, last, measure) {
            if (measure != "KL") {
                out <- measureVal < last
            } else {
                out <- abs(0 - measureVal) > abs(0 - last)
            }
            return(out)
        }

        # Main decision making code
        if (i > 1) {
            move <- which(params[i, 1:(ncol(params) - 3)] != params[i-1, 1:(ncol(params) - 3)])
            var_name <- colnames(params)[move]
            if (var_name %in% forceSweep |
                is.factor(params[, var_name]) |
                is.logical(params[, var_name])) {
                if (cond(params[i, paste0(measure)], last, measure)) {
                    params <- params %>% filter(!!sym(var_name) != params[i, var_name])
                } else {
                    last <- params[i, paste0(measure)]
                    params <- params %>% filter(!(!!sym(var_name) == params[i - 1, var_name]))
                }
            } else if (is.numeric(params[, var_name])) {
                if (cond(params[i, paste0(measure)], last, measure)) {
                    params <- params %>% filter(!!sym(var_name) == params[i - 1, var_name])
                } else {
                    last <- params[i, paste0(measure)]
                    params <- params %>% filter(!!sym(var_name) >= params[i, var_name])
                }
            }
            i <- which(is.na(params[[paste0(measure)]]))[1]
        } else {
            last <- params[i, paste0(measure)]
            i <- i + 1
        }

        total_runs <- total_runs + 1
        print(params_record)
        if (total_runs == max_runs | !(T %in% is.na(params[[paste0(measure)]]))) break

    }

    return(list(params=params, params_record=params_record))

}


#' @export
get_KL <- function(sample_data, docName='text', index_name='id', labelsName='label', docsType='bbc',
                   stem=T, ngrams=1, trimPct=0.0001, minDocFreq=2, idfWeight=F, n_extreme=50, iter=100){
  # get KL divergence from the word-distribution from unsupervised EM
  # to the empirical word-distribution

  sample_dfm <- get_dfm(docs=sample_data, docName=docName, index_name=index_name,
                        stem=stem, ngrams=ngrams, trimPct=trimPct,
                        minDocFreq=minDocFreq, idfWeight=idfWeight)
  # empirical distribution of words
  empdist <- get_empirical_dist(docs=sample_data, docsType=docsType, docName=docName,
                                index_name=index_name, labelsName=labelsName,
                                stem=stem, trimPct=trimPct,
                                minDocFreq=minDocFreq, idfWeight=idfWeight,
                                ngram=ngram, n_extreme=n_extreme, dfm=sample_dfm)
  empdist <- as.data.frame(empdist$word_data)
  rownames(empdist) <- empdist$term
  empdist <- empdist[,c('Class_1', 'Class_2')]

  # empirical distribution of classes
  empdist_class <- log(c(mean(sample_data$label), 1- mean(sample_data$label)))

  # posterior from unsupervised EM
  post <- as.data.frame(as.matrix(EM(.D_test = data.matrix(sample_dfm), .supervise = F,
                                     .word_prob=as.matrix(empdist), .class_prob=empdist_class)$eta))
  colnames(post) <- c('Class_1', 'Class_2')

  # KL divergence from the posterior to the empirical distribution for class 1
  kl_post <- KL(empdist$Class_1, post$Class_1)

  kl_prior_vec <- rep(0, iter)
  for (i in 1:iter){
    # prior
    prior <- data.frame(Class_1 = as.vector(gtools::rdirichlet(n=1, alpha = rep(2, ncol(sample_dfm)))),
                        Class_2 = as.vector(gtools::rdirichlet(n=1, alpha = rep(2, ncol(sample_dfm)))))
    prior <- log(prior)
    kl_prior_vec[i] <- KL(empdist$Class_1, prior$Class_1)
  }

  res <- list(kl_post = kl_post, kl_prior = kl_prior_vec)

  return(res)
}


#' @export
check_extreme_words_inlabel <- function(data, docsType, docName, index_name, labelsName, dfm,
                                        n_extreme=10, iter=1000, sample_size=100, sampling_method=NULL){


  #' @title Check how often the extreme words are in the labeled documents.
  #' @description  check the distribution of the extreme words in the labeled documents of size 100
  # the extreme wrods are the words whose difference in the class porterior between one class
  # and the other class is the biggest. They are defined using the entire data.
  # this function returns the dataframe of the
  #'
  #' @param n_extreme number of extreme words
  #' @param iter the number of sampling
  #'

  # get extreme words
  emp_dist <- get_empirical_dist(docs=data,
                                 docsType=docsType, docName=docName, index_name=index_name, labelsName=labelsName,
                                 n_extreme=n_extreme, dfm=dfm)
  extreme_neg <- emp_dist$extreme_words_neg
  extreme_pos <- emp_dist$extreme_words_pos


  if(is.null(sampling_method)){

    # storage
    sample_extreme_neg_df <- data.frame(matrix(NA, nrow=iter, ncol=length(extreme_neg)))
    colnames(sample_extreme_neg_df) <- extreme_neg
    sample_extreme_pos_df <- data.frame(matrix(NA, nrow=iter, ncol=length(extreme_pos)))
    colnames(sample_extreme_pos_df) <- extreme_pos

    # sample 100 labeled documents and get the distribution of extreme word
    for (i in 1:iter){
      select <- sample(x=1:nrow(dfm), size=sample_size, replace=F)
      sample_word_dist <- colSums(dfm[select,])
      sample_extreme_neg_df[i,] <- sample_word_dist[extreme_neg]
      sample_extreme_pos_df[i,] <- sample_word_dist[extreme_pos]

      if (i %% 100 == 0){
        print(paste0(i, ' sampled...'))
      }
    }
  }else if(sampling_method=='bootstrap'){

    if (sample_size != nrow(dfm)){
      warning('it is recommended that the sample size is the same as the data size')
    }

    # storage
    sample_extreme_neg_df <- data.frame(matrix(NA, nrow=iter, ncol=length(extreme_neg)))
    colnames(sample_extreme_neg_df) <- extreme_neg
    sample_extreme_pos_df <- data.frame(matrix(NA, nrow=iter, ncol=length(extreme_pos)))
    colnames(sample_extreme_pos_df) <- extreme_pos

    for (i in 1:iter){
      select <- sample(x=1:nrow(dfm), size=sample_size, replace=T)
      sample_word_dist <- colSums(dfm[select,])
      sample_extreme_neg_df[i,] <- sample_word_dist[extreme_neg]
      sample_extreme_pos_df[i,] <- sample_word_dist[extreme_pos]

      if (i %% 100 == 0){
        print(paste0(i, ' sampled...'))
      }
    }


  }else if(sampling_method=='jackknife'){

    if (iter != nrow(dfm)){
      warning("the number of iteration is set to  the same as the data size in jackknife,
              regardless of the 'iter' value")
    }
    iter <- nrow(dfm)

    # storage
    sample_extreme_neg_df <- data.frame(matrix(NA, nrow=iter, ncol=length(extreme_neg)))
    colnames(sample_extreme_neg_df) <- extreme_neg
    sample_extreme_pos_df <- data.frame(matrix(NA, nrow=iter, ncol=length(extreme_pos)))
    colnames(sample_extreme_pos_df) <- extreme_pos

    for (i in 1:nrow(dfm)){
      select <- 1:nrow(dfm)[-i]
      sample_word_dist <- colSums(dfm[select,])
      sample_extreme_neg_df[i,] <- sample_word_dist[extreme_neg]
      sample_extreme_pos_df[i,] <- sample_word_dist[extreme_pos]

      if (i %% 100 == 0){
        print(paste0(i, ' sampled...'))
      }
    }

  }



  return(list(pos=sample_extreme_pos_df, neg=sample_extreme_neg_df))
}


#' @export
plot_hist_extreme <- function(extreme_df, docsType){
  # use with check_extreme_words_inlabel function

  pos_hist1 <- ggplot(extreme_df) +
    geom_histogram(aes_string(x = colnames(extreme_df)[1])) +
    ggtitle(paste0(docsType, ': Extreme word (positive) 1'))
  pos_hist2 <- ggplot(extreme_df) +
    geom_histogram(aes_string(x = colnames(extreme_df)[2])) +
    ggtitle('Extreme word (positive) 2')
  pos_hist3 <- ggplot(extreme_df) +
    geom_histogram(aes_string(x = colnames(extreme_df)[3])) +
    ggtitle('Extreme word (positive) 3')
  pos_hist4 <- ggplot(extreme_df) +
    geom_histogram(aes_string(x = colnames(extreme_df)[4])) +
    ggtitle('Extreme word (positive) 4')
  res <- plot_grid(pos_hist1, pos_hist2, pos_hist3, pos_hist4)
  return(res)
}


#' @export 
balancing_data <- function(data, pos_ratio, size = nrow(data), seed=1990){
  
  #' @title Create a (un)balanced dataset with specified size and ratio of positive documents
  #'
  #' @param pos_ratio ratio of positive document ()
  #' @param iter the number of sampling

  if (pos_ratio == "true_ratio") {
    return(data)
  } else {
    
    pos_ratio <- as.numeric(pos_ratio)

    pos <- data[data$label == 1,]
    neg <- data[data$label == 0,]

    pos_size <- floor(size * pos_ratio)
    neg_size <- size - pos_size

    while (neg_size > nrow(neg) | pos_size > nrow(pos)) {
      size <- size - 1
      pos_size <- floor(size * pos_ratio)
      neg_size <- size - pos_size
    }

    if (pos_size <= nrow(pos) & neg_size <= nrow(neg)) {
      set.seed(seed)
      pos_sample <- pos %>% sample_n(pos_size)
      set.seed(seed)
      neg_sample <- neg %>% sample_n(neg_size)
    } else {
      warning('Cannot produce sample data with specified size and ratio.')
    }

    return(rbind(pos_sample, neg_sample))
  }
}

#' @export
sample_data <- function(data, n=NULL, pos_prob){
  #' @title get a sample of dataset with a specific proportion of positive documents.

  # If the requested proportion of positive dosc are larger than the 
  # proportion of positive docs in the population data
  # then use all the  positive docs and reduce the negative docs
  if (is.null(n)){
    pop_pos_prob  <- mean(data$label)
    if (pos_prob > pop_pos_prob){
      pos_n <- sum(data$label)
      neg_n <- (1-pos_prob)/pos_prob * pos_n
    }else{ 
      neg_n <- nrow(data) - sum(data$label)
      pos_n <- pos_prob / (1-pos_prob) * neg_n
    }
  }else{
    pos_n <- n * pos_prob
    neg_n <- n * (1-pos_prob)
  }


  pos_data <- data %>% filter(label == 1)
  neg_data <- data %>% filter(label == 0)
  pos_sample <- pos_data %>% sample_n(size=pos_n)
  neg_sample <- neg_data %>% sample_n(size=neg_n)
  out <- rbind(pos_sample, neg_sample)

  # permutate the order of the rows to avoid stacking pos and neg docs
  perm_row <- sample(x=seq(1, nrow(out)), size=nrow(out))
  out <- out[perm_row,]
  return(out)
  
}



#' @export
reshapedocs <- function(output){
  #' @title reshape docs object such that the prediction at each iteration is stored
  #' @param "output" is the entire output (not "output" element of the list) from the run_models_fast function

  return(lapply(output$output, reshapedocs_helper))
}

reshapedocs_helper <- function(output_model){
  #' @title (helper) reshape docs object such that the prediction at each iteration is stored


  output_model_docs <- output_model$docs
  output_model_outdocs <- output_model$out_docs

  # create a storage
  alldocs_acrossiter <- tibble(
    text=c(output_model_docs[[1]]$text, output_model_outdocs[[1]]$text),
    id=c(output_model_docs[[1]]$id, output_model_outdocs[[1]]$id),
    label=c(output_model_docs[[1]]$label, output_model_outdocs[[1]]$label),
    train=c(rep(1,nrow(output_model_docs[[1]])), rep(0, nrow(output_model_outdocs[[1]])))
  )

  # store predicted probability of being a pos doc for each iteration
  n_iter <- length(output_model_docs)
  n_docs <- nrow(alldocs_acrossiter)

  iter <- matrix(0, nrow=n_docs, ncol=n_iter)
  colnames(iter) <- paste0("iter", seq(1,ncol(iter)))
  for (i in 1:length(output_model_docs)){
    iter[,i] <- c(output_model_docs[[i]]$Class_2, output_model_outdocs[[i]]$Class_2)
  }

  alldocs_acrossiter <- as_tibble(cbind(alldocs_acrossiter, iter))
  
  #For hand labeled documents, replace the predicted probability with -1 
  for (i in 1:n_iter){
    if (i == 1){
      handlabel <- output_model$handLabeledIndex[[1]]
    }else{
      handlabel <- setdiff(output_model$handLabeledIndex[[i]],
                           output_model$handLabeledIndex[[i-1]])
    }
    # flag which docs are labeled at each iteration
    handlabel_logical <- alldocs_acrossiter$id %in% handlabel
    alldocs_acrossiter[handlabel_logical,paste0("iter",seq(i,n_iter))] <- -1 

  }

  return(alldocs_acrossiter)
}

#' @export
get_results_matrix <- function(out, out_sample = FALSE, model_name) {
#' @title Get Results Matrix for active_EM Output.
#' @description Translates activeEM() output object
#' into result matrix for graphing.
#' @param out list; Output from `active_EM` function.
#' @param out_sample logical; Indicates whether to get output matrix
#' out of sample test. Will only work if 'active_EM' exports an `out_docs`
#' object.
#' @param model_name string; Name of model.
#' @return mtx; Object for feeding into `get_figure_grid` function.
#' @note Ensure that 'exportAll = T' in `active_EM` function.

  obj <- tibble(
    model_name = rep(model_name, length(out$docs)),
    iter = 0,
    docs_labeled = 0,
    accuracy_in = 0,
    precision_in = 0,
    recall_in = 0,
    F1_in = 0,
    pos_pi = 0
  )

  if (out_sample) {
    obj <- bind_cols(
      obj,
      tibble(
        accuracy_out = rep(0, length(out$docs)),
        precision_out = 0,
        recall_out = 0,
        F1_out = 0
      )
    )
  }

  for (i in seq_along(out$docs)) {
    cf <- get_conf_matrix(
      docs = out$docs[[i]],
      labeledIndex = out$handLabeledIndex[[i]],
      useLabeled = FALSE
    )
    obj[i, ]$iter <- i - 1
    obj[i, ]$docs_labeled <- length(out$handLabeledIndex[[i]])
    obj[i, ]$accuracy_in <- get_classification_accuracy(cf)
    obj[i, ]$precision_in <- get_precision_binary(cf)
    obj[i, ]$recall_in <- get_recall_binary(cf)
    obj[i, ]$F1_in <- get_F1_binary(cf)
    ## obj[i, ]$pos_pi <- exp(out$EMoutput[[i]]$pi)[2]
  }

  if (out_sample) {
     for (i in seq_along(out$docs)) {
       cf_out <- get_conf_matrix(
         docs = out$out_docs[[i]],
         labeledIndex = out$handLabeledIndex[[i]],
         useLabeled = FALSE
       )
       obj[i, ]$accuracy_out <- get_classification_accuracy(cf_out)
       obj[i, ]$precision_out <- get_precision_binary(cf_out)
       obj[i, ]$recall_out <- get_recall_binary(cf_out)
       obj[i, ]$F1_out <- get_F1_binary(cf_out)
     }
  }

  metadata <- as_tibble(out$metadata, .rows = nrow(obj))
  obj <- bind_cols(obj, metadata)

  return(obj)
}

#' @export
get_label_prop_fig <- function(out, dataset_name, index_name = "id") {
#' @title Get Label Proportion Figure
#' @description Produces a figure that proportion of positive to negative
#' queries across active iteration.
  results <- tibble(iter = NA, pos_label = NA, neg_label = NA)
  i <- 1
  for (ind in out$handLabeledIndex) {
    results[i, "pos_label"] <- out$docs[[1]] %>%
      filter(!!dplyr::sym(index_name) %in% ind) %>%
      pull(label) %>% sum()
    results[i, "neg_label"] <- length(ind) - results[i, "pos_label"]
    results[i, "iter"] <- i - 1
    i <- i + 1
  }
  results <- reshape2::melt(data = results, id.vars = "iter",
                            measure.vars = c("pos_label", "neg_label"))
  plot <- ggplot(results) +
    aes(x = iter, y = value, fill = variable) +
    geom_area() +
    scale_x_continuous(breaks = integer_breaks()) +
    ggtitle("Proportion of Labeled Documents by Active Iteration")

  cap_obj <- process_metadata(out$metadata)
  plot <- plot + labs(subtitle=paste("Dataset:", dataset_name),
                      caption = cap_obj)

  return(list(plot = plot, results = results))
}

#' @export
process_metadata <- function(md) {
#' @title Process Model Metadata
#' @description Exports a caption object for figures based on
#' model specifications.
  cap_obj <- paste(
    "Corpus Size: ", md$corpusSize, "\n",
    "Train Size: ", md$trainSize, "\n",
    "Num Clusters: ", md$n_cluster, "\n",
    "Lambda Weight: ", md$lambda, "\n",
    "Pos Class Ratio: ", round(md$pos_ratio, 2), "\n",
    "# of Docs to Initialize: ", md$initSize, "\n",
    "Max Docs Labeled per Iter: ", md$maxQuery
  )
  if (md$queryType != "log_ratio") {
    cap_obj <- paste(cap_obj, "\n", "Sampling Type: ", md$queryType)
  } else {
    cap_obj <- paste(cap_obj, "\n",
                     "Sampling Type: ", md$queryType,
                     "Mu: ", md$mu,
                     "Tau: ", md$tau)
  }
  if (md$lambda_decay) {
    cap_obj <- paste(cap_obj, "\n", "Lambda Decay of",
                     md$ld_rate, "on weighted models.")
  }
  return(cap_obj)
}

#' @export
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

#' @export
get_predict_prop_fig <- function(out, dataset_name, index_name = "id",
                                 doc_labels = FALSE, metadata = TRUE,
                                 arrange_y_axis = FALSE, out_sample = TRUE,
                                 labeled_only = FALSE) {
#' @title Get Prediction Probability Figure
#' @description Produces a figure that charts changes to out of sample
#' classification probability across active iterations.

  i <- 1
  results <- tibble()
  if (out_sample) {
    docs_obj <- out$out_docs
  } else {
    docs_obj <- out$docs
  }

  obj_len <- length(docs_obj)
  last_hl_index <- out$handLabeledIndex[[obj_len]]
  for (docs in docs_obj) {
    if (i == 1) {
      new_labeled <- out$handLabeledIndex[[i]]
    } else {
      new_labeled_ind <- which(!(out$handLabeledIndex[[i]] %in%
                                  out$handLabeledIndex[[i - 1]]))
      new_labeled <- out$handLabeledIndex[[i]][new_labeled_ind]
    }
    if (labeled_only & !out_sample) {
      docs <- docs %>%
        dplyr::filter(!!dplyr::sym(index_name) %in% last_hl_index)
    }
    res <- docs %>%
      dplyr::select(id, label, Class_2) %>%
      dplyr::mutate(labeled = if_else(
                      !!dplyr::sym(index_name) %in% new_labeled, 1, 0
                    )) %>%
      dplyr::mutate(iter = i - 1)
    results <- results %>% dplyr::bind_rows(res)
    i <- i + 1
  }

  if (arrange_y_axis) {
    results <- results %>%
      mutate(id = as.numeric(id)) %>%
      arrange(iter, id) %>%
      mutate(id = as.factor(id))
  }
  results <- results %>% rename(`Pos Class Prob` = Class_2)

  docs_length <- results %>% filter(iter == 1) %>% nrow()
  p1 <- ggplot(results) +
    aes(x = iter, y = id)
    if (docs_length > 1000) {
      p1 <- p1 + geom_raster(aes(fill = `Pos Class Prob`))
    } else {
      p1 <- p1 + geom_tile(aes(fill = `Pos Class Prob`), colour = "white")
    }
  p1 <- p1 + scale_fill_gradient(low = "red",
                                 high = "blue") +
    scale_y_discrete(guide = guide_axis(n.dodge = 1, angle = 0)) +
    scale_x_continuous(breaks = integer_breaks()) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 5))

  if (out_sample) {
    p1 <- p1 + ggtitle("Out Sample Label Prediction Across Active Iteration")
  } else {
    p1 <- p1 + ggtitle("In Sample Label Prediction Across Active Iteration")
  }
  if (!doc_labels) {
    p1 <- p1 + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  if (!out_sample) {
    if (docs_length > 1000) {
      p1 <- p1 + geom_raster(aes(alpha = labeled), fill = "yellow")
    } else {
      p1 <- p1 + geom_tile(aes(alpha = labeled
                               ),
                           fill = "yellow")
    }
  }

  p2 <- ggplot(results) +
    aes(x = 1, y = id)
  if (length(unique(results$label)) == 1) {
    if (docs_length > 1000) {
      p2 <- p2 + geom_raster(aes(fill = "red"))
    } else {
      p2 <- p2 + geom_tile(aes(fill = "red"), colour = "white")
    }
  } else {
    if (docs_length > 1000) {
      p2 <- p2 + geom_raster(aes(fill = label))
    } else {
      p2 <- p2 + geom_tile(aes(fill = label), colour = "white")
    }
    p2 <- p2 + scale_fill_gradient(low = "red",
                                   high = "blue")
  }
  p2 <- p2 + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    ) +
    ggtitle("True \n Labels")

  cap_obj <- process_metadata(out$metadata)

  plot <- p2 + p1 + patchwork::plot_layout(widths = c(0.1, 1))
  plot <- plot + labs(subtitle=paste("Dataset:", dataset_name))
  if (metadata) {
    plot <- plot + labs(caption = cap_obj)
  }
  return(plot)
}

#' @export
read_doc <- function(docs, index_name = "id", id_val) {
  doc <- docs %>%
    filter(!!dplyr::sym(index_name) == id_val) %>%
    pull(text)
  return(doc)
}
