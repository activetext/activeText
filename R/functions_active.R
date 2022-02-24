##################################################
## Project: active
## Script purpose: Functions for Active Learning EM
## Date: 2019/6/29
## Author: Mitchell Bosley
##################################################

#' @importFrom dplyr "%>%"
## Active Learning EM Implementation
#' @export   
active_EM <- function(docs, labels=c(0, 1), docName = "test", indexName = "id",
                      labelsName = NULL, weight = 1, n_class = 2, n_cluster = 2,
                      initIndex = NULL, handLabel = T, bound = 0, maxActive = 5,
                      initSize = 10, maxQuery = 10, lazy_eval = F, forceList = F,
                      counter_on = T, queryType = "basic_entropy",
                      whichOutTest = NULL, seed = NULL, fixed_words = NULL,
                      dfm = NULL, exportAll = F, log_ratio_threshold = 0.001,
                      log_ratio_conv_type = "maximand", mu = 0.0001, tau = 0.0001,
                      regions = "both", lambda_decay = F, ld_rate = 0.2,
                      tune_lambda = FALSE, tune_lambda_prop_init = 0.1,
                      tune_lambda_range = seq(0, 1, 0.1), tune_lambda_k = 10,
                      tune_lambda_parallel = TRUE, NB_init = TRUE,
                      export_val_stats_only = FALSE, model_name = NULL,
                      ...) {

#' @title        Active Learning EM Algorithm
#' @description  Active learning for weighted-EM algorithm. After initial EM algorithm converges,
#'               oracle is queried for labels to documents that the EM algorithm was most unsure of.
#'               This process iterates until max iterations are reached, or there are no
#'               documents in the window of uncertainty.
#'
#' @param docs         [matrix]     Matrix of labeled and unlabeled documents, where each
#'                                  row has index values and a nested Matrix of word tokens.
#' @param labels       [vector]     Vector of character strings indicating classification
#'                                  options for labeling.
#' @param weight       [numeric]    Numeric value between 0 and 1. Used to weight unlabeled
#'                                  documents.
#' @param n_class      [numeric]    Number of classes to be considered.
#' @param docName      [character]  Character string indicating the variable in 'docs'
#'                                  that denotes the text of the documents to be classified.
#' @param indexName    [character]  Character string indicating the variable in 'docs'
#'                                  that denotes the index value of the document to be classified.
#' @param labelsName   [character]  Character string indicating the variable in \code{docs}
#'                                  that denotes the already known labels of the documents.
#'                                  By default, value is set to \code{NULL}.
#' @param initIndex    [vector]     Vector that indicates which documents to use to initialize the
#'                                  algorithm. By default set to \code{NULL}, which causes a random
#'                                  subset of the documents to be selected.
#' @param handLabel    [logical]    Boolean logical value indicating whether to initiate user-input script.
#'                                  If set to \code{FALSE}, and if \code{labelsName} is provided, the script
#'                                  queries the document label directly from the column denoted by \code{labelsName}.
#' @param bound        [numeric]    Minimum bound of entropy to call for additional labelling.
#' @param maxActive    [numeric]    Value of maximum allowed active learning iterations.
#' @param initSize     [numeric]    Value of maximum allowed iterations within the EM algorithm.
#' @param maxQuery     [numeric]    Maximum number of documents queried in each EM iteration.
#' @param lazy_eval    [logical]    If \code{lazy_eval == T}, convergence is measured by comparing changes
#'                                  in log likelihood across model iterations rather than directly computing
#'                                  maximand.
#' @param forceList    [logical]    Switch indicating whether to force the filtering of documents with
#'                                  no entropy. Set to \code{FALSE} by default.
#' @param counter_on   [logical]    Switch indicating whether the progress of each sequence of the EM algorithm
#'                                  is reported. By default set to \code{TRUE}.
#' @param whichOutTest [vector]     Vector of document index labels used to identify documents to be used for
#'                                  out of sample validation of the learned model. Set to \code{NULL} by default.
#'                                  If a vector of labels is provided, the function outputs an additional argument
#'                                  containing classification likelihoods for all documents identified by the vector.
#' @param exportAll    [logical]    Switch indicating whether to export model predictions from each stage of the algorithm.
#'                                  If true, the function exports a list of lists containing all predictions.
#' @param export_val_stats_only     Boolean, indicating whether to export validation stats only from model runs.
#' @param queryType    [string]     String indicating which type of uncertainty sampling to use. Options are \code{"standard_entropy"}
#'                                  or \code{"normalized_entropy"}, \code{"tiered_entropy"}, or \code{"tiered_entropy_weighted"}.
#' @param quantileBreaks [vector]   Vector of break points to distinguish entropy zones. The first value is
#'                                  the break point between the first and second tier, the second is the
#'                                  break point between the second and third tier.
#' @param sampleProps    [vector]   Vector of sampling proportions for each entropy zone. The first value is
#'                                  the proportion of \code{maxQuery} to be sampled from the high entropy region,
#'                                  the second value is the proportion to be sampled from the middle entropy region,
#'                                  and the third value is the proportion to be sampled from the lowest entropy region.
#' @param seed         [numeric]    Sets seed for model.
#' @param fixed_words  [matrix]     Matrix of fixed words with class probabilities, where ncol is the number of classes.
#' @param supervise    [logical]    T if supervised. F is unsupervised.
#' @param dfm          [matrix]     Option to manually supply a dfm from quanteda.
#' @param log_ratio_threshold [numeric] Threshold at which convergence is declared when using `queryType="log_ratio"`.
#' @param log_ratio_conv_type [string] If `queryType="log_ratio"`, this supplies the way that convergence is estimated.
#'                                  Set to `maximand` by default.
#' @param mu                        Parameters for error acceptance with `queryType=log_ratio`.
#' @param tau                       Parameters for error acceptaance with `queryType=log_ratio`.
#' @param regions      [string]     Can be set to "both", "pos", or "neg" to sample from certain regions during log ratio sampling.
#' @param lambda_decay [logical]    Determines whether lambda value decays over active learning iterations or not.
#' @param ld_rate      [float]      If `lambda_decay == TRUE`, sets the rate at which decay occurs.
#' @param tune_lambda  [logical]    Logical value indictating whether to tune lambda values with cross validation over
#'                                  active learning iterations.
#' @param tune_lambda_prop_init [numeric] Float value indicating the proportion of documents to label supply
#' rather than label with EM during lambda tuning.
#' @param tune_lambda_range [vector] Vector of float values, indicating the range of lambda values to search
#' over when tuning lambda at each active iteration.
#' @param tune_lambda_k [integer] Integer value indicating what k-fold level to cross validate at when
#' tuning lambda.
#' @param NB_init [boolean] Indicates whether each active iteration should start with a naive step in the EM
#' or whether to initialize with model predictions from previous active iteration.
#' @param model_name [string] Model name string for exporting when `export_val_stats_only == TRUE`.
#' @param ...                       Additional parameters to pass to `get_dfm` and `EM()` and `get_uncertain_docs()`.
#'
#'
#' @return             [list]       List containing labeled document matrix, prior weights, word likelihoods, and a vector
#'                                  of user-labeled documents ids.

  ## Messages
  ## --------------------------------------------------------------------------

  if (n_cluster == 2 & queryType %in% c("margin_cluster", "basic_entropy_cluster")) {
    queryType <- "basic_entropy"
    message("Cluster sampling only works with greater than two clusters.
Defaulting to basic_entropy sampling scheme.")
  }

  ## Setup
  ## --------------------------------------------------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }

  ## Create containers
  handLabeledIndex <- list()
  EM_docs <- list()
  out_docs <- list()
  output <- list()

  if (is.null(initIndex)) {
    ## Randomly sample documents
    if (!is.null(whichOutTest)) {
      rdmDocs <- docs %>%
        dplyr::filter(!(!!dplyr::sym(indexName)) %in% whichOutTest) %>%
        dplyr::pull(!!dplyr::sym(indexName)) %>%
        sample(initSize)
    } else {
      rdmDocs <- docs %>%
        dplyr::pull(!!dplyr::sym(indexName)) %>%
        sample(initSize)
    }

    ## Execute initial query
    docs <- query_label(
      docs, rdmDocs, n_class, labels,
      docName, indexName, labelsName,
      activeIter=count, maxIter=maxActive,
      handLabel
    )

    ## Get hand labeled index
    handLabeledIndex[[1]] <- rdmDocs

  } else {
    docs <- query_label(
      docs, initIndex, n_class, labels, docName,
      indexName, labelsName,
      activeIter=count,
      maxIter=maxActive, handLabel
    )
    handLabeledIndex[[1]] <- initIndex
  }
  ## Get vector of class labels
  classes <- get_classes(n_class)

  ## Create dfm containing all documents and terms in corpus
  if (is.null(dfm)) {
    dfm <- get_dfm(docs=docs, docName=docName, indexName=indexName, ...)
  }

  ## Split dfm into in-sample and out-sample sections
  dfms <- split_dfm(dfm, whichOutTest)

  ## Define which docs are to be used for in-sample validation
  count <- 1
  EM_docs[[count]] <- docs %>%
    dplyr::filter(!(!!dplyr::sym(indexName)) %in% whichOutTest)
  out_docs_static <- docs %>%
    dplyr::filter(!!dplyr::sym(indexName) %in% whichOutTest)

  ## configure max active steps appropriately given data size
  maxActive <- min(
    floor(nrow(EM_docs[[count]]) / maxQuery + 1 * (initSize/maxQuery)) - 1, maxActive
  )

  ## create a vector of n cluster values
  n_cluster_vec <- rep(n_cluster, maxActive + 1)

  ## create a vector of lambda values if decay is set to true
  if (lambda_decay) {
    decay <- function(rate, iters) {
      data <- c()
      for (i in 0:iters) {
        data[i + 1] <- 1 / (1 + rate * i)
      }
      return(data)
    }
    lambda_vec <- decay(rate = ld_rate, iters = maxActive)
  } else {
    lambda_vec <- rep(weight, maxActive + 1)
  }
  
  ## Main loop
  ## --------------------------------------------------------------------------
  repeat {

    ## Tunes lambda value at each iteraction, if enabled
    if (tune_lambda) {
      tuning_docs <- docs %>%
        dplyr::filter(!!dplyr::sym(indexName)
                      %in% handLabeledIndex[[count]])

      tune_out <- tune_lambda(
        docs = tuning_docs,
        n_clusters = n_cluster,
        lambdas = tune_lambda_range,
        prop_init = tune_lambda_prop_init,
        parallel = tune_lambda_parallel,
        k = tune_lambda_k,
        seed = seed
      )
      lambda_vec[count] <- tune_out$best$lambda
      n_cluster_vec[count] <- tune_out$best$n_cluster

    }

    ## Get EM_dfms and classification matrices
    EM_dfms <- split_dfm(
      dfm=dfms$second_section,
      splitIndex=handLabeledIndex[[count]]
    )

    ## Grab previous model parameters if necessary
    if (count == 1) {
      choose_NB_init <- TRUE
      prev_word_prob <- prev_class_prob <- NULL
    } else if (count > 1 & NB_init == FALSE) {
      choose_NB_init <- FALSE
      prev_word_prob <- output[[count - 1]]$eta
      prev_class_prob <- output[[count - 1]]$pi
    }

    ## Get initializing class prob matrix for EM algorithm.
    class_mtx <- get_class_matrix(
      EM_docs[[count]], n_class,
      handLabeledIndex[[count]],
      docName, indexName
    )

    ## Run EM Algorithm
    output[[count]] <- multiEM(
      .D_train = EM_dfms$first_section,
      .C_train=class_mtx,
      .D_test=EM_dfms$second_section,
      .n_class=n_class,
      .lambda=lambda_vec[count],
      .n_cluster=n_cluster_vec[count],
      .counter_on=counter_on,
      .active_iter=count,
      .maxactive_iter=maxActive,
      .fixed_words=fixed_words,
      .supervise=T,
      .choose_NB_init = choose_NB_init,
      .labeled_docs = handLabeledIndex[[count]],
      .prev_word_prob = prev_word_prob,
      .prev_class_prob = prev_class_prob
    )

    ## quick and dirty function to translate from multi
    if (n_cluster > 2) {
      EM_classlik <- matchCluster2Class(
        output[[count]]$classLik,
        count, n_cluster, n_class
      )
    } else {
      EM_classlik <- output[[count]]$classLik
    }

    ## Match EM output to document matrix by class
    ## NOTE: This only works when n_cluster > 2 when n_class = 2
    EM_docs[[count]] <- match_EM_to_docs(
      EM_docs[[count]],
      EMoutput=EM_classlik,
      classes, docName, indexName, labelsName
    )

    ## Match EM output to document matrix by cluster
    if (n_cluster > 2) {
      EM_docs[[count]] <- match_clusters_to_docs(
        docs = EM_docs[[count]],
        EMoutput = output[[count]]$classLik,
        indexName = indexName,
        n_cluster = n_cluster
      )
    }

    ## get out of sample prediction
    if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {
      out_prediction <- E_step(
        .D_test=dfms$first_section,
        .class_prob=output[[length(output)]]$pi,
        .word_prob=output[[length(output)]]$eta
      )
      if (n_cluster > 2) {
        EM_out_classlik <- matchCluster2Class(
          exp(out_prediction),
          count, n_cluster, n_class
        )
      } else {
        EM_out_classlik <- exp(out_prediction)
      }
      out_docs[[count]] <- match_EM_to_docs(
        out_docs_static,
        EMoutput=EM_out_classlik,
        classes, docName, indexName,
        labelsName
      )
    }

    ## check for convergence if using log-ratio sampling
    if (queryType == "log_ratio" & count > 1) {
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

    ## Check for maximum iterations reached
    if (count == maxActive) {
      if (counter_on) {
        message("\n Stopping: Reached maximum iterations.")
      }
      break
    }

    ## Label uncertain documents
    uncertainClass <- get_uncertain_docs(
      docs=EM_docs[[count]],
      bound=bound,
      maxQuery=maxQuery,
      indexName=indexName,
      handLabeledIndex=handLabeledIndex[[count]],
      forceList=forceList,
      queryType=queryType,
      ## dfm=dfms$second_section,
      ## mu=mu, tau=tau,
      ## EM_out=output[[count]],
      ## regions = regions,
      n_cluster = n_cluster
    )

    ## Check for convergence
    if (length(uncertainClass) == 0) {
      if (counter_on) {
        message("\n Stopping: Exhausted uncertain documents.")
      }
      break
    }

    ## Label uncertain documents if algorithm hasn't stopped
    EM_docs[[count + 1]] <- query_label(
      EM_docs[[count]], uncertainClass,
      n_class, labels, docName,
      indexName, labelsName,
      activeIter=count, maxIter=maxActive,
      handLabel
    )

    ## Update hand-labeled index
    handLabeledIndex[[count + 1]] <- c(handLabeledIndex[[count]],
                                       uncertainClass)

    ## Update counter
    count <- count + 1

  }

  ## End of main loop
  ## --------------------------------------------------------------------------

  ## Only export last active results if desired
  if (!exportAll) {
    handLabeledIndex <- handLabeledIndex[[length(handLabeledIndex)]]
    EM_docs <- EM_docs[[length(EM_docs)]]
    if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {
      out_docs <- out_docs[[length(out_docs)]]
    }
  }

  ## get vector of term sparsity values from dfm for exporting
  term_sparsity <- get_term_sparsity(dfm)

  ## add model metadata to output
  metadata <- list(
    trainSize = nrow(EM_docs[[1]]), corpusSize = nrow(docs),
    lambda = weight, initSize = initSize, maxQuery = maxQuery,
    queryType = queryType, pos_ratio = sum(docs[[labelsName]])/nrow(docs),
    mu = mu, tau = tau, lambda_decay = lambda_decay, ld_rate = ld_rate,
    n_cluster = n_cluster, NB_init = NB_init
  )

  ## export model output, depending on parameters choices
  if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {
    export_obj <- list(
      EMoutput=output, out_docs=out_docs, docs=EM_docs,
      handLabeledIndex=handLabeledIndex, term_sparsity=term_sparsity,
      lambda_vec=lambda_vec, n_cluster_vec=n_cluster_vec,
      metadata = metadata
    )
  } else {
    export_obj <- list(
      EMoutput=output, docs=EM_docs, handLabeledIndex=handLabeledIndex,
      term_sparsity=term_sparsity, lambda_vec=lambda_vec,
      n_cluster_vec = n_cluster_vec, metadata = metadata
    )
  }

  ## export only validation stats; useful if running many models
  if (exportAll & export_val_stats_only) {
    export_obj <- get_results_matrix(
      export_obj,
      out_sample = !(is.null(whichOutTest)),
      model_name = model_name
    )

  }

  ## return results
  return(export_obj)

}

#' @export
active_SVM <- function(docs, labels=c(0, 1), docName, indexName, labelsName=NULL, weight=1, n_class=2, initIndex=NULL,
                       handLabel=T, bound=0, maxActive=5, initSize=10, maxQuery=10, lazy_eval=F, forceList=F, counter_on=T, 
                       queryType="basic_entropy", whichOutTest=NULL, seed=NULL, fixed_words=NULL, dfm=NULL, exportAll=F, ...) {
  
  #' @title        Active Learning with SVM
  #'
  #' @param docs         [matrix]     Matrix of labeled and unlabeled documents, where each
  #'                                  row has index values and a nested Matrix of word tokens.
  #' @param labels       [vector]     Vector of character strings indicating classification
  #'                                  options for labeling.
  #' @param weight       [numeric]    Numeric value between 0 and 1. Used to weight unlabeled
  #'                                  documents.
  #' @param n_class      [numeric]    Number of classes to be considered.
  #' @param docName      [character]  Character string indicating the variable in 'docs'
  #'                                  that denotes the text of the documents to be classified.
  #' @param indexName    [character]  Character string indicating the variable in 'docs'
  #'                                  that denotes the index value of the document to be classified.
  #' @param labelsName   [character]  Character string indicating the variable in \code{docs}
  #'                                  that denotes the already known labels of the documents.
  #'                                  By default, value is set to \code{NULL}.
  #' @param initIndex    [vector]     Vector that indicates which documents to use to initialize the
  #'                                  algorithm. By default set to \code{NULL}, which causes a random
  #'                                  subset of the documents to be selected.
  #' @param handLabel    [logical]    Boolean logical value indicating whether to initiate user-input script.
  #'                                  If set to \code{FALSE}, and if \code{labelsName} is provided, the script
  #'                                  queries the document label directly from the column denoted by \code{labelsName}.
  #' @param bound        [numeric]    Minimum bound of entropy to call for additional labelling.
  #' @param maxActive    [numeric]    Value of maximum allowed active learning iterations.
  #' @param initSize     [numeric]    Value of maximum allowed iterations within the EM algorithm.
  #' @param maxQuery     [numeric]    Maximum number of documents queried in each EM iteration.
  #' @param lazy_eval    [logical]    If \code{lazy_eval == T}, convergence is measured by comparing changes
  #'                                  in log likelihood across model iterations rather than directly computing
  #'                                  maximand.
  #' @param forceList    [logical]    Switch indicating whether to force the filtering of documents with
  #'                                  no entropy. Set to \code{FALSE} by default.
  #' @param counter_on   [logical]    Switch indicating whether the progress of each sequence of the EM algorithm
  #'                                  is reported. By default set to \code{TRUE}.
  #' @param whichOutTest [vector]     Vector of document index labels used to identify documents to be used for
  #'                                  out of sample validation of the learned model. Set to \code{NULL} by default.
  #'                                  If a vector of labels is provided, the function outputs an additional argument
  #'                                  containing classification likelihoods for all documents identified by the vector.
  #' @param exportAll    [logical]    Switch indicating whether to export model predictions from each stage of the algorithm.
  #'                                  If true, the function exports a list of lists containing all predictions.
  #' @param queryType    [string]     String indicating which type of uncertainty sampling to use. Options are \code{"standard_entropy"}
  #'                                  or \code{"normalized_entropy"}, \code{"tiered_entropy"}, or \code{"tiered_entropy_weighted"}.
  #' @param quantileBreaks [vector]   Vector of break points to distinguish entropy zones. The first value is
  #'                                  the break point between the first and second tier, the second is the
  #'                                  break point between the second and third tier.
  #' @param sampleProps    [vector]   Vector of sampling proportions for each entropy zone. The first value is
  #'                                  the proportion of \code{maxQuery} to be sampled from the high entropy region,
  #'                                  the second value is the proportion to be sampled from the middle entropy region,
  #'                                  and the third value is the proportion to be sampled from the lowest entropy region.
  #' @param seed         [numeric]    Sets seed for model.
  #' @param fixed_words  [matrix]     Matrix of fixed words with class probabilities, where ncol is the number of classes.
  #' @param supervise    [logical]    T if supervised. F is unsupervised.
  #' @param dfm          [matrix]     Option to manually supply a dfm from quanteda.
  #' @param ...                       Additional parameters to pass to `get_dfm` and `EM()` and `get_uncertain_docs()`.
  #'
  #'
  #' @return             [list]       List containing labeled document matrix, prior weights, word likelihoods, and a vector
  #'                                  of user-labeled documents ids.
  
  
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Create containers
  handLabeledIndex <- list()
  EM_docs <- list()
  out_docs <- list()
  output <- list()
  
  # Set count at 0 for initialization
  count <- 0
  
  if (is.null(initIndex)) {
    # Randomly sample tweets
    if (!is.null(whichOutTest)) {
      rdmDocs <- docs %>%
        dplyr::filter(!(!!dplyr::sym(indexName)) %in% whichOutTest) %>%
        dplyr::sample_n(initSize)
    } else {
      rdmDocs <- dplyr::sample_n(docs, initSize)
    }
    
    # Execute initial query
    docs <- query_Label(docs, rdmDocs, n_class, labels, docName, indexName, labelsName,
                        activeIter=count, maxIter=maxActive, handLabel)
    
    # Get hand labeled index
    handLabeledIndex[[1]] <- get_index(rdmDocs, indexName)
  } else {
    handLabeledIndex[[1]] <- initIndex
  }
  
  # Get vector of class labels
  classes <- get_classes(n_class)
  
  # Set count at 1 for first run of the algorithm
  count <- 1
  
  # Create dfm containing all documents and terms in corpus
  if (is.null(dfm)) {
    dfm <- get_dfm(docs=docs, docName=docName, indexName=indexName, ...)
  }
  
  
  # Split dfm into in-sample and out-sample sections
  dfms <- split_dfm(dfm, whichOutTest)
  
  # Define which docs are to be used for in-sample validation
  EM_docs[[count]] <- docs %>% dplyr::filter(!(!!dplyr::sym(indexName)) %in% whichOutTest)
  out_docs_static <- docs %>% dplyr::filter(!!dplyr::sym(indexName) %in% whichOutTest)
  
  # Main loop
  repeat {
    
    # Get EM_dfms and classification matrices
    EM_dfms <- split_dfm(dfm=dfms$second_section, splitIndex=handLabeledIndex[[count]])
    class_mtx <- get_class_matrix(EM_docs[[count]], n_class, handLabeledIndex[[count]], docName, indexName)
    
    # -----SVM replaced EM in this section ------------------------------------------------------------------
    # Class 2 as positive class
    class_vec <- class_mtx[,2]
    
    
    # NOTE: Tuning function sometimes returns error. Manually impelment tuing?
    
    svm_tune <- tune(svm, train.x=EM_dfms$first_section, 
                     train.y=class_vec, 
                     kernel="linear", ranges=list(cost=sort(rexp(n=10, rate=50))),
                     tunecontrol=tune.control(sampling='cross',cross=5))
    
    # Fit SVM
    svm_fit <- svm(x=as.matrix(EM_dfms$first_section),
                   y=ifelse(as.matrix(class_mtx)[,2] == 1, 1, 0),
                   type='C-classification', kernel='linear',
                   cost=svm_tune$best.parameters, scale=F, probability=T)
    
    # predicted classes
    # have to merge with the labeled and predicted label! 
    pred_class <- predict(svm_fit, as.matrix(EM_dfms$second_section), probability=T)
    
    # predicted probabilities
    pred_prob <- rbind(class_mtx, attr(pred_class, "probabilities"))
    
    # match SVM output to the document matrix
    EM_docs[[count]] <- match_EM_to_docs(EM_docs[[count]], 
                                         EMoutput=pred_prob,
                                         classes, docName, indexName, labelsName)
    
    
    
    # get out of sample prediction
    if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {
      out_pred <- predict(svm_fit, as.matrix(dfms$first_section), probability=T)
      out_prediction <- attr(out_pred, "probabilities")
      out_docs[[count]] <- match_EM_to_docs(out_docs_static, 
                                            EMoutput=out_prediction,
                                            classes, docName, indexName, labelsName)
    }
    # --------------------------------------------------------------------------------------------------------
    
    if (counter_on) {
      Sys.sleep(.02)
    }
    
    # Check for maximum iterations reached
    if (count == maxActive + 1) {
      if (counter_on) {
        message("\n Stopping: Reached maximum iterations.")
      }
      break
    }
    
    
    
    # Get uncertain documents
    uncertainClass <- get_uncertain_docs(docs=EM_docs[[count]], bound=bound, maxQuery=maxQuery,
                                         indexName=indexName, handLabeledIndex=handLabeledIndex[[count]],
                                         forceList=forceList, queryType=queryType)
    
    
    # Check for convergence
    if (nrow(uncertainClass) == 0) {
      if (counter_on) {
        message("\n Stopping: Exhausted uncertain documents.")
      }
      break
    }
    
    # Query oracle for uncertain document labels
    EM_docs[[count + 1]] <- query_Label(EM_docs[[count]], uncertainClass, n_class, labels, docName,
                                        indexName, labelsName, activeIter=count, maxIter=maxActive, handLabel)
    
    # Update hand-labeled index
    handLabeledIndex[[count + 1]] <- c(handLabeledIndex[[count]], uncertainClass[[paste0(indexName)]])
    
    # Update counter
    count <- count + 1
  }
  
  if (exportAll == F) {
    handLabeledIndex <- handLabeledIndex[[length(handLabeledIndex)]]
    EM_docs <- EM_docs[[length(EM_docs)]]
    if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {
      out_docs <- out_docs[[length(out_docs)]]
    }
  }
  
  # get vector of term sparsity values from dfm for exporting
  term_sparsity <- get_term_sparsity(dfm)
  
  # return results
  return(
    if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {
      list(EMoutput=output, out_docs=out_docs, docs=EM_docs,
           handLabeledIndex=handLabeledIndex, term_sparsity=term_sparsity)
    } else {
      list(EMoutput=output, docs=EM_docs, handLabeledIndex=handLabeledIndex,
           term_sparsity=term_sparsity)
    }
  )
  
}



# test function with multi cluster EM
#' @export
active_multiEM <- function(docs, labels=c(0, 1), docName, indexName,
                           labelsName=NULL, weight=1, n_class=2, n_cluster=2,
                           initIndex=NULL, handLabel=T, bound=0, maxActive=5, 
                           initSize=10, maxQuery=10, lazy_eval=F, forceList=F, 
                           counter_on=T, queryType="basic_entropy", 
                           whichOutTest=NULL, seed=NULL, fixed_words=NULL, 
                           dfm=NULL, exportAll=F, log_ratio_threshold=0.001,
                           log_ratio_conv_type="maximand", mu=0.0001, tau=0.0001, ...) {
  
  #' @title        Active Learning EM Algorithm
  #' @description  Active learning for weighted-EM algorithm. After initial EM algorithm converges,
  #'               oracle is queried for labels to documents that the EM algorithm was most unsure of.
  #'               This process iterates until max iterations are reached, or there are no
  #'               documents in the window of uncertainty.
  #'
  #' @param docs         [matrix]     Matrix of labeled and unlabeled documents, where each
  #'                                  row has index values and a nested Matrix of word tokens.
  #' @param labels       [vector]     Vector of character strings indicating classification
  #'                                  options for labeling.
  #' @param weight       [numeric]    Numeric value between 0 and 1. Used to weight unlabeled
  #'                                  documents.
  #' @param n_class      [numeric]    Number of classes to be considered.
  #' @param docName      [character]  Character string indicating the variable in 'docs'
  #'                                  that denotes the text of the documents to be classified.
  #' @param indexName    [character]  Character string indicating the variable in 'docs'
  #'                                  that denotes the index value of the document to be classified.
  #' @param labelsName   [character]  Character string indicating the variable in \code{docs}
  #'                                  that denotes the already known labels of the documents.
  #'                                  By default, value is set to \code{NULL}.
  #' @param initIndex    [vector]     Vector that indicates which documents to use to initialize the
  #'                                  algorithm. By default set to \code{NULL}, which causes a random
  #'                                  subset of the documents to be selected.
  #' @param handLabel    [logical]    Boolean logical value indicating whether to initiate user-input script.
  #'                                  If set to \code{FALSE}, and if \code{labelsName} is provided, the script
  #'                                  queries the document label directly from the column denoted by \code{labelsName}.
  #' @param bound        [numeric]    Minimum bound of entropy to call for additional labelling.
  #' @param maxActive    [numeric]    Value of maximum allowed active learning iterations.
  #' @param initSize     [numeric]    Value of maximum allowed iterations within the EM algorithm.
  #' @param maxQuery     [numeric]    Maximum number of documents queried in each EM iteration.
  #' @param lazy_eval    [logical]    If \code{lazy_eval == T}, convergence is measured by comparing changes
  #'                                  in log likelihood across model iterations rather than directly computing
  #'                                  maximand.
  #' @param forceList    [logical]    Switch indicating whether to force the filtering of documents with
  #'                                  no entropy. Set to \code{FALSE} by default.
  #' @param counter_on   [logical]    Switch indicating whether the progress of each sequence of the EM algorithm
  #'                                  is reported. By default set to \code{TRUE}.
  #' @param whichOutTest [vector]     Vector of document index labels used to identify documents to be used for
  #'                                  out of sample validation of the learned model. Set to \code{NULL} by default.
  #'                                  If a vector of labels is provided, the function outputs an additional argument
  #'                                  containing classification likelihoods for all documents identified by the vector.
  #' @param exportAll    [logical]    Switch indicating whether to export model predictions from each stage of the algorithm.
  #'                                  If true, the function exports a list of lists containing all predictions.
  #' @param queryType    [string]     String indicating which type of uncertainty sampling to use. Options are \code{"standard_entropy"}
  #'                                  or \code{"normalized_entropy"}, \code{"tiered_entropy"}, or \code{"tiered_entropy_weighted"}.
  #' @param quantileBreaks [vector]   Vector of break points to distinguish entropy zones. The first value is
  #'                                  the break point between the first and second tier, the second is the
  #'                                  break point between the second and third tier.
  #' @param sampleProps    [vector]   Vector of sampling proportions for each entropy zone. The first value is
  #'                                  the proportion of \code{maxQuery} to be sampled from the high entropy region,
  #'                                  the second value is the proportion to be sampled from the middle entropy region,
  #'                                  and the third value is the proportion to be sampled from the lowest entropy region.
  #' @param seed         [numeric]    Sets seed for model.
  #' @param fixed_words  [matrix]     Matrix of fixed words with class probabilities, where ncol is the number of classes.
  #' @param supervise    [logical]    T if supervised. F is unsupervised.
  #' @param dfm          [matrix]     Option to manually supply a dfm from quanteda.
  #' @param log_ratio_threshold [numeric] Threshold at which convergence is declared when using `queryType="log_ratio"`.
  #' @param log_ratio_conv_type [string] If `queryType="log_ratio"`, this supplies the way that convergence is estimated.
  #'                                  Set to `maximand` by default.
  #' @param mu                        Parameters for error acceptance with `queryType=log_ratio`.
  #' @param tau                       Parameters for error acceptaance with `queryType=log_ratio`.
  #' @param ...                       Additional parameters to pass to `get_dfm` and `EM()` and `get_uncertain_docs()`.
  #'
  #'
  #' @return             [list]       List containing labeled document matrix, prior weights, word likelihoods, and a vector
  #'                                  of user-labeled documents ids.
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Create containers
  handLabeledIndex <- list()
  EM_docs <- list()
  out_docs <- list()
  output <- list()
  
  if (is.null(initIndex)) {
    # Randomly sample tweets
    if (!is.null(whichOutTest)) {
      rdmDocs <- docs %>%
        dplyr::filter(!(!!dplyr::sym(indexName)) %in% whichOutTest) %>%
        dplyr::sample_n(initSize)
    } else {
      rdmDocs <- dplyr::sample_n(docs, initSize)
    }
    
    # Execute initial query
    docs <- query_Label(docs, rdmDocs, n_class, labels, docName, indexName, labelsName,
                        activeIter=count, maxIter=maxActive, handLabel)
    
    # Get hand labeled index
    handLabeledIndex[[1]] <- get_index(rdmDocs, indexName)
  } else {
    handLabeledIndex[[1]] <- initIndex
  }
  
  # Get vector of class labels
  classes <- get_classes(n_cluster)
  
  # Create dfm containing all documents and terms in corpus
  if (is.null(dfm)) {
    dfm <- get_dfm(docs=docs, docName=docName, indexName=indexName, ...)
  }
  
  # Split dfm into in-sample and out-sample sections
  dfms <- split_dfm(dfm, whichOutTest)
  

  # Define which docs are to be used for in-sample validation
  count <- 1
  EM_docs[[count]] <- docs %>% dplyr::filter(!(!!dplyr::sym(indexName)) %in% whichOutTest)
  out_docs_static <- docs %>% dplyr::filter(!!dplyr::sym(indexName) %in% whichOutTest)
  
  # Main loop
  repeat {
    
    # Get EM_dfms and classification matrices
    EM_dfms <- split_dfm(dfm=dfms$second_section, 
                         splitIndex=handLabeledIndex[[count]])
    class_mtx <- get_class_matrix(EM_docs[[count]], n_class, 
                                  handLabeledIndex[[count]], 
                                  docName, indexName)
    
    # Run EM Algorithm
    output[[count]] <- multiEM(.D_train=EM_dfms$first_section, 
                               .C_train=class_mtx, 
                               .D_test=EM_dfms$second_section,
                               .n_class=n_class,
                               .n_cluster=n_cluster,
                               .lambda=weight, 
                               .counter_on=counter_on, 
                               .active_iter=count,
                               .maxactive_iter=maxActive, 
                               .fixed_words=fixed_words)
    
    
    # Match EM output to document matrix
    EM_docs[[count]] <- match_EM_to_docs(EM_docs[[count]], 
                                         EMoutput=output[[count]]$classLik, 
                                         classes, docName, indexName, labelsName)
    
    # get out of sample prediction
    if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {


      out_prediction <- E_step(.D_test=dfms$first_section,
                               .class_prob=output[[length(output)]]$pi,
                               .word_prob=output[[length(output)]]$eta)

      out_docs[[count]] <- match_EM_to_docs(out_docs_static, 
                                            EMoutput=exp(out_prediction), 
                                            classes, docName, indexName, 
                                            labelsName)
    }
    
    if (counter_on) {
      Sys.sleep(.02)
    }
    
    # check for convergence if using log-ratio sampling
    if (queryType == "log_ratio" & count > 1) {
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
    
    # Check for maximum iterations reached
    if (count == maxActive) {
      if (counter_on) {
        message("\n Stopping: Reached maximum iterations.")
      }
      break
    }
    
    # Label uncertain documents
    uncertainClass <- get_uncertain_docs(docs=EM_docs[[count]],
                                         bound=bound,
                                         maxQuery=maxQuery,
                                         indexName=indexName,
                                         handLabeledIndex=handLabeledIndex[[count]],
                                         forceList=forceList,
                                         queryType=queryType,
                                         dfm=dfms$second_section,
                                         mu=mu, tau=tau,
                                         EM_out=output[[count]],
                                         seed=seed)
    # Check for convergence
    if (nrow(uncertainClass) == 0) {
      if (counter_on) {
        message("\n Stopping: Exhausted uncertain documents.")
      }
      break
    }

    # Label uncertain documents if algorithm hasn't stopped
    EM_docs[[count + 1]] <- query_Label(EM_docs[[count]], uncertainClass, 
                                        n_class, labels, docName,
                                        indexName, labelsName, 
                                        activeIter=count, maxIter=maxActive,
                                        handLabel)
    
    # Update hand-labeled index
    handLabeledIndex[[count + 1]] <- c(handLabeledIndex[[count]], 
                                       uncertainClass[[paste0(indexName)]])
    
    # Update counter
    count <- count + 1
    
  }
  
  if (exportAll == F) {
    
    handLabeledIndex <- handLabeledIndex[[length(handLabeledIndex)]]
    EM_docs <- EM_docs[[length(EM_docs)]]
    if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {
      out_docs <- out_docs[[length(out_docs)]]
    }
  }
  
  # get vector of term sparsity values from dfm for exporting
  term_sparsity <- get_term_sparsity(dfm)
  
  # return results
  return(
    
    if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {
      
      list(EMoutput=output, out_docs=out_docs, docs=EM_docs,
           handLabeledIndex=handLabeledIndex, term_sparsity=term_sparsity)
    } else {
      
      list(EMoutput=output, docs=EM_docs, handLabeledIndex=handLabeledIndex,
           term_sparsity=term_sparsity)
    }
  )
  
}

