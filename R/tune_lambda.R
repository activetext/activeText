##################################################
## Project: active
## Script purpose: Functions for Tuning Lambda
## Date: 2020/2/23
## Author: Mitchell Bosley
##################################################

#' @export
tune_lambda_helper <- function(i, docs,
                               docName,
                               indexName,
                               n_cluster,
                               labelsName,
                               initSize,
                               initIndex,
                               lambda,
                               maxQuery,
                               whichOutTest,
                               maxActive,
                               seed,
                               verbose = TRUE) {

#' @title Helper Function for `Tune_Lambda()`
#' @description helper function for tune_lambda()
#'
#' @param i Iteration number.
#' @param docs Dataframe for courpus.
#' @param docName String; name of variable for document text.
#' @param indexName String; name of variable for corpus index.
#' @param n_cluster Int; number of clusters.
#' @param labelsName String; name of variable for document labels.
#' @param initSize Int; number of docs to initialize model with.
#' @param initIndex vec; Value to initialize model with.
#' @param lambda float; lambda value for weighting.
#' @param maxQuery int; max number of documents ot query per iteration.
#' @param whichOutTest vec; Index values to use for out-of-sample
#' validation.
#' @param seed int; Seed value for randomization.
#' @param verbose logical; Denotes whether model progress is printed.
#'
#' @return Statistics summary for active_EM() iteration.


  if (verbose) {
    print(paste0('model is running with lambda = ', lambda,
                 ' with n_clusters = ', n_cluster))
  }
  res <- active_EM(docs=docs,
                   docName=docName,
                   n_cluster=n_cluster,
                   indexName=indexName,
                   labelsName=labelsName,
                   weight=lambda,
                   handLabel=F,
                   initIndex=initIndex,
                   seed=seed,
                   forceList=T,
                   maxQuery=maxQuery,
                   counter_on=FALSE,
                   whichOutTest=whichOutTest,
                   maxActive=maxActive)
  res_em <- res$EMoutput

  # estimated pi
  pi_est <- exp(res_em[[length(res_em)]]$pi)[2]

  # confusion matrix
  res_docs <- res$out_docs
  res_docs$pred_class <- ifelse(res_docs$Class_2 >= 0.5, 1, 0)
  cf <- get_conf_matrix(docs = res$out_docs,
                      labelsName = "label", indexName = "id")
  TP <- cf[2,2]
  TN <- cf[1,1]
  FP <- cf[1,2]
  FN <- cf[2,1]
  F1 <- 2*TP / (2*TP + FN + FP)
  if (is.nan(F1)) {
    F1 <- 0
  }
  
  # check the proportion of pos/neg in the labeled documents
  res_hand <- res$handLabeledIndex
  pi_labeled <- mean(res$docs[res$docs$id %in% res_hand,]$label)

  return(c(lambda, n_cluster, TP, TN, FP, FN, F1, pi_est, pi_labeled))
}

find_best_ncluster_lambda <- function(result){
  best <- as_tibble(result) %>%
  group_by(n_cluster, lambda) %>%
  summarize(mean_F1 = mean(F1)) %>%
  top_n(n=1, wt=mean_F1) %>%
  select(n_cluster, lambda)
  return(best)
}

#' @export
tune_lambda <- function(docs,
                        docName = "text",
                        indexName = "id",
                        labelsName ='label',
                        n_clusters = c(2),
                        maxQuery = 1,
                        active = FALSE,
                        lambdas, seed,
                        prop_init = 0.1,
                        parallel = TRUE,
                        k = 5) {

#' @title Tune Lambda Weight and Number of Clusters
#' @description Experiment function for tuning lambda.
#'
#' @param docs Dataframe for courpus.
#' @param docName String; name of variable for document text.
#' @param indexName String; name of variable for corpus index.
#' @param labelsName String; name of variable for document labels.
#' @param n_clusters vec; sequence of number of clusters to be tested.
#' @param lambdas vec; sequence of lambda values to test.
#' @param maxQuery int; max number of documents ot query per iteration.
#' @param seed int; Seed value for randomization.
#' @param prop_init; float; Proportion of docs in each training set
#' to directly label; the rest are used by the EM.
#' @param parallel logical; Whether or not to parallelize the tuning process.
#' @param k int; number for cross-validation.
#'
#' @return Results from experiment.

  # storage
  l <- length(lambdas) * length(n_clusters) * k
  result <- data.frame(
    lambda = rep(0, l),
    n_cluster = rep(0, l),
    pi_est = rep(0, l),
    TP = rep(0, l),
    TN = rep(0, l),
    FP = rep(0, l),
    FN = rep(0, l),
    F1 = rep(0, l),
    pi_labeled = rep(0, l)
  )

  # cross validation exercise
  # split data up into k chunks
  docs_split_ls <- docs %>%
      group_by((row_number() - 1) %/% (n() / k)) %>%
      nest %>% pull(data)

  # get possible combinations of testing sets
  # NOTE: make sure total docs is divisible by k evenly
  test_sets <- combn(k, k - 1)
  train_cont <- test_cont <- vector(mode = "list", length = k)
  for (col in 1:ncol(test_sets)) {
    for (row in 1:nrow(test_sets)) {
      val <- test_sets[row, col]
      train_cont[[col]] <- rbind(train_cont[[col]], docs_split_ls[[val]])
    }
    test_val <- which(!(1:k %in% test_sets[, col]))
    test_cont[[col]] <- docs_split_ls[[test_val]]
  }

  # get proper list sizes for loop
  obj_train <- obj_test <- list()
  for (n_cluster in n_clusters) {
    for (lambda in lambdas) {
        obj_train <- append(obj_train, train_cont)
        obj_test <- append(obj_test, test_cont)
    }
  }
  train_cont <- obj_train
  test_cont <- obj_test

  # define armount of labeled data to use
  init_num <- nrow(train_cont[[1]]) * prop_init

  # define number of active iterations NOTE: Buggy
  if (active == TRUE) {
    active_num <- (nrow(train_cont[[1]]) - init_num) / maxQuery
  } else {
    active_num <- 0
  }

  # define lambda values
  lambda_seq <- c()
  for (lambda in lambdas) {
    lambda_seq <- c(lambda_seq, rep(lambda, k))
  }
  obj <- c()
  for (n_cluster in n_clusters) {
    obj <- c(obj, lambda_seq)
  }
  lambda_seq <- obj

  # define n_cluster values
  n_cluster_seq <- c()
  for (n_cluster in n_clusters) {
    for (lambda in lambdas) {
      n_cluster_seq <- c(n_cluster_seq, rep(n_cluster, k))
    }
  }

  # serial
  if (!parallel) {
    result <- foreach::foreach(i=1:(k * length(lambdas) * length(n_clusters)),
                      .combine = rbind, .export=c("tune_lambda_helper"),
                      .packages = c("dplyr", "activeEMtext", "parallel",
                                    "quanteda")) %do% {

                                      tune_lambda_helper(i = i,
                                                         docs = rbind(test_cont[[i]],
                                                                      train_cont[[i]]),
                                                         docName = docName,
                                                         indexName = indexName,
                                                         labelsName = labelsName,
                                                         maxQuery = maxQuery,
                                                         maxActive = active_num,
                                                         lambda = lambda_seq[i],
                                                         seed = seed,
                                                         n_cluster = n_cluster_seq[i],
                                                         initIndex = train_cont[[i]]$id[1:init_num],
                                                         whichOutTest = test_cont[[i]]$id)

                                    }
  } else {
    cl <- parallel::makeCluster(parallel::detectCores() - 1, outfile = "")
    doParallel::registerDoParallel(cl)
    result <-
      foreach::foreach(
                 i=1:(k * length(lambdas) * length(n_clusters)),
                 .combine = rbind, .export = c("tune_lambda_helper"),
                 .packages = c("dplyr", "activeEMtext", "parallel",
                               "quanteda")) %dopar% {
                                 tune_lambda_helper(i = i,
                                                    docs = rbind(test_cont[[i]],
                                                                 train_cont[[i]]),
                                                    docName = docName,
                                                    indexName = indexName,
                                                    labelsName = labelsName,
                                                    maxQuery = maxQuery,
                                                    maxActive = active_num,
                                                    lambda = lambda_seq[i],
                                                    seed = seed,
                                                    n_cluster = n_cluster_seq[i],
                                                    initIndex = train_cont[[i]]$id[1:init_num],
                                                    whichOutTest = test_cont[[i]]$id) }
    parallel::stopCluster(cl)
  }

  #' parallel
  colnames(result) <- c('lambda', 'n_cluster', 'TP', 'TN',
                        'FP','FN', 'F1', 'pi_est', 'pi_labeled')
  rownames(result) <- NULL
  #' select best # of clusters and lambda in terms of the mean F1
  best <- find_best_ncluster_lambda(tibble::as_tibble(result))

  return(list(result = as_tibble(result), pi_pop = mean(docs$label),
              best = best))
}

#' @export
find_best_ncluster_lambda <- function(out, show_all = FALSE){
  if (!show_all) {
    best <- out %>%
      dplyr::group_by(n_cluster, lambda) %>%
      dplyr::summarize(mean_F1 = mean(F1)) %>%
      top_n(n = 1, wt = mean_F1) %>%
      dplyr::slice(1) %>%
      dplyr::select(n_cluster, lambda)
  } else {
    best <- out %>%
      dplyr::group_by(n_cluster, lambda) %>%
      dplyr::summarize(mean_F1 = mean(F1)) %>%
      top_n(n = 1, wt = mean_F1)
  }
  return(best)
}

#' @export
get_tune_lambda_fig <- function(out, k) {
#' @title Gets Tune Lambda Figure
#' @description Gets visualization from output of `tune_lambda` function.
  plot <- out$result %>%
    dplyr::mutate(n_cluster = as.factor(n_cluster)) %>%
    dplyr::group_by(n_cluster, lambda) %>%
    dplyr::summarize(mean_F1 = mean(F1),
              var_F1 = var(F1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(upper_ci = mean_F1 + 1.96 * var_F1 / sqrt(k),
          lower_ci = mean_F1 - 1.96 * var_F1 / sqrt(k)) %>%
    ggplot2::ggplot(ggplot2::aes(x = lambda, y = mean_F1,
                                 fill = n_cluster, color = n_cluster)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                color = "grey70", alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::ggtitle("Finding Optimal Lambda with Cross-Validation")
  return(plot)
}

#' @export
get_cross_validate_key <- function(docs, k) {
  rsample::vfold_cv()
  docs_split_ls <- docs %>%
      group_by((row_number() - 1) %/% (n() / k)) %>%
      nest %>% pull(data)
}
