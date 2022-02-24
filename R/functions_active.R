###############################################################################
## TITLE : Active EM Dev
## PROJECT : Active
## NAME : Mitchell Bosley
## DATE : 2020-07-16
###############################################################################

#' @export
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
#' @param lambda       [numeric]    Numeric value between 0 and 1. Used to weight unlabeled
#'                                  documents.
#' @param n_class      [numeric]    Number of classes to be considered.
#' @param doc_name      [character]  Character string indicating the variable in 'docs'
#'                                  that denotes the text of the documents to be classified.
#' @param index_name    [character]  Character string indicating the variable in 'docs'
#'                                  that denotes the index value of the document to be classified.
#' @param labels_name   [character]  Character string indicating the variable in \code{docs}
#'                                  that denotes the already known labels of the documents.
#'                                  By default, value is set to \code{NULL}.
#' @param initIndex    [vector]     Vector that indicates which documents to use to initialize the
#'                                  algorithm. By default set to \code{NULL}, which causes a random
#'                                  subset of the documents to be selected.
#' @param handlabel    [logical]    Boolean logical value indicating whether to initiate user-input script.
#'                                  If set to \code{FALSE}, and if \code{labels_name} is provided, the script
#'                                  queries the document label directly from the column denoted by \code{labels_name}.
#' @param bound        [numeric]    Minimum bound of entropy to call for additional labelling.
#' @param max_active    [numeric]    Value of maximum allowed active learning iterations.
#' @param init_size     [numeric]    Value of maximum allowed iterations within the EM algorithm.
#' @param max_query     [numeric]    Maximum number of documents queried in each EM iteration.
#' @param lazy_eval    [logical]    If \code{lazy_eval == T}, convergence is measured by comparing changes
#'                                  in log likelihood across model iterations rather than directly computing
#'                                  maximand.
#' @param force_list   [logical]    Switch indicating whether to force the filtering of documents with
#'                                  no entropy. Set to \code{FALSE} by default.
#' @param counter_on   [logical]    Switch indicating whether the progress of each sequence of the EM algorithm
#'                                  is reported. By default set to \code{TRUE}.
#' @param which_out_test [vector]     Vector of document index labels used to identify documents to be used for
#'                                  out of sample validation of the learned model. Set to \code{NULL} by default.
#'                                  If a vector of labels is provided, the function outputs an additional argument
#'                                  containing classification likelihoods for all documents identified by the vector.
#' @param export_all    [logical]    Switch indicating whether to export model predictions from each stage of the algorithm.
#' @param export_all_em [logical]    Switch indicating whether to export model
#'                                  If true, the function exports a list of lists containing all predictions.
#' @param export_val_stats_only     Boolean, indicating whether to export validation stats only from model runs.
#' @param query_type    [string]     String indicating which type of uncertainty sampling to use. Options are \code{"standard_entropy"}
#'                                  or \code{"normalized_entropy"}, \code{"tiered_entropy"}, or \code{"tiered_entropy_weighted"}.
#' @param quantileBreaks [vector]   Vector of break points to distinguish entropy zones. The first value is
#'                                  the break point between the first and second tier, the second is the
#'                                  break point between the second and third tier.
#' @param sampleProps    [vector]   Vector of sampling proportions for each entropy zone. The first value is
#'                                  the proportion of \code{max_query} to be sampled from the high entropy region,
#'                                  the second value is the proportion to be sampled from the middle entropy region,
#'                                  and the third value is the proportion to be sampled from the lowest entropy region.
#' @param seed         [numeric]    Sets seed for model.
#' @param fixed_words  [matrix]     Matrix of fixed words with class probabilities, where ncol is the number of classes.
#' @param supervise    [logical]    T if supervised. F is unsupervised.
#' @param dfms          [matrix]     Option to manually supply a dfm from quanteda.
#' @param log_ratio_threshold [numeric] Threshold at which convergence is declared when using `query_type="log_ratio"`.
#' @param log_ratio_conv_type [string] If `query_type="log_ratio"`, this supplies the way that convergence is estimated.
#'                                  Set to `maximand` by default.
#' @param mu                        Parameters for error acceptance with `query_type=log_ratio`.
#' @param tau                       Parameters for error acceptaance with `query_type=log_ratio`.
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
#' @param agg_type [string] Indicating how to aggregate model predictions.
#' @param n_cluster_collapse_type [string] Indicates how to collapse multiple clusters into binary class.
#' By default, set to "simple", which takes the negative class probablity as the difference between the
#' positive class probability and 1. Can also be set to "max_neg", which calculates the normalized ratio of
#' positive cluster to the largest negative cluster.
#' @param beta [numeric] prior parameter for eta
#' @param active_eta_query [boolean] Indicates whether to query oracle for eta tuning.
#' @param cont_metadata_varnames Vector of continuous metadata varnames
#' @param binary_metadata_varnames Vector of binary metadata varnames
#' @param contextual_metadata_varnames Vector of contextual metadata varnames
#
#' @param ... Additional parameters to pass to `get_dfm` and `EM()` and `get_uncertain_docs()`.
#'
#'
#' @return             [list]       List containing labeled document matrix, prior weights, word likelihoods, and a vector
#'                                  of user-labeled documents ids.
active_EM_dev <- function(docs,
                          labels=c(0, 1),
                          doc_name = "text",
                          index_name = "id",
                          labels_name = NULL,
                          lambda = 1,
                          n_class = 2,
                          n_cluster = 2,
                          init_index = NULL,
                          handlabel = T,
                          bound = 0,
                          max_active = 5,
                          init_size = 10,
                          max_query = 10,
                          lazy_eval = FALSE,
                          force_list = FALSE,
                          counter_on = TRUE,
                          query_type = "basic_entropy",
                          which_out_test = NULL,
                          seed = NA,
                          fixed_words = NULL,
                          dfms = NULL,
                          export_all_em = FALSE,
                          export_all = FALSE,
                          log_ratio_threshold = 0.001,
                          log_ratio_conv_type = "maximand",
                          mu = 0.0001,
                          tau = 0.0001,
                          regions = "both",
                          lambda_decay = FALSE,
                          ld_rate = 0.2,
                          tune_lambda = FALSE,
                          tune_lambda_prop_init = 0.1,
                          tune_lambda_range = seq(0, 1, 0.1),
                          tune_lambda_k = 10,
                          tune_lambda_parallel = TRUE,
                          NB_init = TRUE,
                          export_val_stats_only = FALSE,
                          model_name = "Model",
                          agg_type = "best",
                          n_cluster_collapse_type = "simple",
                          beta = NA,
                          active_eta_query = FALSE,
                          keywords_list = NA,
                          keywords_scheme = NA,
                          true_eta = NA,
                          gamma = NA,
                          validation_mode = FALSE,
                          cont_metadata_varnames = NA,
                          binary_metadata_varnames = NA,
                          contextual_varnames = NA,
                          mc_iter = NA,
                          ...) {


  ## Messages
  ## --------------------------------------------------------------------------
  active_initial_messages(n_cluster, query_type)

  ## Setup
  ## --------------------------------------------------------------------------
  if (!is.na(seed)) {
    set.seed(seed)
  }

  ## initialize count
  count <- 1

  ## get initial documents
  hand_labeled_index <- to_label_ind <- get_initial_documents(
    docs, init_index, index_name, init_size, which_out_test
  )

  ## label initial documents
  docs <- query_label(
    docs, to_label_ind, n_class, labels,
    doc_name, index_name, labels_name,
    active_iter = 0, maxIter = max_active,
    handlabel = handlabel,
    metadata_vars = c(contextual_varnames,
                      cont_metadata_varnames,
                      binary_metadata_varnames)
  )

  ## Create dfm containing all documents and terms in corpus
  if (is.null(dfms)) {
    dfms <- list()
    dfms[[1]] <- get_dfm(docs=docs, doc_name=doc_name, index_name=index_name, ...)
  } else if (!is.list(dfms)) {
    dfms <- list(dfms)
  }

  ## initialize beta table
  if (is.na(beta)) {
    beta_tbl <- initialize_beta_tbl(dfms, n_class, keywords_list, gamma)
  }

  ## Define which docs are to be used for in-sample validation
  if (!is.null(which_out_test)) {
    in_docs <- docs %>%
      dplyr::filter(!(!!dplyr::sym(index_name)) %in% which_out_test)
    out_docs <- docs %>%
      dplyr::filter(!!dplyr::sym(index_name) %in% which_out_test)
  } else {
    in_docs <- docs
    out_docs <- NULL
  }

  ## configure max active steps appropriately given data size
  max_active <- min(
    floor(nrow(in_docs[[count]]) / max_query + 1 * (init_size/max_query)) - 1,
    max_active
  )

  ## create a vector of n cluster and lambda values
  n_cluster_vec <- rep(n_cluster, max_active + 1)
  lambda_vec <- generate_lambda_vec(
    lambda_decay, lambda,
    rate, iters, max_active
  )

  ## create keywords list export object
  keywords_export <- list()
  keywords_export[[1]] <- keywords_list

  ## adds metadata and results table if in validation mode
  if (validation_mode) {
    ## add model metadata to output
    metadata <- list(
      train_size = nrow(in_docs), corpus_size = nrow(docs),
      lambda = lambda, init_size = init_size, max_query = max_query,
      query_type = query_type, pos_ratio = sum(docs[[labels_name]])/nrow(docs),
      mu = mu, tau = tau, lambda_decay = lambda_decay, ld_rate = ld_rate,
      n_cluster = n_cluster, NB_init = NB_init, n_dfms = length(dfms),
      agg_type = agg_type, n_cluster_collapse_type = n_cluster_collapse_type,
      gamma = gamma, mc_iter = mc_iter, seed = seed,
      keywords_scheme = keywords_scheme
    )

    ## get results summary table
    res_tbl <- gen_results_tbl(
      !is.null(which_out_test),
      metadata, max_active + 1, model_name
    )
  }

  ## creates container for exporting results
  if (export_all) {
    em_res_list <- list()
  }

  ## Main loop
  ## --------------------------------------------------------------------------

  # start global timer here
  repeat {
    # measure time
    start_time <- Sys.time()

    ## Tunes lambda value at each iteraction, if enabled
    if (tune_lambda) {
      tl_out <- tune_lambda_in_active(
        docs, index_name, hand_labeled_index, n_cluster,
        tune_lambda_range, tune_lambda_prop_init,
        tune_lambda_parallel, tune_lambda_k, seed
      )
      lambda_vec[count] <- tl_out$best$lambda
      n_cluster_vec[count] <- tl_out$best$n_cluster
    }

    ## Get initializing class prob matrix for EM algorithm.
    class_mtx <- get_class_matrix(
      in_docs, n_class, hand_labeled_index,
      doc_name, index_name
    )

    ## actively update eta
    if (count > 1 && active_eta_query){
      if (length(prev_model_output) > 1) {
        warning("using multiple dfms not supported for active word prob updating.")
      }
      active_eta_out <- active_eta_update(
        beta_tbl = beta_tbl,
        prev_model_output = prev_model_output,
        n_class = n_class, n_query = 10, gamma = gamma,
        handlabel = handlabel, true_eta = true_eta,
        dfm = dfms[[1]], keywords_list = keywords_list,
        update_scheme = "update_list"
      )
      beta_tbl <- active_eta_out$beta_tbl
      keywords_list <- active_eta_out$keywords_list
      keywords_export[[count]] <- keywords_list
    }

    ## run EM until convergence for each DFM specification
    i <- 1
    model_outputs <- list()

    for (dfm in dfms) {

      choose_NB_init <- ifelse(count > 1 && NB_init, TRUE, FALSE)

      model_outputs[[i]] <- multiEM(
        .D_train = dfm[rownames(dfm) %in% hand_labeled_index, ],
        .C_train = class_mtx,
        .D_test = dfm[!(rownames(dfm) %in% which_out_test) &
                      !(rownames(dfm) %in% hand_labeled_index), ],
        .D_test_out = dfm[rownames(dfm) %in% which_out_test, ],
        .n_class = n_class, .lambda = lambda_vec[count],
        .n_cluster = n_cluster_vec[count], .counter_on = counter_on,
        .active_iter = count, .maxactive_iter = max_active,
        .fixed_words = fixed_words, .supervise = TRUE,
        .choose_NB_init = ifelse(count == 1 | choose_NB_init, TRUE, FALSE),
        .labeled_docs = hand_labeled_index,
        .prev_word_prob = `if`(choose_NB_init, NA, prev_model_output[[i]]$eta),
        .prev_class_prob = `if`(choose_NB_init, NA, prev_model_output[[i]]$pi),
        .prev_mu = `if`(
          choose_NB_init | !is.na(cont_metadata_varnames),
          NA, prev_model_output[[i]]$mu
        ),
        .prev_psi = `if`(
          choose_NB_init | !is.na(binary_metadata_varnames),
          NA, prev_model_output[[i]]$psi
        ),
        .beta = beta_tbl,
        .cont_metadata_varnames = cont_metadata_varnames,
        .binary_metadata_varnames = binary_metadata_varnames
      )

      i <- i + 1
    }

    ## export em results if needed
    if (export_all) {
      em_res_list[[count]] <- model_outputs
    }

    ## aggregate model parameters from multiple dfm versions
    agg_output <- aggregate_model_predictions(
      model_out_to_tbl(model_outputs),
      dfms = dfms,
      models = model_outputs,
      val_data = in_docs %>%
        filter(!!sym(index_name) %in% hand_labeled_index),
      n_cluster = n_cluster,
      agg_type = agg_type,
      n_cluster_collapse_type = n_cluster_collapse_type
    )

    in_docs <- update_docs(in_docs, agg_output$in_agg)
    out_docs <- update_docs(out_docs, agg_output$out_agg)

    ## update results table
    if (validation_mode) {
      res_tbl <- update_results(
        !is.null(which_out_test),
        res_tbl, in_docs, out_docs,
        hand_labeled_index, count,
        Sys.time() - start_time
      )
    }

    ## check for convergence if using log-ratio sampling
    if (query_type == "log_ratio" && count > 1) {
      check_lr_convergence(
        output, count, log_ratio_threshold, log_ratio_conv_type
      )
    }

    ## Check for maximum iterations reached
    if (count == max_active + 1) {
      if (counter_on) {
        message("\n Stopping: Reached maximum iterations.")
      }
      break
    }

    ## Label uncertain documents
    uncertain_class <- get_uncertain_docs(
      docs = in_docs, bound = bound, max_query = max_query,
      index_name = index_name, hand_labeled_index = hand_labeled_index,
      force_list = force_list, query_type = query_type, n_cluster = n_cluster
    )

    ## Check for convergence
    if (length(uncertain_class) == 0) {
      if (counter_on) {
        message("\n Stopping: Exhausted uncertain documents.")
      }
      break
    }

    ## Label uncertain documents if algorithm hasn't stopped
    in_docs <- query_label(
      in_docs, uncertain_class,
      n_class, labels, doc_name,
      index_name, labels_name,
      active_iter = count, maxIter = max_active,
      handlabel = handlabel,
      metadata_vars = c(contextual_varnames,
                        cont_metadata_varnames,
                        binary_metadata_varnames)
    )

    ## Update hand-labeled index
    hand_labeled_index <- c(hand_labeled_index, uncertain_class)

    ## impute previous model output
    prev_model_output <- model_outputs

    ## Update counter
    count <- count + 1

  }

  ## End of main loop
  ## --------------------------------------------------------------------------

  ## export model output, depending on parameters choices
  export_obj <- list(
    EMoutput = `if`(export_all, em_res_list, model_outputs),
    docs = in_docs,
    hand_labeled_index = hand_labeled_index,
    lambda_vec = lambda_vec,
    n_cluster_vec = n_cluster_vec
  )
  if (validation_mode) {
    metadata[["neg_keywords"]] <-  `if`(
      is.na(keywords_list[[1]]) && length(keywords_list[[1]]) == 1,
      NA, paste(keywords_list[[1]], collapse = " ")
    )
    metadata[["pos_keywords"]] <- `if`(
      is.na(keywords_list[[2]]) && length(keywords_list[[1]]) == 1,
      NA, paste(keywords_list[[2]], collapse = " ")
    )
    export_obj[["metadata"]] <- metadata
  }
  if (!is.null(which_out_test) && length(which_out_test) != 0) {
    export_obj[["out_docs"]] <- out_docs
  }
  export_obj[["keywords_export"]] <- keywords_export

  ## measure time (end)
  time_elapsed <- Sys.time() - start_time

  ## export only validation stats; useful if running many models
  if (validation_mode) {
    if (export_val_stats_only) {
      return(res_tbl)
    } else {
      export_obj[["res_tbl"]] <- res_tbl
      return(export_obj)
    }
  } else {
    ## return all results
    return(export_obj)
  }
}
