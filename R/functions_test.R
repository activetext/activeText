# test script for varying lambda
#' @export
active_EM_varylambda2 <- function(docs, labels=c(0, 1), docName, indexName,
                                labelsName=NULL, weight=1, n_class=2,
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
    classes <- get_classes(n_class)
    
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
        
        #use difference lambda depending on the proportion of documents labeles
        label_ratio <- nrow(EM_dfms$first_section)/nrow(dfms$second_section)
        # try with the simple decreasing function
        #vary_lambda <- 1 - label_ratio
        # or uniformly decreasing
        vary_lambda <- seq(from=1, to=0, length.out=maxActive)[count]
        
        lambda_seq <- seq(0, 1, .1) 
        
        # Run EM Algorithm
        lambda_trials <- list()
        for (i in 1:length(lambda_seq)) {
            
            lambda_trials[[i]] <- EM(.D_train=EM_dfms$first_section, 
                                   .C_train=class_mtx, 
                                   .D_test=EM_dfms$second_section,
                                   .n_class=n_class, .lambda=lambda_seq[i], 
                                   .lazy_eval=lazy_eval, 
                                   .counter_on=counter_on, 
                                   .active_iter=count,
                                   .maxactive_iter=maxActive, 
                                   .fixed_words=fixed_words, 
                                   .supervise=T)
        }
        
        
        cfs <- list()
        F1_scores <- c()
        lambda_trials_output <- list()
        for (i in 1:length(lambda_trials)) {
            lambda_trials_output[[i]] <- match_EM_to_docs(EM_docs[[count]], 
                                                          EMoutput=lambda_trials[[i]]$classLik, 
                                                          classes, docName, indexName, labelsName)
            cfs[[i]] <- get_conf_matrix(docs=lambda_trials_output[[i]], labelsName=labelsName,
                                       indexName=indexName, labeledIndex=handLabeledIndex[[count]],
                                       n_class=n_class, levels=c(0, 1), useLabeled=F)
            
            F1_scores[i] <- get_F1_binary(cfs[[i]])
        }
        
        which_max <- which(F1_scores %in% max(F1_scores))
        if (length(which_max) > 1) {
            max_F1_score_ind <- sample(which(F1_scores %in% max(F1_scores)), size = 1)
        } else {
            max_F1_score_ind <- which(F1_scores %in% max(F1_scores))
        }
        
        output[[count]] <- lambda_trials[[max_F1_score_ind]]
        
        
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
                                             EM_out=output[[count]])
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


#' @export
run_models_fast_test <- function(docs, docName, indexName, labelsName, n_class=2, iters=1,
                            trainSize, labelSize, querySize=100, lambda=0.1, seed=NULL, useLabeled=F,
                            queryType="basic_entropy", fixed_words=NULL, dfm=NULL, verbose=T,
                            forceList=T, bound=0, weight_fpfn=T, log_ratio_conv_type="none",
                            mu=0.0001, tau=0.0001, ...) {
    
    #' @title Run Models Fast
    #' @description Runs 6 model specifications in parallel in order to quickly compare results across models.
    #' 
    #' @param docs Matrix of documents.
    #' @param docName String indicating column name containing document text.
    #' @param indexName String indicating column name containing document index.
    #' @param labelsName String indicating column name containing document labels.
    #' @param n_class Integer indicating number of classes to use.
    #' @param trainSize Integer indicating number of documents to set aside as training data.
    #' @param labelSize Integer indicating total number of documents to actively label.
    #' @param querySize Integer indicating the number of documents to be labeled on each active iteration. 
    #' @param lambda Numeric value indicating lambda value to use for "wEM" and "awEM" models.
    #' @param seed Optional integer value to specify seed.
    #' @param useLabeled Boolean used to indicate whether to calculate in-sample statistics using labeled docs or not. F by default.
    #' @param queryType String indicating which type of uncertainty sampling scheme to use. "basic_entropy" by default.
    #' @param fixed_words List of words which remain fixed at their empirical probabilities across all active iteration. Experimental feature.
    #' @param dfm Optional argument to supply an externally generaged dfm.
    #' @param verbose Boolean used to indicate whether to show progress over iterations.
    #' @param ... Optional arguments to pass to `get_dfm()` if a dfm object is not already supplied.
    #' 
    #' @return List of matrices containing model validation scores across each active iteration.
    #' Outputs accuracy, precision, recall, and F1 scores for in-sample and out-of-sample validation.
    
    
    # ---------------------------- WARNINGS -------------------------------- #
    if (labelSize %% querySize != 0) {
        stop("labelSize must be divisible evenly by querySize.")
    }
    
    if (labelSize == querySize) {
        stop("labelSize cannot equal querySize.")
    }
    
    if (!(is.null(seed)) & iters!=length(seed)) {
        stop("Number of supplied seeds must equal the number of iterations.")
    }
    
    if (forceList & queryType == "log_ratio" & log_ratio_conv_type == "maximand") {
        stop("forceList must be FALSE in order for log_ratio maximand convergence
             to work.")
    }
    # ---------------------------- END WARNINGS ------------------------------ #
    
    # ----------------------------  SETUP  ----------------------------------- #
    models <- c("NB", "EM", "wEM", "aNB", "aEM", "awEM")
    n_model <- length(models)
    label_seq <- seq(querySize, labelSize, querySize)
    
    # create arrays if necessary
    if (iters > 1) {
        acc_in <- prec_in <- rec_in <- F1_arr_in <- array(data=NA, dim=c(n_model, length(label_seq), iters))
        if (nrow(docs) != trainSize) {
            acc_out <- prec_out <- rec_out <- F1_arr_out <- array(data=NA, dim=c(n_model, length(label_seq), iters))
        } 
    }
    # --------------------------  END SETUP  --------------------------------- #
    
    
    # --------------------------  MAIN LOOP  --------------------------------- #
    if (verbose) {
        message("Starting model")
    }
    
    for (k in 1:iters) {
        if (iters > 1) message(paste("Starting iteration", k))
        
        # set seed
        if (is.null(seed) == F) {
            set.seed(seed[k])
        }
        
        # seperate training from testing 
        docs_train <- docs %>%
            dplyr::sample_n(trainSize)
        test_index <- docs %>%
            dplyr::filter(!(!!dplyr::sym(indexName)) %in% docs_train[[paste0(indexName)]]) %>%
            dplyr::pull(!!dplyr::sym(indexName))
        
        # Create containers for results
        # Use FBM for parallelization
        accuracy_in <- bigstatsr::FBM(n_model, length(label_seq), init=NA)
        precision_in <- bigstatsr::FBM(n_model, length(label_seq), init=NA)
        recall_in <- bigstatsr::FBM(n_model, length(label_seq), init=NA)
        F1_in <- bigstatsr::FBM(n_model, length(label_seq), init=NA)
        accuracy_out <- bigstatsr::FBM(n_model, length(label_seq), init=NA)
        precision_out <- bigstatsr::FBM(n_model, length(label_seq), init=NA)
        recall_out <- bigstatsr::FBM(n_model, length(label_seq), init=NA)
        F1_out <- bigstatsr::FBM(n_model, length(label_seq), init=NA)
        
        # Create DFM if necessary
        if (is.null(dfm)) {
            dfm <- get_dfm(docs=docs, docName=docName, indexName=indexName, ...)
        }
        
        docsSample <- docs_train %>% dplyr::slice(1:querySize)
        initIndex <- get_index(docsSample, paste0(indexName))
        
        # Initial documents 
        active_init <- query_Label(docs=docs, toLabel=docsSample[1:querySize, ], n_class=n_class,
                                   labels=labels, docName=docName, indexName=indexName,
                                   labelsName=labelsName, handLabel=F)
        
        # parallelized dopar loop
        cl <- parallel::makeCluster(parallel::detectCores() - 1)
        doParallel::registerDoParallel(cl)
        foreach::foreach (model_type = models, .combine=c,
                          .export = c("get_conf_matrix_weight", "get_precision_binary_weighted",
                                      "get_recall_binary_weighted"),
                          .packages=c("dplyr", "activeEMtext", "parallel", "quanteda")) %dopar% {
                              
                              # Specify weight
                              weight <- dplyr::case_when(
                                  model_type == "aNB" | model_type == "NB" ~ 0,
                                  model_type == "aEM" | model_type == "EM" ~ 1,
                                  model_type == "awEM" | model_type == "wEM" ~ lambda
                              )
                              
                              # specify querytype
                              query <- dplyr::case_when(
                                  model_type %in% c("NB", "EM", "wEM") ~ "random",
                                  model_type %in% c("aNB", "aEM", "awEM") ~ queryType
                              )
                              
                              # Run model
                              if (model_type == "awEM") {
                                  output <- active_EM_varylambda2(docs=active_init, labels=c(0, 1), weight=weight,
                                                                 n_class=2, bound=bound, maxQuery=querySize, 
                                                                 maxActive=labelSize/querySize,
                                                                 docName=docName, indexName=indexName,
                                                                 labelsName=labelsName, handLabel=F,
                                                                 initIndex=initIndex, forceList=forceList, counter_on=F,
                                                                 whichOutTest=test_index,
                                                                 queryType=query, exportAll=T, 
                                                                 fixed_words=fixed_words, dfm=dfm,
                                                                 log_ratio_conv_type=log_ratio_conv_type,
                                                                 mu=mu, tau=tau)
                              } else {
                                  output <- active_EM(docs=active_init, labels=c(0, 1), weight=weight,
                                                      n_class=2, bound=bound, maxQuery=querySize, 
                                                      maxActive=labelSize/querySize,
                                                      docName=docName, indexName=indexName,
                                                      labelsName=labelsName, handLabel=F,
                                                      initIndex=initIndex, forceList=forceList, counter_on=F,
                                                      whichOutTest=test_index,
                                                      queryType=query, exportAll=T, 
                                                      fixed_words=fixed_words, dfm=dfm,
                                                      log_ratio_conv_type=log_ratio_conv_type,
                                                      mu=mu, tau=tau)
                              }
                              
                              
                              
                              # Assign results
                              j <- case_when(
                                  model_type == "NB" ~ 1,
                                  model_type == "EM" ~ 2,
                                  model_type == "wEM" ~ 3,
                                  model_type == "aNB" ~ 4,
                                  model_type == "aEM" ~ 5,
                                  model_type == "awEM" ~ 6
                              )
                              
                              # weighted FP and FN to compute precision and recall
                              if (weight_fpfn == T){
                                  for (i in 1:length(output$docs)) {
                                      
                                      in_conf_matrix <- get_conf_matrix_weight(docs=output$docs[[i]], labelsName=labelsName, indexName=indexName,
                                                                               labeledIndex=output$handLabeledIndex[[i]], levels=c(0,1), n_class=n_class,
                                                                               useLabeled=useLabeled)
                                      
                                      
                                      if (nrow(docs) != trainSize) {
                                          out_conf_matrix <- get_conf_matrix_weight(docs=output$out_docs[[i]], labelsName=labelsName,
                                                                                    levels=c(0, 1), n_class=n_class)
                                      }
                                      
                                      # Get in-sample validation statistics
                                      accuracy_in[j, i] <- get_classification_accuracy(in_conf_matrix$conf_matrix)
                                      precision_in[j, i] <- get_precision_binary_weighted(in_conf_matrix$conf_matrix, in_conf_matrix$FP_w)
                                      recall_in[j, i] <- get_recall_binary_weighted(in_conf_matrix$conf_matrix, in_conf_matrix$FN_w)
                                      F1_in[j, i] <- get_F1_binary(in_conf_matrix$conf_matrix, precision_in[j,i], recall_in[j,i])
                                      
                                      # Get out-of-sample validation statistics
                                      if (nrow(docs) != trainSize) {
                                          accuracy_out[j, i] <- get_classification_accuracy(out_conf_matrix$conf_matrix)
                                          precision_out[j, i] <- get_precision_binary_weighted(out_conf_matrix$conf_matrix, out_conf_matrix$FP_w)
                                          recall_out[j, i] <- get_recall_binary_weighted(out_conf_matrix$conf_matrix, out_conf_matrix$FN_w)
                                          F1_out[j, i] <- get_F1_binary(out_conf_matrix$conf_matrix, precision_out[j,i], recall_out[j,i])
                                      }
                                      
                                  }
                              }else{
                                  for (i in 1:length(output$docs)) {
                                      
                                      in_conf_matrix <- get_conf_matrix(docs=output$docs[[i]], labelsName=labelsName, indexName=indexName,
                                                                        labeledIndex=output$handLabeledIndex[[i]], levels=c(0,1), n_class=n_class,
                                                                        useLabeled=useLabeled)
                                      
                                      
                                      if (nrow(docs) != trainSize) {
                                          out_conf_matrix <- get_conf_matrix(docs=output$out_docs[[i]], labelsName=labelsName,
                                                                             levels=c(0, 1), n_class=n_class)
                                      }
                                      
                                      # Get in-sample validation statistics
                                      accuracy_in[j, i] <- get_classification_accuracy(in_conf_matrix)
                                      precision_in[j, i] <- get_precision_binary(in_conf_matrix)
                                      recall_in[j, i] <- get_recall_binary(in_conf_matrix)
                                      F1_in[j, i] <- get_F1_binary(in_conf_matrix)
                                      
                                      # Get out-of-sample validation statistics
                                      if (nrow(docs) != trainSize) {
                                          accuracy_out[j, i] <- get_classification_accuracy(out_conf_matrix)
                                          precision_out[j, i] <- get_precision_binary(out_conf_matrix)
                                          recall_out[j, i] <- get_recall_binary(out_conf_matrix)
                                          F1_out[j, i] <- get_F1_binary(out_conf_matrix)
                                      }
                                      
                                  }
                              }
                              NULL
                              
                          }
        
        parallel::stopCluster(cl)
        
        # --------------------------  END MAIN LOOP  --------------------------------- #
        
        
        # --------------------------  RESHAPE DATA  ---------------------------------- #
        
        # If running monte carlo, impute results to array
        if (iters > 1) {
            if (nrow(docs) != trainSize) {
                acc_in[, , k] <- accuracy_in[]
                prec_in[, , k] <- precision_in[]
                rec_in[, , k] <- recall_in[]
                F1_arr_in[, , k] <- F1_in[]
                acc_out[, , k] <- accuracy_out[]
                prec_out[, , k] <- precision_out[]
                rec_out[, , k] <- recall_out[]
                F1_arr_out[, , k] <- F1_out[]
            } else {
                acc_in[, , k] <- accuracy_in[]
                prec_in[, , k] <- precision_in[]
                rec_in[, , k] <- recall_in[]
                F1_arr_in[, , k] <- F1_in[]
            }
        }
    }
    
    if (iters > 1) {
        if (nrow(docs) != trainSize) {
            stats <- list(acc_in=acc_in, prec_in=prec_in, rec_in=rec_in, F1_arr_in=F1_arr_in,
                          acc_out=acc_out, prec_out=prec_out, rec_out=rec_out, F1_arr_out=F1_arr_out)
        } else {
            stats <- list(acc_in=acc_in, prec_in=prec_in, rec_in=rec_in, F1_arr_in=F1_arr_in)
        }
        
        for (i in 1:length(stats)) {
            stats[[i]] <- apply(stats[[i]], 1:2, mean, na.rm=TRUE)
        }
    } else {
        if (nrow(docs) != trainSize) {
            stats <- list(accuracy_in=accuracy_in[], precision_in=precision_in[], recall_in=recall_in[],
                          F1_in=F1_in[], accuracy_out=accuracy_out[], precision_out=precision_out[],
                          recall_out=recall_out[], F1_out=F1_out[])
        } else {
            stats <- list(accuracy_in=accuracy_in[], precision_in=precision_in[], recall_in=recall_in[],
                          F1_in=F1_in[])
        }
        
    }
    
    # clean output
    stats <- lapply(stats, function(x) {rownames(x) <- c("NB", "EM", "wEM", "aNB", "aEM", "awEM"); x})
    stats <- lapply(stats, function(x) {colnames(x) <- label_seq; x})
    
    # Remove last column for in-label validation if hand labeled docs are not used for validation
    if (useLabeled == F & trainSize == labelSize) {
        for (i in 1:4) {
            stats[[i]] <- stats[[i]][, -ncol(stats[[i]])]
        }
    }
    
    # --------------------------  END RESHAPE DATA  ------------------------------------------------- #
    
    return(stats)
    
}

