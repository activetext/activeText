##################################################
## Project: active
## Script purpose: stable Functions for EM Algorithm (2 cluster)
## Date: 2020/2/7
## Author: Saki Kuzushima
##################################################

# Note thaat this code is commit 48d2526
# https://github.com/ksaki/active/blob/48d2526877122a7fef7c77aa562bd60d66c316db/package/R/functions_EM.R
# This file is a stable version of EM function. 
# Later, I modified M step so that it take the log scaled output of E step,
# and it computes most things in log scale, in order to avoid numerical underflow. (under development)


# E step
E_step <- function(.D_test, .class_prob, .word_prob){
    #' E step of the EM algorithm
    #' @param .D_test document term matrix for the test documents
    #' @param .class_prob log of the class probability 
    #' @param .word_prob log of the word probability 
    #' @return a matrix of the log probability of each doc having each class 
    #' dimension is length(D_test) by 2
    
  
    num <-  sweep(.D_test %*% .word_prob, 2, .class_prob, '+')
    den <- apply(num, 1, matrixStats::logSumExp)
    out <- sweep(num, 1, den, "-")  
    return(out)
}

# M step 
get_word_prob <- function(.D, .E, .fixed_words=NULL){
    #' get log likelihood of a document
    #' @param .D document term matrix
    #' @param .E labels for all documents
    #' @param .fixed_words matrix of fixed words with class probabilities,
    #' where ncol is the number of classes.
    #' @return a matrix of word probability(nrow = 2, ncol = V (size of vocab))
    
    # size of vocabulary
    V <- ncol(.D)
    
    # Eq. (5), took log
    # numerator: 1 + the number of words t in thetraining docs with class j 
    # denominator: V + the number of words in the training docs with class j 
    
    num <- log(1 + Matrix::t(.D) %*% .E) 
    den <- log(V + Matrix::colSums(Matrix::t(.D) %*% .E))
    out <- sweep(num, 2, den, "-")
    
    if (!is.null(.fixed_words)) {
        out2 <- rbind(out[!(rownames(out) %in% rownames(.fixed_words)), ], .fixed_words) 
        out <- out2[match(rownames(out), rownames(out2)), ]
    }
    
    return(out)
}

# Weighted M step 
get_class_prob <- function(.D_train, .D_test, .E, .n_class, .lambda, .supervise){
    
    #' get weighted log of class probability
    #' @param .D_train document term matrix for the training documents
    #' @param .D_test document term matrix for the test documents
    #' @param .E output of E step for all documents
    #' @param .n_class  number of classes
    #' @param .lambda weight
    #' @param supervise T if supervise, F is unsupervise
    #' @return a matrix of log likelihood of a class(nrow = 2, ncol = # of document)
    
    # number of the labeled and unlabeled documents, combined
    test_len <- ifelse(is.null(nrow(.D_test)), 0, nrow(.D_test))
    if (.supervise == T) {
      train_len <- nrow(.D_train)
    } else {
      train_len <- 0
    }
    
    out <- log(1 + Matrix::colSums(.E)) - log(.n_class + train_len + .lambda * test_len)
    return(out)
}

get_maximand <- function(.D_train, .D_test, .E, .class_prob, .word_prob, supervise){
    #' get maximand. Use this to check if the EM keeps incrasing maximand
    #' @param .D_train document term matrix for the training documents
    #' @param .D_test document term matrix for the test documents
    #' @param .E output of E step for all documents
    #' @param .class_prob  output of get_class_prob()
    #' @param .word_prob output of get_word_prob()
    #' @param .supervise T if supervised. F is unsupervised
    #' @return a scalar of the value of maximand in each EM iteration
     
    
    #llk
    if (supervise == T) {
      # train doc
      train_doc_prob <- sum(.E[1:nrow(.D_train), ] * sweep(.D_train %*% .word_prob, 2, .class_prob, '+'))
    } else {
      train_doc_prob <- 0
    }
   
  
    # test doc
    insidelog <- sweep(.D_test %*% .word_prob, 2, .class_prob, '+')
    test_doc_prob <- sum(apply(insidelog, 1, matrixStats::logSumExp))
    
    doc_llk <- train_doc_prob + test_doc_prob
    
    #prior
    prior <- sum(.class_prob + Matrix::colSums(.word_prob))
    
    
    return(doc_llk + prior) 
}

get_model_diff <- function(.class_prob, .class_prob_prev, .word_prob, .word_prob_prev) {
    #' @title Get Difference Between Model Parameters
    #'
    #' @description Calculates total difference in parameter log likelihood between two runs of EM model.
    #' 
    #' @param class_prob Vector of current class probabilities
    #' @param class_prob_prev Vector of class probabilities from previous iteration
    #' @param word_prob Matrix of current word probabilities
    #' @param word_prob_prev Matrix of previous word probabilities
    #' 
    #' @return Total sum difference of model parameters.
    
    class_diff <- sum(abs(.class_prob - .class_prob_prev))
    word_diff <- sum(abs(.word_prob - .word_prob_prev))
    return(class_diff + word_diff)
}


# main ==============================================================
#' @export
EM_2cl <- function(.D_train=NULL, .C_train=NULL, .D_test, .n_class=2, .lambda=0.1, .max_iter=100, .alpha=0.1, .lazy_eval=F,
               .counter_on=T, .active_iter=NULL, .maxactive_iter=NULL, .fixed_words=NULL, .supervise = T, 
               .class_prob=NULL, .word_prob=NULL, .export_all=F){

    #' @title EM algorithm
    #' 
    #' @description Use EM algorithm to maximize the marginal posterior. 
    #' The marginal posterior is the probability of parameters 
    #' given both labeled and unlabeled documents and the labels for the labeled documents
    #' 
    #' @param .D_train document term matrix of the labeled documents
    #' @param .C_train vector of class labels for the labeled documents
    #' @param .D_test document term matrix of the unlabeled documents
    #' @param .n_class  number of classes
    #' @param .lambda  vector of document weights
    #' @param .max_iter maximum number of iteration of the EM algorithm
    #' @param .alpha the threshold of the convergence. If the increase of the maximand becomes less than alpha, 
    #' the iteration stops. 
    #' @param .lazy_eval boolean object. If \code{lazy_eval == T}, convergence is measured by comparing changes
    #' in log likelihood across model iterations rather than directly computing maximand.
    #' @param .counter_on boolean object. If \code{counter_on == T}, displays the progress of the EM
    #' algorithm.
    #' @param .active_iter integer value that tells the EM algorithm which iteration of the active
    #' loop it is in.
    #' @param .maxactive_iter integer value that tells the EM algorithm the maximum allowed active
    #' iterations.
    #' @param .fixed_words matrix of fixed words with class probabilities,
    #' where ncol is the number of classes.
    #' @param .supervise T if supervised. F is unsupervised.
    #' @param .class_prob required if .supervise == T. Starting value of class probability (logged)
    #' @param .word_prob required if .supervise == T. Starting value of word probability (logged)
    #' @param .export_all If T, model parameters from each iteration of the EM algorithm are returned.
    #' If F, only model results from the last iteration are returned.
    #' 
    #' @return maximands is a vector of maximands in each iteration. Each element of the vector 
    #' contains the log maximand in each step.
    #' pi is a vector of log class probabilities. (length = 2)
    #' eta is a matrix of log word probabilities (nrow = the number of all documents, ncol = 2)
    #' 
    #' @details The inputs must conform to the following specifications
    #' D_train: a matrix with dimension: the number of labeled documents * the number of unique words
    #' D_test: a matrix with dimension: the number of labeled documents * the number of unique words
    #' The column length of D_train and D_test must be the same. The elements of the D_train, D_test are
    #' integers (the counts of each unique word appeard in each document)
    #' C_train: vector of labels for the labeled documents. The length must be the same as the row length of D_test
    #' 
    #' @references Active EM in overleaf 
    #' 
    #' 

    # Start counter
    if (.counter_on) {
        if (!is.null(.active_iter) & !is.null(.maxactive_iter)) {
            pb <- progress::progress_bar$new(
                format = paste0("Active Iter: ", .active_iter, "/",
                                .maxactive_iter, " ",
                                "EM Runs: [:bar] :current/:total (max)") , 
                total = .max_iter)
        } else {
            pb <- progress::progress_bar$new(
                format = "EM Runs: [:bar] :current/:total (max)", 
                total = .max_iter)
        }
        pb$tick(0)
        Sys.sleep(0.02)
    }
    
    # initialize the container
    maximands <- rep(0, .max_iter)
    if (.export_all) {
        class_probs <- list()
        word_probs <- list()
        class_liks <- list()
    }

    
    # initialize count
    count <- 1
  
    
    if (.supervise == T) {
      # Naive step
      class_prob <- get_class_prob(.D_train=.D_train, .D_test=0, .E=.C_train, 
                                   .n_class=.n_class, .lambda=.lambda, .supervise = T)
      word_prob <- get_word_prob(.D=.D_train, .E=.C_train, .fixed_words=.fixed_words)
      
      # Define total document term matrix
      D_all <- rbind(.D_train, .D_test)
    } else if (is.null(.class_prob) | is.null(.word_prob)) {
      # use the following initial values if class-prob or word_prob were not specified
      
      # set prior by uniform
      tmp <- 0.5
      class_prob <- log(c(tmp, 1 - tmp))
      # set likelihood by Dir(2...2)
      word_prob <- log(t(MCMCpack::rdirichlet(n = 2, alpha = rep(2, ncol(.D_test)))))
      D_all <- .D_test
    } else {
      # use custom initial values as a starting value of EM
      class_prob <- .class_prob
      word_prob <- .word_prob
      
      # no training documents
      D_all <- .D_test
    }
    
    if (.export_all) {
        class_probs[[count]] <- class_prob
        word_probs[[count]] <- word_prob
    }
   
    
    # tentatively set max iteration
    while (count <= .max_iter){
        # E step
        E_out <- E_step(.D_test=.D_test, .class_prob=class_prob, .word_prob=word_prob)
        
        # combine known label and estimated label prob with lambda weighting for M step
        # Use E_w for weighting of class and word probabilities in M step, but report E
        
        if (.supervise) {
            E <- rbind(.C_train, exp(E_out))
            E_w <- rbind(.C_train, .lambda * exp(E_out))
        } else {
            E <- E_w<- exp(E_out)
        }

        
        # M step
        class_prob <- get_class_prob(.D_train=.D_train, .D_test=.D_test, .E=E_w, 
                                     .n_class=.n_class, .lambda=.lambda, .supervise=.supervise) 
        word_prob <- get_word_prob(.D=D_all, .E=E_w, .fixed_words=.fixed_words)
        

        # calculate maximand (Question (MB): does the maximand need to be weighted?)
        maximands[count] <- get_maximand(.D_train=.D_train, .D_test=.D_test, .E=E_w, 
                                         .class_prob=class_prob, .word_prob=word_prob, 
                                         supervise=.supervise)
        
        if (.lazy_eval == T) {
            # Lazy convergence evaluation
            if (count > 1) {
                model_diff <- get_model_diff(.class_prob=class_prob, .class_prob_prev=class_prob_prev,
                                             .word_prob=word_prob, .word_prob_prev=word_prob_prev)
                if (model_diff < .alpha) {
                    break
                }
            }
            # save current model parameters for convergence evaluation
            class_prob_prev <- class_prob
            word_prob_prev <- word_prob
        } else {
          # check convergence by maximands      
          if (count > 1){
            if (maximands[count] - maximands[count - 1] < .alpha){
              break
            }
          }
        }
        
        if (.counter_on == T) {
            # update counter
            pb$tick()
            Sys.sleep(.075)
        }

        # update count
        count <- count + 1
        
        if (.export_all) {
            class_probs[[count]] <- class_prob
            word_probs[[count]] <- word_prob
            class_liks[[count - 1]] <- E
        }
    
    }
    
    # calculate weights
    if(.export_all) {
        ratio <- word_probs[[count]][, 2] - word_probs[[count]][, 1]
    } else {
        ratio <- word_prob[, 2] - word_prob[ , 1]
    }
    lr_weights <- D_all %*% as.matrix(ratio)
   
    
    if (.export_all) {
        rownames(word_probs[[1]]) <- rownames(word_probs[[length(word_probs)]])
        return(list(maximands=maximands[1:count],
                    pi=class_probs,
                    eta=word_probs,
                    classLik=class_liks,
                    lr_weights=lr_weights))
    } else {
        return(list(maximands=maximands[1:count],
                    pi=class_prob,
                    eta=word_prob,
                    classLik=E,
                    lr_weights=lr_weights))
    }
    

}