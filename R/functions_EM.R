##################################################
## Project: active
## Script purpose: Functions for EM Algorithm
## Date: 2019/6/29
## Author: Saki Kuzushima
##################################################

EM <- function(.D_train = NULL, .C_train = NULL, .D_test, .n_class = 2,
               .lambda = 0.1, .max_iter = 100, .alpha = 0.1, .lazy_eval = F,
               .counter_on = T, .active_iter = NULL, .maxactive_iter = NULL,
               .fixed_words = NULL, .supervise = T,
               .class_prob = NULL, .word_prob = NULL, .export_all = F){

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
      class_prob <- get_class_prob_NB(.D_train = .D_train, .C_train = .C_train)
      word_prob <- get_word_prob_NB(.D = .D_train, .C_train = .C_train, .beta=beta)

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
      E_log <- E_step(
        .D_test = .D_test, .class_prob = class_prob, .word_prob = word_prob)

        E_log_w <- rbind(log(.C_train), log(.lambda) + E_log)

        # For the output
        if (.supervise) {
            E <- rbind(.C_train, exp(E_log))
        } else {
            E <- exp(E_log)
        }


        # M step
        # DO NOT USE! _LOG FUNCTIONS ARE NOT TESTED YET
        class_prob <- get_class_prob_log(.C_train=.C_train,
                                         .D_train=.D_train,
                                         .D_test=.D_test,
                                         .E=E_log_w,
                                         .n_class=.n_class,
                                         .lambda=.lambda,
                                         .supervise=T)

        word_prob <- get_word_prob_log(.C_train=.C_train,
                                       .D_train=.D_train,
                                       .D_test=.D_test,
                                       .E=E_log_w,
                                       .n_class=.n_class,
                                       .lambda=.lambda)


        # calculate maximand (Question (MB): does the maximand need to be weighted?)
        maximands[count] <- get_maximand(.D_train=.D_train, .D_test=.D_test, .E=E,
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




#' @export
multiEM <- function(.D_train=NULL, .C_train=NULL, .D_test,
                    .D_test_out = NULL,
                    .n_class=2, .n_cluster=2,
                    .lambda=0.1, .max_iter=100,
                    .alpha=0.1, .labeled_docs = NULL,
                    .counter_on=T, .active_iter=NULL,
                    .maxactive_iter=NULL, .fixed_words=NULL,
                    .export_all=F, .supervise=T,
                    .choose_NB_init = FALSE,
                    .prev_word_prob = NULL, .prev_class_prob = NULL,
                    .prev_mu = NA, .prev_psi = NA,
                    .beta=NULL, .binary_metadata_varnames = NA,
                    .cont_metadata_varnames = NA) {

  #' @title EM algorithm
  #'
  #' @description Use EM algorithm to maximize the marginal posterior.
  #' The marginal posterior is the probability of parameters
  #' given both labeled and unlabeled documents and the labels for the labeled documents
  #'
  #' @param .D_train document term matrix of the labeled documents
  #' @param .C_train vector of class labels for the labeled documents
  #' @param .D_test document term matrix of the unlabeled documents
  #' @param .D_test_out document term matrix for out of sample validation
  #' @param .n_class  number of classes
  #' @param .n_cluster number of clusters
  #' @param .lambda  vector of document weights
  #' @param .max_iter maximum number of iteration of the EM algorithm
  #' @param .alpha the threshold of the convergence. If the increase of the maximand becomes less than alpha,
  #' the iteration stops.
  #' @param .labeled_docs Optional vector of index values for labeled documents.
  #' Used if `.choose_NB_init == FALSE`
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
  #' @param .choose_NB_init boolean object. By default true, and EM starts with Naive Bayes step.
  #' If false, and if an appropriate `.C_train` is provided, the initial M step is performed with
  #' document class probabilities from both labeled and unlabeled documents, as weighted by the chosen
  #' `.lambda` value.
  #' @param .cont_metadata_varnames vector of strings indicating variable names of continuous metadata
  #' @param .binary_metadata_varnames vector of stricts indicating variable names of binary metadata
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


  ## Start counter
  if (.counter_on) {
    if (!is.null(.active_iter) & !is.null(.maxactive_iter)) {
      pb <- progress::progress_bar$
        new(format = paste0(
              "Active Iter: ",
              .active_iter, "/",
              .maxactive_iter, " ",
              "EM Runs: [:bar] :current/:total (max)"),
            total = .max_iter)
    } else {
      pb <- progress::progress_bar$
        new(
          format = "EM Runs: [:bar] :current/:total (max)",
          total = .max_iter)
    }
    pb$tick(0)
  }

  ## initialize the container
  maximands <- rep(0, .max_iter)
  if (.export_all) {
    class_probs <- list()
    word_probs <- list()
    class_liks <- list()
    if (!is.na(.binary_metadata_varnames)) {
      psi_vals <- list()
    }
    if (!is.na(.cont_metadata_varnames)) {
      mu_vals <- list()
      sig_vals <- list()
    }
  }

  ## initialize count
  count <- 1

  ## For labeled documents that has negative labels,
  ## cluster probability are uniformly assigned
  ## TODO: Investigate if this creates instability across iterations.
  if (.choose_NB_init) {
    if (.n_class == 2) {
      cluster_train <- matrix(0, nrow = nrow(.C_train), ncol = .n_cluster)
      for (i in 1:nrow(cluster_train)){
        if (.C_train[i, 2] == 1) {
          ## if positive label, we flag the last cluster
          cluster_train[i, .n_cluster] <- 1
        } else {
          ## if negative, we randomly chose a cluster (but not the last cluster)
          pick <- sample(seq(1, .n_cluster - 1), 1)
          cluster_train[i, pick] <- 1
        }
      }
    } else {
      cluster_train <- matrix(0, nrow = nrow(.C_train), ncol = .n_cluster)
      # for (i in 1:nrow(cluster_train)){
      #   for (j in 1:n_cluster) {
      #     if (.C_train[i, j] == 1) {
      #       ## if positive label, we flag the last cluster
      #       cluster_train[i, j] <- 1
      #     }
      #   }
      # }

      cluster_train <- .C_train
    }
  }

  ## Define total document term matrix
  D_all <- rbind(.D_train, .D_test)

  ## extract metadata from DFMs if needed
  if (!is.na(.binary_metadata_varnames)) {
    .X_b_train <- as.matrix(.D_train@docvars[.binary_metadata_varnames])
    .X_b_test <- as.matrix(.D_test@docvars[.binary_metadata_varnames])
    .X_b_all <- rbind(.X_b_train, .X_b_test)
  }
  if (!is.na(.cont_metadata_varnames)) {
    .X_c_train <- as.matrix(.D_train@docvars[.cont_metadata_varnames])
    .X_c_test <- as.matrix(.D_test@docvars[.cont_metadata_varnames])
    .X_c_all <- rbind(.X_c_train, .X_c_test)
  }
  ## Convert beta for multiple cluster case
  if (.n_cluster > 2){
    .binary_beta <- .beta
    .beta <- matrix(0, nrow=nrow(.beta), ncol=.n_cluster)
    ## Copy the value of beta for negative class to 
    ## all clusters linked to the negative class
    .beta[, seq(1, .n_cluster-1)] <- .binary_beta[,1] 
    ## Copy the value of beta for positive class to 
    ## the last cluster linked to the positive class 
    .beta[, .n_cluster] <- .binary_beta[,2] 
  }

  ## Initialize with naive step, or use previous active step's params
  if (.choose_NB_init) {
    class_prob <- get_class_prob_NB(
      .D_train = .D_train,
      .C_train = cluster_train
    )
    word_prob <- get_word_prob_NB(
      .D = .D_train,
      .C_train = cluster_train,
      .beta=.beta
    )
    if (!is.na(.binary_metadata_varnames)) {
      psi <- get_psi(
        .X_b = .X_b_train,
        .E = cluster_train
      )
    }
    if (!is.na(.cont_metadata_varnames)) {
      mu <- get_mu(
        .X_c = .X_c_train,
        .E = cluster_train
      )
      sig <- get_sig(
        .X_c = .X_c_train,
        .E = cluster_train,
        .mu = mu
      )
    }
  } else if (!.choose_NB_init) {
    class_prob <- .prev_class_prob
    word_prob <- .prev_word_prob
    if (!is.na(.prev_mu)) {
      mu <- .prev_mu
      sig <- .prev_sig
    }
    if (!is.na(.prev_psi)) {
      psi <- .prev_psi
    }
  }

  if (.export_all) {
    class_probs[[count]] <- class_prob
    word_probs[[count]] <- word_prob
    if (!is.na(.cont_metadata_varnames)) {
      mu_vals[[count]] <- mu
    }
    if (!is.na(.binary_metadata_varnames)) {
      psi_vals[[count]] <- psi
    }
  }

  ## tentatively set max iteration
  while (count <= .max_iter){

    ## E step
    E_out <- E_step_multi(
      .C_train = .C_train, .D_train = .D_train,
      .X_c = `if`(is.na(.cont_metadata_varnames[1]), NA, .X_c_all),
      .X_b = `if`(is.na(.binary_metadata_varnames[1]), NA, .X_b_all),
      .D_test = .D_test, .class_prob = class_prob,
      .word_prob = word_prob,
      .mu = `if`(is.na(.cont_metadata_varnames[1]), NA, mu),
      .sig = `if`(is.na(.cont_metadata_varnames[1]), NA, sig),
      .psi = `if`(is.na(.binary_metadata_varnames[1]), NA, psi),
      .n_class = .n_class
    )

    ## combine known label and estimated label
    ## prob with lambda weighting for M step.
    ## NOTE: Use E_w for weighting of class and
    ## word probabilities in M step, but report E
    E_w <- rbind(E_out$E[!E_out$unlab, ],
                 E_out$E[E_out$unlab, ] + log(.lambda))
    E <- exp(E_out$E)

    ## M step
    class_prob <- get_class_prob(
      .D_train = .D_train, .D_test = .D_test, .E = exp(E_w),
      .n_class = .n_class, .lambda = .lambda,
      .supervise = .supervise
    )

    word_prob <- get_word_prob(
      .D = D_all, .E = exp(E_w), .beta = .beta,
      .fixed_words = .fixed_words
    )

    if (!is.na(.cont_metadata_varnames[1])) {
      mu <- get_mu(
        .X_c = .X_c_all,
        .E = exp(E_w)
      )
      sig <- get_sig(
        .X_c = .X_c_all,
        .E = exp(E_w),
        .mu = mu
      )
    }
    if (!is.na(.binary_metadata_varnames[1])) {
      psi <- get_psi(
        .X_b = .X_b_all,
        .E = exp(E_w)
      )
    }

    ## calculate maximand
    maximands[count] <- get_maximand(
      .D_train = .D_train, .D_test = .D_test,
      .X_c_train = .X_c_train, .X_c_test = .X_c_test,
      .X_b_train = .X_b_train, .X_b_test = .X_b_test,
      .E = exp(E_w),
      .class_prob = class_prob, .word_prob = word_prob,
      .mu = `if`(is.na(.cont_metadata_varnames)[1], NA, mu),
      .psi = `if`(is.na(.binary_metadata_varnames)[1], NA, psi),
      .sig = `if`(is.na(.cont_metadata_varnames)[1], NA, sig),
      supervise = .supervise
    )

    ## check convergence by maximands
    if (count > 1){
      if (maximands[count] - maximands[count - 1] < .alpha){
        break
      }
    }

    ## update counter
    if (.counter_on == T) {
      pb$tick()
    }

    ## update count
    count <- count + 1

    if (.export_all) {
      class_probs[[count]] <- class_prob
      word_probs[[count]] <- word_prob
      class_liks[[count - 1]] <- E
      if (!is.na(.cont_metadata_varnames)) {
        mu_vals[[count]] <- mu
        sig_vals[[count]] <- sig
      }
      if (!is.na(.binary_metadata_varnames)) {
        psi_vals[[count]] <- psi
      }
    }
  }

  ## calculate out of sample prediction
  if (!is.null(.D_test_out) && nrow(.D_test_out) > 0) {
    out_prediction <- E_step(
      .D_test = .D_test_out,
      .X_c = as.matrix(.D_test_out@docvars[.cont_metadata_varnames]),
      .X_b = as.matrix(.D_test_out@docvars[.binary_metadata_varnames]),
      .class_prob = class_prob,
      .word_prob = word_prob,
      .mu = `if`(is.na(.cont_metadata_varnames), NA, mu),
      .psi = `if`(is.na(.binary_metadata_varnames), NA, psi),
      .sig = `if`(is.na(.cont_metadata_varnames), NA, sig)
    ) %>% exp()
  } else {
    out_prediction <- NULL
  }

  ## calculate weights
  if (.export_all) {
    ## TODO: Change 2 to number of cluster
    ratio <- word_probs[[count]][, 2] - word_probs[[count]][, 1]
  } else {
    ratio <- word_prob[, 2] - word_prob[ , 1]
  }
  lr_weights <- D_all %*% as.matrix(ratio)

  if (.export_all) {
    rownames(word_probs[[1]]) <- rownames(word_probs[[length(word_probs)]])
    return(list(
      maximands = maximands[1:count],
      pi = class_probs,
      eta = word_probs,
      mu = `if`(is.na(.cont_metadata_varnames), NA, mu_vals),
      psi = `if`(is.na(.binary_metadata_varnames), NA, psi_vals),
      sig = `if`(is.na(.cont_metadata_varnames), NA, sig_vals),
      classLik = class_liks,
      out_prediction = out_prediction,
      lr_weights = lr_weights
    ))
  } else {
    return(list(
      maximands = maximands[1:count],
      pi = class_prob,
      eta = word_prob,
      mu = `if`(is.na(.cont_metadata_varnames), NA, mu),
      psi = `if`(is.na(.binary_metadata_varnames), NA, psi),
      sig = `if`(is.na(.cont_metadata_varnames), NA, sig),
      classLik = E,
      out_prediction = out_prediction,
      lr_weights = lr_weights))
  }
}
