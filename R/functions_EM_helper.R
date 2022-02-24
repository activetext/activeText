##################################################
## Project: active
## Script purpose: Helper functions for EM
## Date: 2020/2/1
## Author: Saki Kuzushima
##################################################


# E step
E_step <- function(.D_test, .X_c = NA, .X_b = NA,
                   .class_prob, .word_prob,
                   .mu = NA, .sig = NA, .psi = NA){
    #' E step of the EM algorithm
    #' @param .D_test document term matrix for the test documents
    #' @param .class_prob log of the class probability 
    #' @param .word_prob log of the word probability 
    #' @return a matrix of the log probability of each doc having each class 
    #' dimension is length(D_test) by 2

  lk <- .D_test %*% .word_prob
  if (!is.null(nrow(.mu))) {
    cont_llk <- matrix(0, nrow = nrow(lk), ncol = ncol(lk))
    for (j in 1:nrow(.mu)) {
      for (k in 1:ncol(.mu)) {
        cont_llk[, k] <- cont_llk[, k] +
          dnorm(
            .X_c[, j, drop = FALSE], .mu[j, k],
            sqrt(.sig[j, k]), log = TRUE
          )
      }
      ## cont_llk <- cont_llk + log(apply(
      ##   .mu[j, , drop = FALSE], 2, function(mu) dnorm(.X_c[ , j, drop = FALSE], mu)
      ## ))
    }
    lk <- lk + cont_llk
  }
  if (!is.null(nrow(.psi))) {
    binary_llk <- 0
    for (j in 1:nrow(.psi)) {
      binary_llk <- binary_llk + apply(
        .psi[j, , drop = FALSE], 2, function(psi) dbinom(.X_b[ , j, drop = FALSE], 1, psi, log = TRUE)
      )
    }
    lk <- lk + binary_llk
  }
  
  num <- sweep(lk, 2, .class_prob, '+')
  den <- apply(num, 1, matrixStats::logSumExp)
  out <- sweep(num, 1, den, "-")
  return(out)
}

#' @export
E_step_multi <- function(.C_train, .D_train, .D_test,
                         .X_c = NA, .X_b = NA,
                         .class_prob, .word_prob,
                         .mu = NA, .psi = NA, .sig = NA){
  
  #' E step of the EM algorithm with multiple clusters
  #' it has to take .D_train (labeled documents) too because
  #' for negative label documents, we have to estimate cluster probability
  #' TO BE UPDATED
  #' @title E step with multiple cluster
  #' @param .C_train: class matrix of the training data
  #' @param .D_train: document term matrix of the training data
  #' @param .D_test document term matrix for the test data
  #' @param .class_prob log of the class probability 
  #' @param .word_prob log of the word probability
  #' @param .mu log of continuous metadata probability
  #' @param .psi log of binary metadata probability
  #' @return a matrix of the log probability of each doc having each class 
  #' dimension is length(D_test) by 2

  # Notation
  # N: number of documents
  # V: size of vocabulary
  # K: number of latent cluster
  # C: number of classes

  # binary vector of positive and negative docs, and unlabelled docs
  poslab <- c(as.logical(.C_train[,2]), rep(F, nrow(.D_test))) # T if positive F ow
  neglab <- c(as.logical(.C_train[,1]), rep(F, nrow(.D_test))) # T if negative F ow
  unlab <- c(rep(F, length.out=nrow(.D_train)), rep(T, length.out=nrow(.D_test))) # T if labeled F ow 
  .D_all <- rbind(.D_train, .D_test)

  lk <- .D_all %*% .word_prob

  mtx_logsumexp <- function(llk_mtx_one, llk_mtx_two) {
    #' Vectorize matrices, stack them, apply logsumexp,
    #' then reshape into matrix form.

    m <- nrow(lk)
    n <- ncol(lk)

    llk_vec_one <- as.vector(llk_mtx_one)
    llk_vec_two <- as.vector(llk_mtx_two)

    llk_vec_sum <- apply(
      rbind(llk_vec_one, llk_vec_two), 2, matrixStats::logSumExp
    )

    llk_mtx_sum <- matrix(llk_vec_sum, nrow = m, ncol = n, byrow = TRUE)
    return(llk_mtx_sum)
  }

  if (!is.null(nrow(.mu))) {
    cont_llk <- matrix(0, nrow = nrow(lk), ncol = ncol(lk))
    for (j in 1:nrow(.mu)) {
      for (k in 1:ncol(.mu)) {
        cont_llk[, k] <- cont_llk[, k] +
          dnorm(
            .X_c[, j, drop = FALSE], .mu[j, k],
            sqrt(.sig[j, k]), log = TRUE
          )
      }
      ## cont_llk <- cont_llk + log(apply(
      ##   .mu[j, , drop = FALSE], 2, function(mu) dnorm(.X_c[ , j, drop = FALSE], mu)
      ## ))
    }
    lk <- lk + cont_llk
  }

  if (!is.null(nrow(.psi))) {
    for (j in 1:nrow(.psi)) {
      binary_llk <- apply(
        .psi[j, , drop = FALSE], 2,
        function(psi) dbinom(.X_b[ , j , drop = FALSE], 1, psi, log = TRUE)
      )
      lk <- lk + binary_llk
      ## lk <- mtx_logsumexp(lk, binary_llk)
    }
  }

  num <- sweep(lk, 2, .class_prob, '+')

  # This is the numerator
  # .D_all = N * V, .word_prob = V * K, .class_prob = 1 * K
  # So, .D_all %*% .word_prob is the document log likelihood. 
  # We will multyply prior, .class_prob, to each document (addition in log scale) 
  # Note that .word_prob and .class_prob are in log scale. 

  den <- apply(num, 1, matrixStats::logSumExp)
  # Take a sum across latent clsuters for each document 
  # This is a denominator

  out <- sweep(num, 1, den, "-")  
  #normalization
  
  # for positive label documents, set probability of positive class to be 1 and others 0
  # Because they are in log scale, I plug a very small value instead of negative infinity

  out[poslab, ncol(out)] <- 0   # last cluster is associated with positive class
  out[poslab, -ncol(out)] <- -Inf   # the other clusters are with negative class
  
  # for negative label documents, 
  # set the cluster associated with positive class to be zero
  # to make exp(x) = 0, I set an extremely small value
  # then normalize

  neglab_out <- num[neglab, -ncol(num), drop=F] # create a temporary object


  den <- apply(neglab_out, 1, matrixStats::logSumExp)
  neglab_out <- sweep(neglab_out, 1, den, "-")
  # normalize
  
  out[neglab, -ncol(out)] <- neglab_out # bring it back
  out[neglab, ncol(out)] <- -Inf
  # for unlabeled documents, usual E step is enough but add weight
  
  return(list(E=out, unlab=unlab))
}




# M step 
get_word_prob <- function(.D, .E, .beta=NA, .fixed_words=NULL){


   # THIS FUNCTION ASSMES THAT .E IS ALREADY WEIGHTED, AND IN LINEAR SCALE

    #' get log likelihood of a document
    #' @param .D document term matrix
    #' @param .E labels for all documents
    #' @param .beta prior parameter for eta
    #' @param .fixed_words matrix of fixed words with class probabilities,
    #' where ncol is the number of classes.
    #' @return a matrix of word probability(nrow = 2, ncol = V (size of vocab))
    
    # size of vocabulary
    V <- ncol(.D)
    
    # Eq. (5), took log
    # numerator: 1 + the number of words t in thetraining docs with class j 
    # denominator: V + the number of words in the training docs with class j

    # if .beta is not provided, use (2,2,...) as beta
    if (is.na(.beta)){
      .beta <- rep(2, ncol(.D))

      # If we want to use multi cluster, fix this vector and make it to a matrix
      .beta_neg <- rep(2, ncol(.D)) # for the negative class, add almost uninformative prior
      .beta <- cbind(.beta_neg, .beta)
      colnames(.beta) <- NULL # make sure .beta does not show up in the col names
    }

    num <- log(.beta - 1 + Matrix::t(.D) %*% .E)
    den <- log(V + Matrix::colSums(Matrix::t(.D) %*% .E))
    out <- sweep(num, 2, den, "-")
    
    if (!is.null(.fixed_words)) {
        out2 <- rbind(out[!(rownames(out) %in% rownames(.fixed_words)), ], .fixed_words) 
        out <- out2[match(rownames(out), rownames(out2)), ]
    }
    
    return(out)
}

#' @export
get_word_prob_NB <- function(.D_train, .C_train, .beta = NA){
  #' @title Get Word Probability (Naive Bayes)
	#' @description get word probability for NB step

	#' @param .D_train DTM for training docs
	#' @param .C_train class matrix for training docs
  #' @param .beta prior parameter for eta (positive column)
  # if .beta is not provided, use (2,2,...) as beta


  if (is.na(.beta)){
    .beta <- rep(2, ncol(.D_train))
    # NOTE:
    # If we want to use multi cluster, fix this vector and make it to a matrix
    .beta_neg <- rep(2, ncol(.D_train)) # for the negative class, add almost uninformative prior
    .beta <- cbind(.beta_neg, .beta)
    colnames(.beta) <- NULL # make sure .beta does not show up in the col names
  }

	num <- log(.beta - 1 + Matrix::crossprod(.D_train, .C_train))
	den <- apply(num, 2, matrixStats::logSumExp)
	return(sweep(num, 2, den, "-"))
}


get_mu <- function(.X_c, .E) {
  num <- Matrix::t(.X_c) %*% .E
  den <- 1 + apply(.E, 2, sum)
  ## den <- apply(.E, 2, sum)
  return(sweep(num, 2, den, "/"))
}

get_sig <- function(.X_c, .E, .mu) {

  m <- ncol(.X_c)
  n <- ncol(.mu)

  calc_sig <- function(.X_c_col, .mu_col, .E_col ) {
    sig <- (Matrix::t(.X_c_col) - .mu_col)^2 %*% .E_col + 1 #+ .mu_col^2 + 1
    return(sig)
  }

  num <- matrix(NA, nrow = m, ncol = n)
  for (j in 1:m) {
    for (i in 1:n) {
      num[j, i] <- calc_sig(.X_c[, m], .mu[j, i], .E[, i])
    }
  }
  den <- 1 + apply(.E, 2, sum)

  return(sweep(num, 2, den, "/"))

}

get_psi <- function(.X_b, .E) {
  num <- 1 + Matrix::t(.X_b) %*% .E
  den <- 2 + apply(.E, 2, sum)
  return(sweep(num, 2, den, "/"))
}


split_sparsemat <- function(.C_train, .D_train){
  #' Split DTM into positive and negative label areas
  #' Extract sparse matrix triplet (row, col, value)


  # sparse representation
  # Split labeled DTM into positive and negative labeled ones
  pos <- as.logical(.C_train[,ncol(.C_train)])
  neg <- !pos

  .D_train_pos <- .D_train[pos,]
  .D_train_neg <- .D_train[neg,]

  # positive label docs
  dtm_train_pos_i <- .D_train_pos@i + 1 # change zero base to 1 base 
  dtm_train_pos_t <- rep(seq_along(diff(.D_train_pos@p)),diff(.D_train_pos@p))
  dtm_train_pos_x <- log(.D_train_pos@x)

  # negative label docs
  dtm_train_neg_i <- .D_train_neg@i + 1 # change zero base to 1 base 
  dtm_train_neg_t <- rep(seq_along(diff(.D_train_neg@p)),diff(.D_train_neg@p))
  dtm_train_neg_x <- log(.D_train_neg@x)

  # unlabeled docs
  dtm_test_i <- .D_test@i + 1 # change zero base to 1 base 
  dtm_test_t <- rep(seq_along(diff(.D_test@p)),diff(.D_test@p))
  dtm_test_x <- log(.D_test@x)

  # create factor variable for the splits later
  dtm_t_lev <- 1:ncol(.D_train) # vocabulary size

  # split x, i vector such that 
  # each element of x_ls contains a vector of the number of words
  # (that are not zero) in log scale
  # each element of i_ls contains a vector of the row numbers
  # that haas at least one word t in the document
  dtm_train_pos_x_ls <- split(dtm_train_pos_x, factor(dtm_train_pos_t, levels=dtm_t_lev))
  dtm_train_pos_i_ls <- split(dtm_train_pos_i, factor(dtm_train_pos_t, levels=dtm_t_lev))

  dtm_train_neg_x_ls <- split(dtm_train_neg_x, factor(dtm_train_neg_t, levels=dtm_t_lev))
  dtm_train_neg_i_ls <- split(dtm_train_neg_i, factor(dtm_train_neg_t, levels=dtm_t_lev))

  dtm_test_x_ls <- split(dtm_test_x, factor(dtm_test_t, levels=dtm_t_lev))
  dtm_test_i_ls <- split(dtm_test_i, factor(dtm_test_t, levels=dtm_t_lev))

  return(list(dtm_train_pos_x_ls=dtm_train_pos_x_ls,
              dtm_train_pos_i_ls=dtm_train_pos_i_ls,
              dtm_train_neg_x_ls=dtm_train_neg_x_ls,
              dtm_train_neg_i_ls=dtm_train_neg_i_ls,
              dtm_test_x_ls=dtm_test_x_ls,
              dtm_test_i_ls=dtm_test_i_ls))

}

get_word_prob_log <- function(.C_train, .D_train, .D_test, .E, .n_class, .lambda){

  #' get log likelihood of a document 
  #' @param .D document term matrix
  #' @param .E output from the E step. Log scale. 
  #' @param .fixed_words matrix of fixed words with class probabilities,
  #' where ncol is the number of classes.
  #' @return a matrix of word probability(nrow = 2, ncol = V (size of vocab))

  
  # storage
  num <- matrix(0.0, nrow=ncol(.D_train), ncol=ncol(.E))


  # indicator for positive label, negative label, unlabel
  pos <- as.logical(.C_train[,2])
  neg <- !pos
  D_poslab <- .D_train[pos,]
  D_neglab <- .D_train[neg,]
  E_neglab <- .E[c(neg, rep(F, nrow(.D_test))),]
  E_unlab <- .E[c(rep(F, nrow(.D_train)), rep(T, nrow(.D_test))),]


  # note; for each? (later)
  for (v in 1:ncol(.D_train)){
    for (k in 1:ncol(.E)){

      # the last cluster is the total the number of word t occuring in the positvely labeled docs.
      if (k == ncol(.E)){
        #labeled_pos <- matrixStats::logSumExp(log(D_poslab[,v]))
        labeled_pos <- log(sum(D_poslab[,v]))
      # The rest of the clusters are zero in linear scale, so set values close to -Inf in log scale.
      }else{
        labeled_pos <- -Inf
      }

      # labeled (negative)
      labeled_neg <- matrixStats::logSumExp(E_neglab[,k] + log(D_neglab[,v]))

      # unlabeled
      unlabeled <- matrixStats::logSumExp(E_unlab[,k] + log(.D_test[,v]))
    
      
      num[v,k] <- matrixStats::logSumExp(c(log(1), labeled_pos, labeled_neg, log(.lambda) + unlabeled))
    }
  }

  # normalize
  den <- apply(num, 2, matrixStats::logSumExp)
  out <- sweep(num, 2, den, "-")

  return(out)

}


get_class_prob <- function(.D_train, .D_test, .E, .n_class, .lambda, .supervise){

    # THIS FUNCTION ASSMES THAT .E IS ALREADY WEIGHTED, AND IN LINEAR SCALE

    
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


get_class_prob_NB <- function(.D_train, .C_train){

	#' NB step, class probability
	#' @param .D_train DTM for training docs
	#' @param .C_train class matrix for training docs
	return(log(1 + Matrix::colSums(.C_train)) - log(2 + nrow(.D_train)))

}

get_class_prob_log <- function(.C_train, .D_train, .D_test, .E, .n_class, .lambda, .supervise){
    
    #' get weighted log of class probability
    #' @param .D_train document term matrix for the training documents
    #' @param .D_test document term matrix for the test documents
    #' @param .E output of E step for all documents (log scale)
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

    temp_lab <- Matrix::colSums(.C_train) # labeled 
    temp_unlab <- apply(.E, 2, matrixStats::logSumExp) # unlabeled
    num <- log(1 + temp_lab + .lambda * exp(temp_unlab))
    out <- num - matrixStats::logSumExp(num)
    
    return(out)
}



get_maximand <- function(.D_train, .D_test, .E, .X_c_train, .X_c_test, .X_b_train, .X_b_test,
                         .class_prob, .word_prob, .mu = NA, .psi = NA, .sig = NA, supervise){
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
    tr_lik <- .D_train %*% .word_prob

    if (!is.null(nrow(.mu))) {
      cont_llk <- matrix(0, nrow = nrow(tr_lik), ncol = ncol(tr_lik))
      for (j in 1:nrow(.mu)) {
        for (k in 1:ncol(.mu)) {
          cont_llk[, k] <- cont_llk[, k] +
            dnorm(
              .X_c_train[, j, drop = FALSE], .mu[j, k],
              sqrt(.sig[j, k]), log = TRUE
            )
        }
      }
      tr_lik <- tr_lik + cont_llk
    }

    if (!is.null(nrow(.psi))) {
      for (j in 1:nrow(.psi)) {
        binary_llk <- apply(
          .psi[j, , drop = FALSE], 2,
          function(psi) dbinom(.X_b_train[ , j , drop = FALSE], 1, psi, log = TRUE)
        )
        tr_lik <- tr_lik + binary_llk
      }
    }

    train_doc_prob <- sum(
      .E[1:nrow(.D_train), ] *
      sweep(tr_lik, 2, .class_prob, '+')
    )
  } else {
    train_doc_prob <- 0
  }
   
  
    # test doc

  tst_lik <- .D_test %*% .word_prob

  if (!is.null(nrow(.mu))) {
    cont_llk <- matrix(0, nrow = nrow(tst_lik), ncol = ncol(tst_lik))
    for (j in 1:nrow(.mu)) {
      for (k in 1:ncol(.mu)) {
        cont_llk[, k] <- cont_llk[, k] +
          dnorm(
            .X_c_test[, j, drop = FALSE], .mu[j, k],
            sqrt(.sig[j, k]), log = TRUE
          )
      }
    }
    tst_lik <- tst_lik + cont_llk
  }

  if (!is.null(nrow(.psi))) {
    for (j in 1:nrow(.psi)) {
      binary_llk <- apply(
        .psi[j, , drop = FALSE], 2,
        function(psi) dbinom(.X_b_test[ , j , drop = FALSE], 1, psi, log = TRUE)
      )
      tst_lik <- tst_lik + binary_llk
    }
  }

  insidelog <- sweep(tst_lik, 2, .class_prob, '+')
  test_doc_prob <- sum(apply(insidelog, 1, matrixStats::logSumExp))
    
  doc_llk <- train_doc_prob + test_doc_prob

  #prior

  prior <- .class_prob + Matrix::colSums(.word_prob)

  prior_sum <- sum(prior)

  mu_prior_lk <- 0
  if (!is.null(nrow(.mu))) {
    for (k in 1:ncol(.mu)) {
      mu_prior_lk <- mu_prior_lk + sum(dnorm(.mu[, k], log = TRUE))
    }
  }
  prior_sum <- prior_sum + mu_prior_lk
  if (!is.null(nrow(.psi))) {
    prior_sum <- prior_sum + sum(Matrix::colSums(.psi))
  }


  return(doc_llk + prior_sum)
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

