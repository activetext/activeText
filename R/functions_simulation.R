#' @export

activeEM_sim <- function (dtm, class_mtx, labels, labeled_index, n_class, weight, whichOutTest=NULL,
                          maxEM=100, EMthreshold=.1, maxActive=10, lazy_eval=F, counter_on=T) {
    
    # Set count at 0 for initialization
    count <- 1
    
    # Split dtm into in-sample and out-sample sections
    dtms <- split_dtm(dtm, whichOutTest)
    
    # Main loop
    output <- list()
    repeat {
        
        # Get EM_dtms
        EM_dtms <- split_dtm(dtm=dtms$second_section, splitIndex=labeled_index)
        
        # Run EM Algorithm
        output[[count]] <- EM(.D_train=EM_dtms$first_section, .C_train=class_mtx, .D_test=EM_dtms$second_section,
                              .n_class=n_class, .lambda=weight, .max_iter=maxEM, .alpha=EMthreshold, .lazy_eval=lazy_eval, 
                              .counter_on=counter_on)
        
        # Update class matrix
        class_mtx <- output[[count]]$classLik
        
        # Check for maximum iterations reached
        if (count == maxActive) break
        
        # Get uncertain documents
        class_ent <- cbind(class_mtx, get_entropy(class_mtx))
        uncertain_docs <- rownames(class_ent[order(class_ent[, ncol(class_mtx) + 1],
                                                   decreasing=T),][1:maxQuery, ])

        # Check for convergence
        if (length(uncertain_docs) == 0) break
        
        # Query oracle for uncertain document labels
        to_label <- which(rownames(class_mtx) %in% uncertain_docs)
        labels_index <- which(rownames(labels) %in% uncertain_docs)
        for (i in 1:length(to_label)) {
            class_mtx[to_label[i], ] <- 0
            if (labels[labels_index[i]] == 0) {
                class_mtx[to_label[i], 1] <- 1
            } else {
                class_mtx[to_label[i], 2] <- 1
            }
            
        }

        # Update labeled index
        labeled_index <- c(labeled_index, uncertain_docs)
        
        # Update counter
        count <- count + 1
    }
    
    # get out of sample prediction
    if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {
        out_prediction <- E_step(.D_test=dtms$first_section, .class_prob=output[[length(output)]]$pi,
                                 .word_prob=output[[length(output)]]$eta)
    }
    
    # return results
    return(
        if (is.null(whichOutTest) == F & length(whichOutTest) != 0) {
            list(EM_output=output, out_prediction=out_prediction)
        } else {
            list(EMoutput=output)
        }
    )

}
