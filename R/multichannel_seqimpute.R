#' @importFrom stringr str_count
#' @importFrom stringr str_detect
#' @importFrom stringr str_locate
#' @importFrom stringr str_locate_all
#' 
#' @importFrom graphics plot
#' 
#' @importFrom stats as.formula
#' @importFrom stats cutree
#' @importFrom stats lm
#' @importFrom stats predict
#' @importFrom stats rnorm
#' @importFrom stats runif
#' 
#' @importFrom utils capture.output
#' 
#' @importFrom Amelia missmap
#' 
#' @importFrom TraMineR seqdef
#' @importFrom TraMineR seqfplot
#' @importFrom TraMineR seqdplot
#' @importFrom TraMineR seqsubm
#' @importFrom TraMineR seqdist
#' 
#' @importFrom cluster agnes
#' 
#' @importFrom plyr mapvalues
#' 
#' @importFrom dfidx dfidx
#' 
#' @importFrom rms lrm
#' 
#' @importFrom mice as.mids
#' 
#' @importFrom mlr makeClassifTask
#' 
#' @importFrom ranger ranger
#' 
#' @importFrom stats model.matrix
#' 
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' 
#' @importFrom dplyr n_distinct
#' 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom seqimpute seqimpute
#' @export
multichannel_seqimpute <- function(list_datasets, np=1, nf=1, m=5, niter=1, 
                                   timing=FALSE, frame.radius=0, covariates=NULL, time.covariates=NULL, 
                                   regr="multinom", npt=1, nfi=1, ParExec=FALSE, ncores=NULL, 
                                   SetRNGSeed=FALSE, verbose=TRUE, available=TRUE, pastDistrib=FALSE,
                                   futureDistrib=FALSE,...){
  
  
  
  #Setting parallel or sequential backend and  random seed
  if (ParExec & (parallel::detectCores() > 2 & m>1)){
    if(is.null(ncores)){
      Ncpus <- parallel::detectCores() - 1
    }else{
      Ncpus <- min(ncores,parallel::detectCores() - 1)
    }
    cl <- parallel::makeCluster(Ncpus)
    doSNOW::registerDoSNOW(cl) #registerDoParallel doesn't have compatibility with ProgressBar
    if(SetRNGSeed){
      doRNG::registerDoRNG(SetRNGSeed)
    }
    # set progress bar for parallel processing
    pb <- txtProgressBar(max = m, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # condition used to run code part needed for parallel processing
    ParParams = TRUE
  }else{ 
    if (ParExec & m==1){
      if(verbose==T){
        message(paste("/!\\ The number of multiple imputation is 1, parallel processing is only available for m > 1."))
      }
    } else if (ParExec){
      if(verbose==T){
        message(paste("/!\\ The number of cores of your processor does not allow paralell processing, at least 3 cores are needed."))
      }
    }
    if(SetRNGSeed){
      set.seed(SetRNGSeed)
    }
    
    foreach::registerDoSEQ()
    opts = NULL
    
    # condition used to run code part needed for sequential processing
    ParParams  = FALSE 
  }
  ndatasets <- length(list_datasets)
  
  #Beginning of the multiple imputation (imputing "m" times)
  RESULT <- foreach(o=1:m, .inorder = TRUE, .options.snow = opts) %dopar% {
    if (!ParParams){
      # Parallel and sequential execution of foreach don't use the same casting mechanism, this one is used for sequential execution.
      cat("iteration :",o,"/",m,"\n")
    }
    
    imputed_datasets <- list()
    imputed_datasets[[1]] <- list()
    # Initial step #
    for(i in 2:ndatasets){
      imputed_datasets[[1]][[i]] <- seqimpute(list_datasets[[i]],np=np,nf=nf,
            m=1,timing=timing,frame.radius=frame.radius, 
            covariates=covariates, time.covariates = time.covariates, regr=regr,
            npt=npt, nfi=nfi, verbose=FALSE,...)$imp[[1]]
      imputed_datasets[[1]][[i]] <-suppressMessages(transform_factor(imputed_datasets[[1]][[i]]))

    }
    
    
    for(n in 2:(niter+1)){
      imputed_datasets[[n]]<-list()
      for(k in 1:ndatasets){
        if(k==1){
          tmp.t.cov <- imputed_datasets[[n-1]][[2]]
          if(ndatasets>2){
            for(i in 3:ndatasets){
              tmp.t.cov <- cbind(tmp.t.cov,imputed_datasets[[n-1]][[i]])
            }
          }
        }else{
          tmp.t.cov <- imputed_datasets[[n]][[1]]
          
          for(i in 2:ndatasets){
            if(i<k){
              tmp.t.cov <- cbind(tmp.t.cov,imputed_datasets[[n]][[i]])
            }else if(i>k){
              tmp.t.cov <- cbind(tmp.t.cov,imputed_datasets[[n-1]][[i]])
            }else{}
          }
        }
        imputed_datasets[[n]][[k]] <- seqimpute(list_datasets[[k]],np=np,nf=nf,
                                                m=1,timing=timing,frame.radius=frame.radius, 
                                                covariates=covariates, time.covariates = tmp.t.cov, regr=regr,
                                                npt=npt, nfi=nfi, verbose=FALSE,...)$imp[[1]]
        imputed_datasets[[n]][[k]] <-suppressMessages(transform_factor(imputed_datasets[[n]][[k]]))
        
      }
    }
    # 
    return(imputed_datasets)
  }
  
  
  if (ParParams){
    parallel::stopCluster(cl)
  }
  return(RESULT)
}

transform_factor <- function(data){
  alphabet <-  attr(seqdef(data),"alphabet")
  for(j in 1:ncol(data)){
    data[,j] <- factor(data[,j],levels=alphabet)
  }
  return(data)
}

