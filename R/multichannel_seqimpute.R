#' Imputation of missing data in multichannel sequence data
#' 
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom seqimpute seqimpute
#' @importFrom TraMineR seqdef
#'
#' @description This function implements the MICT-multichannel algorithms. This 
#' algorithm extends the MICT and MICT-timing algorithm to the case of 
#' multichannel sequences.
#' 
#'
#' @param channels A list containing either a data frame containing sequences of a categorical
#' variable, where missing data are coded as \code{NA}, or a state sequence
#' object created using the \link[TraMineR]{seqdef} function. If using a
#' state sequence object, any "void" elements will also be treated as missing.
#' See the \code{end.impute} argument if you wish to skip imputing values
#' at the end of the sequences.
#' @param np Number of prior states to include in the imputation model
#' for internal gaps.
#' @param nf Number of subsequent states to include in the imputation model
#' for internal gaps.
#' @param m Number of multiple imputations to perform (default: \code{5}).
#' @param niter Number of iterations of the algorithm.
#' @param timing Logical, specifies the imputation algorithm to use.
#' If \code{FALSE}, the MICT algorithm is applied; if \code{TRUE}, the
#' MICT-timing algorithm is used.
#' @param frame.radius Integer, relevant only for the MICT-timing algorithm,
#' specifying the radius of the timeframe.
#'
#' @param covariates List of the columns of the dataset
#' containing covariates to be included in the imputation model.
#'
#' @param time.covariates List of the columns of the dataset
#'  with time-varying covariates to include in the imputation model.
#'
#' @param regr Character specifying the imputation method. Options include
#' \code{"multinom"} for multinomial models and \code{"rf"} for random forest
#' models.
#'
#' @param npt Number of prior observations in the imputation model for
#' terminal gaps (i.e., gaps at the end of sequences).
#'
#' @param nfi Number of future observations in the imputation model for
#' initial gaps (i.e., gaps at the beginning of sequences).
#'
#' @param ParExec Logical, indicating whether to run multiple imputations
#' in parallel. Setting to \code{TRUE} can improve computation time depending
#' on available cores.
#'
#' @param ncores Integer, specifying the number of cores to use for parallel
#' computation. If unset, defaults to the maximum number of CPU cores minus one.
#'
#' @param SetRNGSeed Integer, to set the random seed for reproducibility in
#' parallel computations. Note that setting \code{set.seed()} alone does not
#' ensure reproducibility in parallel mode.
#'
#'
#' @param verbose Logical, if \code{TRUE}, displays progress and warnings
#' in the console. Use \code{FALSE} for silent computation.
#'
#' @param available Logical, specifies whether to consider already imputed
#' data in the predictive model. If \code{TRUE}, previous imputations are
#' used; if \code{FALSE}, only original data are considered.
#'
#' @param pastDistrib Logical, if \code{TRUE}, includes the past distribution
#' as a predictor in the imputation model.
#'
#' @param futureDistrib Logical, if \code{TRUE}, includes the future
#' distribution as a predictor in the imputation model.
#' @param ... Named arguments that are passed down to the imputation functions.
#'
#' @author Kevin Emery <kevin.emery@@unige.ch>, Andre Berchtold,
#' Anthony Guinchard, and Kamyar Taher
#'
#' @return An object of class \code{seqimp}, which is a list with the following
#' elements:
#' \describe{
#'   \item{\code{data}}{A \code{data.frame} containing the original
#'   (incomplete) data.}
#'   \item{\code{imp}}{A list of \code{m} \code{data.frame} corresponding to
#'   the imputed datasets.}
#'   \item{\code{m}}{The number of imputations.}
#'   \item{\code{method}}{A character vector specifying whether MICT or
#'   MICT-timing was used.}
#'   \item{\code{np}}{Number of prior states included in the imputation model.}
#'   \item{\code{nf}}{Number of subsequent states included in the imputation
#'   model.}
#'   \item{\code{regr}}{A character vector specifying whether multinomial or
#'   random forest imputation models were applied.}
#'   \item{\code{call}}{The call that created the object.}
#' }
#'
#' @references Halpin, B. (2012). Multiple imputation for life-course
#' sequence data. Working Paper WP2012-01, Department of Sociology,
#' University of Limerick. http://hdl.handle.net/10344/3639.
#' @references Halpin, B. (2013). Imputing sequence data: Extensions to
#' initial and terminal gaps, Stata's. Working Paper WP2013-01,
#' Department of Sociology,
#' University of Limerick. http://hdl.handle.net/10344/3620
#' @references Emery, K., Studer, M., & Berchtold, A. (2024). Comparison of
#' imputation methods for univariate categorical longitudinal data.
#' Quality & Quantity, 1-25.
#' https://link.springer.com/article/10.1007/s11135-024-02028-z
#'
#' @export
seqimputemc <- function(channels, np=1, nf=1, m=5, niter=1, 
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
  ndatasets <- length(channels)
  
  #Beginning of the multiple imputation (imputing "m" times)
  o <- NULL
  RESULT <- foreach(o=1:m, .inorder = TRUE, .options.snow = opts) %dopar% {
    if (!ParParams){
      # Parallel and sequential execution of foreach don't use the same casting mechanism, this one is used for sequential execution.
      cat("iteration :",o,"/",m,"\n")
    }
    
    imputed_datasets <- list()
    imputed_datasets[[1]] <- list()
    # Initial step #
    for(i in 2:ndatasets){
      imputed_datasets[[1]][[i]] <- seqimpute(channels[[i]], var=NULL, np=np,nf=nf,
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
        imputed_datasets[[n]][[k]] <- seqimpute(channels[[k]],var=NULL,np=np,nf=nf,
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

