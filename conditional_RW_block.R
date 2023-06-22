sampler_conditional_RW_block <- nimbleFunction(
    name = 'sampler_conditional_RW_block',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptScaleOnly      <- extractControlElement(control, 'adaptScaleOnly',      FALSE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scale               <- extractControlElement(control, 'scale',               1)
        propCov             <- extractControlElement(control, 'propCov',             'identity')
        tries               <- extractControlElement(control, 'tries',               1)
        ## set incorrect defaults so that it throws an error if not explicitly set
        indices             <- extractControlElement(control, 'indices',               1)
        
        
        ## checks on target
        if(length(target) != 1) stop('length of target not equal to one')
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        d <- length(targetAsScalar)
        if(d <= 1) stop('target does not have more than one dimension')
        
        ## checks on indices
        if(length(indices) < 1 | length(indices) >= d) stop('length of indices must be between 1 and the number of dimensions of target')
        if(!all(is.integer(indices))) stop('indices must be integers')
        for(i in 1:length(indices)) {
            if(indices[i] < 1 | indices[i] > d) stop('indices must be within number of dimensions of target')
        }
        
        ## node list generation
        #targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
        if(!is.integer(finalTargetIndex) |
           length(finalTargetIndex) != 1 |
           is.na(finalTargetIndex[1]))
            stop('Problem with target node in sampler_RW_block')
        calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
        calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
        ##calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesDepStage <- model$isStoch(calcNodesDepStage)   ## should be made faster
        calcNodesDepStageDeterm <- calcNodesDepStage[!isStochCalcNodesDepStage]
        calcNodesDepStageStoch <- calcNodesDepStage[isStochCalcNodesDepStage]
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        #d <- length(targetAsScalar)
        d1 <- length(indices)
        scaleHistory  <- c(0, 0)                                                 ## scaleHistory
        acceptanceHistory  <- c(0, 0)                                            ## scaleHistory
        propCovHistory <- if(d1<=10) array(0, c(2,d1,d1)) else array(0, c(2,2,2))   ## scaleHistory
        saveMCMChistory <- if(nimbleOptions('MCMCsaveHistory')) TRUE else FALSE
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d1)
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        chol_propCov_scale <- scale * chol_propCov
        empirSamp <- matrix(0, nrow=adaptInterval, ncol=d1)
        ## nested function and function list definitions
        ##my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
        targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ##my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)   ## old syntax: missing target argument
        my_calcAdaptationFactor <- calcAdaptationFactor(d1, adaptFactorExponent)
        ## checks
        if(!inherits(propCov, 'matrix'))        stop('propCov must be a matrix\n')
        if(!inherits(propCov[1,1], 'numeric'))  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d1))           stop('propCov matrix must have dimension ', d1, 'x', d1, '\n')
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
    },
    run = function() {
        for(i in 1:tries) {
            currentValues <- model[[target]]
            propValueVector <- model[[target]]
            ## extract current values in order
            currValues <- numeric(d1)
            for(j in 1:d1) {
                currValues[j] <- currentValues[indices[j]]
            }
            propValues <- generateProposalVector(currValues)
            for(j in 1:d1) {
                propValueVector[indices[j]] <- propValues[j]
            }
            ##lpMHR <- my_setAndCalculateDiff$run(propValueVector)
            values(model, targetNodesAsScalar) <<- propValueVector
            lpD <- model$calculateDiff(calcNodesProposalStage)
            if(lpD == -Inf) {
                nimCopy(from = mvSaved, to = model,   row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
            ## Drawing a random number is needed during first testing
            ## of this step in order to keep the random numbers identical
            ## to old behavior to see if tests that depend on particular
            ## sample sequences pass.  Rather than calling runif(1, 0, 1) here,
            ## we call decide() to ensure same behavior.
            ## jump <- decide(lpD)
            ## When new behavior is acceptable, we can remove the above line
            ## and uncomment the following:
            jump <- FALSE
            } else {
                ##jump <- my_decideAndJump$run(lpMHR, 0, 0, 0) ## will use lpMHR - 0
                lpD <- lpD + model$calculateDiff(calcNodesDepStage)
                jump <- decide(lpD)
                if(jump) {
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageDeterm, logProb = FALSE)
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesDepStageStoch, logProbOnly = TRUE)
                } else {
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesDepStageDeterm, logProb = FALSE)
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesDepStageStoch, logProbOnly = TRUE)
                }
            }
            if(adaptive)     adaptiveProcedure(jump)
        }
    },
    methods = list(
        generateProposalVector = function(currValues = double(1)) {
            propValueVector <- rmnorm_chol(1, currValues, chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
            returnType(double(1))
            return(propValueVector)
        },
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(!adaptScaleOnly) {
                ## extract current values in order
                currValues <- numeric(d1)
                for(i in 1:d1) {
                    currValues[i] <- values(model, target)[indices[i]]
                }
                empirSamp[timesRan, 1:d1] <<- currValues
            }
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                if(saveMCMChistory) {
                    setSize(scaleHistory, timesAdapted)                 ## scaleHistory
                    scaleHistory[timesAdapted] <<- scale                ## scaleHistory
                    setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
                    acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
                    if(d1 <= 10) {
                        propCovTemp <- propCovHistory                                           ## scaleHistory
                        setSize(propCovHistory, timesAdapted, d1, d1)                             ## scaleHistory
                        if(timesAdapted > 1)                                                    ## scaleHistory
                            for(iTA in 1:(timesAdapted-1))                                      ## scaleHistory
                                propCovHistory[iTA, 1:d1, 1:d1] <<- propCovTemp[iTA, 1:d1, 1:d1]    ## scaleHistory
                        propCovHistory[timesAdapted, 1:d1, 1:d1] <<- propCov[1:d1, 1:d1]            ## scaleHistory
                    }
                }
                adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
                scale <<- scale * adaptFactor
                ## calculate empirical covariance, and adapt proposal covariance
                if(!adaptScaleOnly) {
                    gamma1 <- my_calcAdaptationFactor$getGamma1()
                    for(i in 1:d1)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
                    empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
                    propCov <<- propCov + gamma1 * (empirCov - propCov)
                    chol_propCov <<- chol(propCov)
                }
                chol_propCov_scale <<- chol_propCov * scale
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        getScaleHistory = function() {  ## scaleHistory
            if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
            returnType(double(1))
            return(scaleHistory)
        },          
        getAcceptanceHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
            return(acceptanceHistory)
        },                  
        getPropCovHistory = function() { ## scaleHistory
            if(!saveMCMChistory | d1 > 10)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC and note that to reduce memory use we only save the proposal covariance history for parameter vectors of length 10 or less")
            returnType(double(3))
            return(propCovHistory)
        },
        reset = function() {
            scale   <<- scaleOriginal
            propCov <<- propCovOriginal
            chol_propCov <<- chol(propCov)
            chol_propCov_scale <<- chol_propCov * scale
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            if(saveMCMChistory) {
                scaleHistory  <<- c(0, 0)    ## scaleHistory
                acceptanceHistory  <<- c(0, 0)
                if(d1 <= 10) 
                    propCovHistory <<- nimArray(0, dim = c(2,d1,d1))
            }
            my_calcAdaptationFactor$reset()
        }
    )
)
