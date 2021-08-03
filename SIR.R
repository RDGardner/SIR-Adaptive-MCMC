# Ricca Callis
# Computational Statistics
# EN.625.664.82.SP21
# Course Project
# 4/28/21
library(deSolve)
library(adaptivetau)
library(ggplot2)
library(plyr)
library(tmvtnorm)
library(parallel)
library(doParallel)
library(stringr)
library(lattice)
library(dplyr)
library(truncnorm)
library(emdbook)
library(reshape2)
# Set  working directory
setwd("/Users/riccacallis/Desktop/JHU/Data Science/Computational Statistics/Project")
load(file="callisProjectData.RData")

# Create simulation function
simulateSIR<-function(theta,initalState,times) {
  # Run SIR mode simulation
  modelSIR <- function(time, state, parameters) {
    # 2 State-transition Rates
    beta <- parameters[["R0"]] / parameters[["D"]] # Infectious rate; D = duration of infection
    gamma <- 1 / parameters[["D"]] # Recovery rate
    # 3 Compartments
    S <- state[["S"]] # Susceptible
    I <- state[["I"]] # Infectious
    R <- state[["R"]] # Recovered
    # Population = N
    # Note this is a closed, homogeneous population
    N <- S + I + R
    # 3 Orindary Differential Equations
    dS <- -beta * S * I/N
    dI <- beta * S * I/N - gamma * I
    dR <- gamma * I
    return(list(c(dS, dI, dR)))
  } # End modelSIR()
  # Link model to data
  modelTrajectory <- data.frame(ode(y = initalState,times = times,func = modelSIR, parms = theta,method = "ode45"))
  return(modelTrajectory)
} # End simulateSIR

# KNOWN PARAMETERS: SIR SIMULATION

# theta = {R0, D}
# Set values for R0 (reproductive rate) & D (duration of infection)
# Remember, beta=R0/(N*D) 
# Assume D (duration of infection) is 1 week
theta <-c(R0=3, D=1)
# Initially, only 1 infected individual & 2999 susceptible individuals
initalState<-c(S=2999, I=1, R = 0)
# Discrete time intervals (1 through 100)
times <-1:100
trajectoryOfSimulatedSIR <- simulateSIR(theta, initalState, times)
head(trajectoryOfSimulatedSIR)
# time         S          I           R
# 1    1 2999.0000   1.000000    0.000000
# 2    2 2989.4546   7.357474    3.187943
# 3    3 2921.1745  52.532364   26.293138
# 4    4 2515.7558 308.538572  175.705609
# 5    5 1405.7142 836.552484  757.733295
# 6    6  590.4092 784.372497 1625.218326

# Create function to plot the trajectory
plotModelTrajectory <- function(trajectoryOfSimulatedSIR = NULL, namesOfStates = NULL, data = NULL, timeCol = "time", linesData = FALSE, summary = TRUE, replicateCol = "replicate", nonExtinct = NULL, alpha = 1, plot = TRUE, colour = "red", same = FALSE) {
  if(!is.null(trajectoryOfSimulatedSIR) & !any(duplicated(trajectoryOfSimulatedSIR[[timeCol]]))) {
    trajectoryOfSimulatedSIR[[replicateCol]] <- 1
    if(summary) {
      # Only 1 replicate to summarise: mean, median and CI of
      # the trajectories won't be plotted.
      summary <- FALSE
    } # End nested if
  } # End outer if
  if(is.null(namesOfStates)) {
    numeric.names <- names(trajectoryOfSimulatedSIR)[sapply(names(trajectoryOfSimulatedSIR), function(x) {
      any(class(trajectoryOfSimulatedSIR[[x]]) %in% c("numeric", "integer"))
    })]
    namesOfStates <- setdiff(numeric.names, c(timeCol, replicateCol))
  } # End if
  if (!is.null(trajectoryOfSimulatedSIR)) {
    if(!is.null(nonExtinct)) {
      trajectoryOfSimulatedSIR <- mutate(trajectoryOfSimulatedSIR, infected = eval(parse(text = paste(nonExtinct, collapse = "+")), trajectoryOfSimulatedSIR))
      df.infected <- melt(trajectoryOfSimulatedSIR, measure.vars = "infected", variable.name = "state")
      df.p.ext <- ddply(df.infected, timeCol, function(df) {
        return(data.frame(value = sum(df$value == 0)/nrow(df)))
      })
      df.p.ext$state <- "p.extinction"
      df.p.ext[replicateCol] <- 0
      if(summary) {
        trajectoryOfSimulatedSIR <- subset(trajectoryOfSimulatedSIR, infected>0)
        trajectoryOfSimulatedSIR$infected <- NULL
      } # End nested if
    } # End outer if
    df.trajectoryOfSimulatedSIR <- melt(trajectoryOfSimulatedSIR, measure.vars = namesOfStates, variable.name = "state")
    df.trajectoryOfSimulatedSIR <- subset(df.trajectoryOfSimulatedSIR, !is.na(value))
    if (summary) {
      message("Compute confidence intervals")
      trajectoryOfSimulatedSIR.CI <- ddply(df.trajectoryOfSimulatedSIR, c(timeCol, "state"), function(df) {
        tmp <- as.data.frame(t(quantile(df$value, prob = c(0.025, 0.25, 0.5, 0.75, 0.975))))
        names(tmp) <- c("low_95", "low_50", "median", "up_50", "up_95")
        tmp$mean <- mean(df$value)
        return(tmp)
      }, .progress = "text")
      trajectoryOfSimulatedSIR.CI.line <- melt(trajectoryOfSimulatedSIR.CI[c(timeCol, "state", "mean", "median")], id.vars = c(timeCol, "state"))
      trajectoryOfSimulatedSIR.CI.area <- melt(trajectoryOfSimulatedSIR.CI[c(timeCol, "state", "low_95", "low_50", "up_50", "up_95")], id.vars = c(timeCol, "state"))
      trajectoryOfSimulatedSIR.CI.area$type <- sapply(trajectoryOfSimulatedSIR.CI.area$variable, function(x) {str_split(x, "_")[[1]][1]})
      trajectoryOfSimulatedSIR.CI.area$CI <- sapply(trajectoryOfSimulatedSIR.CI.area$variable, function(x) {str_split(x, "_")[[1]][2]})
      trajectoryOfSimulatedSIR.CI.area$variable <- NULL
      trajectoryOfSimulatedSIR.CI.area <- dcast(trajectoryOfSimulatedSIR.CI.area, paste0(timeCol,"+state+CI~type"))
      p <- ggplot(trajectoryOfSimulatedSIR.CI.area)
      if (!same) {
        p <- p + facet_wrap(~state, scales = "free_y")
      } # End if
      if (is.null(colour)) {
        p <- p + geom_ribbon(data = trajectoryOfSimulatedSIR.CI.area, aes_string(x = timeCol, ymin = "low", ymax = "up", alpha = "CI"))
        p <- p + geom_line(data = trajectoryOfSimulatedSIR.CI.line, aes_string(x = timeCol, y = "value", linetype = "variable"))
      } # End if
      else if (colour == "all") {
        p <- p + geom_ribbon(data = trajectoryOfSimulatedSIR.CI.area, aes_string(x = timeCol, ymin = "low", ymax = "up", alpha = "CI", fill = "state"))
        p <- p + geom_line(data = trajectoryOfSimulatedSIR.CI.line, aes_string(x = timeCol, y = "value", linetype = "variable", colour = "state"))
      } # End else if
      else {
        p <- p + geom_ribbon(data = trajectoryOfSimulatedSIR.CI.area, aes_string(x = timeCol, ymin = "low", ymax = "up", alpha = "CI"), fill = colour)
        p <- p + geom_line(data = trajectoryOfSimulatedSIR.CI.line, aes_string(x = timeCol, y = "value", linetype = "variable"), colour = colour)                
      } # End else
      p <- p + scale_alpha_manual("Percentile", values = c("95" = 0.25, "50" = 0.45), labels = c("95" = "95th", "50" = "50th"))
      p <- p + scale_linetype("Stats")
      p <- p + guides(linetype = guide_legend(order = 1))
    } # End if(summary) 
    else {
      p <- ggplot(df.trajectoryOfSimulatedSIR)
      if (!same) {
        p <- p + facet_wrap(~state, scales = "free_y")
      } # End nested if
      if (is.null(colour)) {
        if (same) {
          p <- p + geom_line(data = df.trajectoryOfSimulatedSIR, aes_string(x = timeCol, y = "value", group = "state", color = "state"), alpha = alpha)
        } # End nested if 
        else {
          p <- p + geom_line(data = df.trajectoryOfSimulatedSIR, aes_string(x = timeCol, y = "value", group = replicateCol), alpha = alpha)
        } # End nested else
      } # End if 
      else if (colour == "all") {
        p <- p + geom_line(data = df.trajectoryOfSimulatedSIR, aes_string(x = timeCol, y = "value", group = replicateCol, color = replicateCol), alpha = alpha)
      } # End elif
      else {
        p <- p + geom_line(data = df.trajectoryOfSimulatedSIR, aes_string(x = timeCol, y = "value", group = replicateCol), alpha = alpha, colour = colour)
      } # End else
    } # End else
    if(!is.null(nonExtinct)) {
      p <- p+geom_line(data = df.p.ext, aes_string(x = timeCol, y = "value"), color = "black", alpha = 1)
    } # End if
  } # End if 
  else {
    p <- ggplot()
  } # End else
  if(!is.null(data)) {
    obs_names <- grep("obs", names(data), value = TRUE)
    if (length(obs_names) == 0) {
      obs_names <- setdiff(names(data), timeCol)
    } # End if
    data <- melt(data, measure.vars = obs_names, variable.name = "state")
    if (linesData) {
      p <- p + geom_line(data = data, aes_string(x = timeCol, y = "value"), colour = "black")
    } # End if 
    else {
      p <- p + geom_point(data = data, aes_string(x = timeCol, y = "value"), colour = "black")
    } # End else
  } # End if
  p <- p + theme_bw() + theme(legend.position = "top", legend.box = "horizontal")
  if(plot) {
    print(p)
  } # End if
  else{
    return(p)
  } # End else
} # End plotModelTrajectory
# Plot the Model Trajectory of each state (S, I, R)
plotModelTrajectory(trajectoryOfSimulatedSIR)
# Evaluate the Prior for parameters (theta)
# Create the density function
# Use uniform priors
priorDensity<-function(theta, log = FALSE) {
  # Prior for each parameter
  R0logPrior <- dunif(theta[["R0"]], min = 1, max = 100, log = TRUE) # U[1,100]
  DlogPrior <- dunif(theta[["D"]], min = 0, max = 10, log = TRUE) # U[0,10]
  logSUM <- R0logPrior + DlogPrior
  # Return the prior density distribution
  return(ifelse(log, logSUM, exp(logSUM)))
} # End priorDensity()
# Find prior for theta = {R0, D}
priorDensity(theta)
# [1] 0.001010101
# Find log(prior) for theta = {R0, D}
priorDensity(theta, log=TRUE)
# [1] -6.897705

# Evaluate the likelihood
# Data -> Trajectory
# Create the density function to determine probability of observing number of cases given current prevalence
pointLikelihoodDensity<-function(dataObservationPoint, modelObservationPoint, theta, log = FALSE){
  # Poisson prevalence centered around "I"
  return(dpois(x = dataObservationPoint[["obs"]],lambda = modelObservationPoint[["I"]], log = log))
} # End pointLikelihoodDensity()
# Test the new likelihood function
# Calculate the likelihood of an observed prevalence of 20 infections, given a current prevalence of 40 infections
pointLikelihoodDensity(dataObservationPoint=c(obs=20), modelObservationPoint=c(I=40), theta, log=TRUE)
# [1] -8.558027

# RUN SIMULATION USING A DATA SET

# Load data for SIRdata1
#SIRdata1 <- read.csv('SIRdata1.csv', header = FALSE)
head(SIRdata1)
#time obs
# 1    1   0
# 2    2   1
# 3    3   3
# 4    4   1
# 5    5   5
# 6    6   9
View(SIRdata1)
colnames(SIRdata1)
# [1] "time" "obs"
# Plot the model trajectory: number of observations over time
plotModelTrajectory(SIRdata1)
# Create function for log(likelihood)
trajectoryObservationsDensity <- function(theta, initalState, data, log = FALSE) {
  # Time sequence (must include initial time)
  times <- c(0,data$time)
  # Simulate model at successive observation times of data (over time)
  trajectoryOfSimulatedSIR <- simulateSIR(theta,initalState,times)
  dens <- 0
  # Compute log(likelihood) by summing the log(likelihood) of ea. data point
  for(i in 1:nrow(data)){
    # Extract data point
    dataObservationPoint <- unlist(data[i,])
    # Extract model point
    # Use i+1 since the first row of trajectoryOfSimulatedSIR contains the initial state
    modelObservationPoint <- unlist(trajectoryOfSimulatedSIR[i+1,])
    # Update marginal log-likelihood
    dens <- dens + pointLikelihoodDensity(dataObservationPoint=dataObservationPoint, modelObservationPoint=modelObservationPoint, theta=theta, log = TRUE)
  } # End for()
  # Return log(p(Data|X_0, theta))
  return(ifelse(log, dens, exp(dens))) # log(likelihood) of chosen parm & initial state
} # End trajectoryObservationsDensity()
# Run the Simulation
# Evaluate Data -> Model Trajectory
trajectoryObservationsDensity(theta, initalState, SIRdata1, log = TRUE)
# [1] -13709.12

# Calculate posterior
# Model Trajectory -> Data
# i.e., take model and generate a randomly sampled data point
pointObservationRandomSample <-function(modelObservationPoint, theta){
  # Poisson prevalence
  observationPoint <- rpois(n = 1, lambda = modelObservationPoint[["I"]])
  return(c(obs = observationPoint))
} # Ends pointObservationRandomSample
# Given true prevalence of 30, generate a randomly observed prevalence
pointObservationRandomSample(modelObservationPoint=c(I=30), theta)
# obs 
# 32

# Instead of generating a single point, calculate the entire posterior/trajectory
# Create the function which will take the pointObservationRandomSample at each time point to generate an entire trajectory
trajectoryObservationsRandomSample <- function(theta, initalState, times) {
  # Simulate model at successive observation times of data
  trajectoryOfSimulatedSIR <- simulateSIR(theta, initalState, times)
  # Generate observations by
  # PointObservation to each row of trajectoryOfSimulatedSIR
  obs <- ddply(trajectoryOfSimulatedSIR, "time" , pointObservationRandomSample, theta = theta)
  trajectoryObservations <- join(trajectoryOfSimulatedSIR,obs, by="time")
  return(trajectoryObservations)
} # End trajectoryObservationsRandomSample()
# Generate the entire trajectory of observations
entireTrajectoryOfSimulatedSIR<-trajectoryObservationsRandomSample(theta, initalState, SIRdata1$time)
head(entireTrajectoryOfSimulatedSIR)
# time         S          I           R obs
# 1    1 2999.0000   1.000000    0.000000   3
# 2    2 2989.4546   7.357474    3.187943   3
# 3    3 2921.1745  52.532364   26.293138  57
# 4    4 2515.7558 308.538572  175.705609 309
# 5    5 1405.7142 836.552484  757.733295 788
# 6    6  590.4092 784.372497 1625.218326 807

# CALCULATE THE POSTERIOR

# Posterior = p(theta, X_0|Data)
logPosteriorDensity <- function(theta, initalState, data) {
  # Calculate the prior
  logPrior <- priorDensity(theta, log = TRUE)
  # Calculate log(likelihood): (p(data}theta))
  # Calculate log(likelihood) of theta
  # Calculate initalState w.r.t the data using trajectoryObservationsDensity
  # Assign to logLikelihood    
  logLikelihood <- trajectoryObservationsDensity(theta, initalState, data, log = TRUE)
  # Calulate log(posterior)
  # Posterior = likelihood*prior
  logPosterior <- logPrior + logLikelihood
  return(logPosterior)
} # End logPosteriorDensity
# Calculate
logPosteriorDensity(theta, initalState, SIRdata1)
# [1] -13716.02

# ASSESS MODEL FIT

# Assess model by running simulation multiple times
# Create function to run multiple simulations
multipleModelSimulations <- function(theta, initalState, times, n, observation=FALSE) {
  stopifnot(n>0)
  if(observation && is.null(pointLikelihoodDensity)){
    stop("Can't generate observation as model doesn't have a ",sQuote("pointLikelihoodDensity")," function.")
  } # End if
  rep <- as.list(1:n)
  names(rep) <- rep
  if (n > 1) {
    progress = "text"
  } # End if
  else {
    progress = "none"
  } # End else
  # Simulate the model under theta & generate observations
  trajectoryOfSimulatedSIR.rep <- ldply(rep,function(x) {
    if(observation){
      trajectoryOfSimulatedSIR <- trajectoryObservationsRandomSample(theta, initalState, times)
    } # End if (observation)
    else {
      trajectoryOfSimulatedSIR <- simulateSIR(theta,initalState,times)
    } # End else
    # Return trajectory
    return(trajectoryOfSimulatedSIR)
  },.progress=progress,.id="replicate")
  # Return model repplicates
  return(trajectoryOfSimulatedSIR.rep)
} # End multipleModelSimulations()

# Plot multiple model simulations against the data
plotModelFit <- function(theta, initalState, data, numberOfReplicates = 1, summary = TRUE, alpha = min(1, 10/numberOfReplicates), allVariables = FALSE, nonExtinct = NULL, observation = TRUE, plot = TRUE) {
  # Can use deterministic model trajectories (using allVariables = TRUE)
  # Simulate multiple replicates (using the numberOfReplicates argument)
  times <- c(0, data$time)
  if (numberOfReplicates > 1) {
    cat("Simulate ", numberOfReplicates, " replicate(s)\n")
  }
  trajectoryOfSimulatedSIR <- multipleModelSimulations(theta = theta, initalState = initalState, times = times, n = numberOfReplicates, observation = observation)
  if(allVariables) {
    namesOfStates <- NULL
  } 
  else {
    namesOfStates <- grep("obs", names(trajectoryOfSimulatedSIR), value = TRUE)
  }
  p <- plotModelTrajectory(trajectoryOfSimulatedSIR = trajectoryOfSimulatedSIR, namesOfStates = namesOfStates, data = data, summary = summary, alpha = alpha, nonExtinct = nonExtinct, plot = FALSE)
  if(plot) {
    print(p)
  } 
  else {
    return(list(trajectoryOfSimulatedSIR = trajectoryOfSimulatedSIR, plot = p))
  }
} # End plotModelFit()

# Assess Model Fit with Plot
# Using the model,generate an observation trajectory from a model using trajectoryObservationsRandomSample
# Plot the observation trajectory (lines) against the data (points)
plotModelFit(theta, initalState, SIRdata1)

# PARAMETER ESTIMATION: R0
# Explore posterior distribution at different values of the parameter R0 (skip)

# MCMC
# Sample from the posterior distribution using MCMC
# Given an initial parameter, evaluate target distribution
# Use Metropolis-Hastings
# Sample from univariate distribution & use a standard Gaussian proposal distribution q(θ′|θ)

# Metropolis-Hastings Sampler
# Univariate
# proposalStDv = sd of Gaussian proposal distribution
# numberOfIterations = sampler runs
univariateMCMCMH <- function(target, initialTheta, proposalStDv, numberOfIterations) {
  # Evaluate target distribution at any parameter value
  targetThetaCurrent <- target(initialTheta)
  # Initialise variables
  # Set current value of theta
  currentTheta <- initialTheta
  # Give current theta to samples
  samples <- currentTheta
  # Set number of accepted runs
  accepted <- 0
  # Run MCMC for each i in sequence
  for (ithIteration in seq_len(numberOfIterations)) {
    # Randomly sample new theta from Gaussian proposal distribution
    # sd = avg(step size of sampler)
    proposedTheta <- rnorm(n = length(currentTheta),mean = currentTheta,sd = proposalStDv)
    # B/c 'rnorm' returns an unnamed vector, set names of proposedTheta = the names of currentTheta
    names(proposedTheta) <- names(currentTheta)
    # Evaluate target at the proposed theta
    targetProposedTheta <- target(proposedTheta)
    # Calculate acceptance
    logAcceptance <- targetProposedTheta - targetThetaCurrent
    # Random sample from [0,1] 
    r <- runif(1)
    # Test acceptance
    # Compare r to p(acceptance)
    # Exponentiate b/c log ratios
    if (r < exp(logAcceptance)) {
      # If accepted:
      # Change current value to proposed value of theta
      currentTheta <- proposedTheta
      # Updated current value of target
      targetThetaCurrent <- targetProposedTheta
      # Increment number of accepted
      accepted <- accepted + 1
    } # End if accepted
    # Add current theta to samples vector
    # rbind for later multivariate use
    samples <- rbind(samples, currentTheta, deparse.level=0)
    # Print update: iteration number, current state of chain, acceptance rate
    message("iteration: ", ithIteration, ", chain:", paste(currentTheta, collapse=" "),
            ", acceptance rate:", accepted / ithIteration)
  } # End for each number of iterations
  # Return the trace of the chain
  # trace = sample theta vector from target distribution
  return(samples)
} # End univariateMCMCMH()
# Check that this new sampler is working as expected
# Use normal distribution
plot(dnorm,xlim = c(-6, 6), ylab = "Probability Density")
# Take the log
logNormalDensity <- function(theta) {
  return(dnorm(x = theta, mean = 0, sd = 1, log = TRUE))
} # End logNormalDensity()
# Test MCMC sampler by sampling from logNormalDensity
# Do randomly generated numbers look like same distribution
startingValue <- 1 
sigma <- 1 # sd of MCMC
iter <- 1000
traceUniMCMC <- univariateMCMCMH(target = logNormalDensity, initialTheta = startingValue,
                   proposalStDv = sigma, numberOfIterations = iter)
# Plot trace of MCMC
plot(traceUniMCMC, type = "l")
# Plot a histogram
hist(traceUniMCMC, freq = FALSE)
curve(dnorm, from = -5, to = 5, col="red", add=TRUE)

# Sample from a posterior distribution
# Use given R0
# Fixed initialState (X0)
logPosteriorDensityR0 <- function(R0) {
  # Create wrapper function: take previous function and change it so it accepts only one parameter
  # Keep infection duration fixed (~ 2 weeks)
  # Return log(posterior)
  return(logPosteriorDensity(theta = c(R0 = R0, D = 1),
                          initalState = c(S = 2999, I = 1, R = 0),
                          data = SIRdata1))
} # End logPosteriorDensityR0
# Run sample from posterior
# Given a value of R0, return a value of posterior
logPosteriorDensityR0(R0 = 1)
# [1] -5177.548

# Run MCMC
startingValue <- 1 
sigma <- 1 # sd of MCMC
iter <- 1000
traceUniMCMC <- univariateMCMCMH(target = logPosteriorDensityR0, initialTheta = startingValue,
                            proposalStDv = sigma, numberOfIterations = iter)
# Plot trace of MCMC
plot(traceUniMCMC, type = "l")
# Makes sense b/c R0 = 1, so disease remains relatively constant
# Descriptive statistics
# Sample mean of R0
mean(traceUniMCMC)
# [1] 1.35343
# Sample median
median(traceUniMCMC)
# [1] 1.346838
# 95% cCI
quantile(traceUniMCMC, probs=c(0.025, 0.975))
# 2.5%    97.5% 
#   1.346838 1.379958 

# Multivariate
plotModelTrajectory(SIRdata3)
# Create log(posterior) fitted to new dataset
# initial conditions of 2999 susceptibles, 1 infected, and 0 recovered.
logPosteriorDensity_SIRdata3 <- function(theta) {
  # Wrapper function to change number of arguments
  return(logPosteriorDensity(theta = theta,initalState = c(S = 2999, I = 1, R = 0),data = SIRdata3))
} # End logPosteriorDensity_SIRdata3()
# Write MCMC function with two unknown parameters
multivariateMCMCMH <- function(target, initialTheta, proposalStDv = NULL,numberOfIterations, covmat = NULL,
                   limits=list(lower = NULL, upper = NULL),startSizeAdaptation = NULL, slowAdaptations = 0.99,
                   startShapeAdaptation = NULL, stopShapeAdaptation = NULL,printResultsEvery = numberOfIterations/100,
                   verbose = FALSE, maxScalingStDv = 50) {
  # Initialise theta
  currentTheta <- initialTheta
  theta.propose <- initialTheta
  # Extract theta of gaussian proposal
  covmat.proposal <- covmat
  lower.proposal <- limits$lower
  upper.proposal <- limits$upper
  # Reorder vector and matrix by names, set to default if necessary
  theta.names <- names(initialTheta)
  if (!is.null(proposalStDv) && is.null(names(proposalStDv))) {
    names(proposalStDv) <- theta.names
  } # End if 
  if (is.null(covmat.proposal)) {
    if (is.null(proposalStDv)) {
      proposalStDv <- initialTheta/10
    } # End nested if
    covmat.proposal <-
      matrix(diag(proposalStDv[theta.names]^2, nrow = length(theta.names)),
             nrow = length(theta.names),
             dimnames = list(theta.names, theta.names))
  } # End if
  else {
    covmat.proposal <- covmat.proposal[theta.names,theta.names]
  } # End else
  if (is.null(lower.proposal)) {
    lower.proposal <- initialTheta
    lower.proposal[] <- -Inf
  } # End if
  else {
    lower.proposal <- lower.proposal[theta.names]
  } # End else
  if (is.null(upper.proposal)) {
    upper.proposal <- initialTheta
    upper.proposal[] <- Inf
  } # End if
  else {
    upper.proposal <- upper.proposal[theta.names]
  } # End else
  # covmat init
  covmat.proposal.init <- covmat.proposal
  adapting.size <- FALSE # will be set to TRUE once we start
  # Adapting the size
  adapting.shape <- 0  # Will be set to the iteration at which
  # Adaptation starts
  # Find estimated theta
  theta.estimated.names <- names(which(diag(covmat.proposal) > 0))
  # Evaluate target at theta init
  targetThetaCurrent <- target(currentTheta)
  if (!is.null(printResultsEvery)) {
    message(Sys.time(), ", Init: ", printNamedVector(currentTheta[theta.estimated.names]),
            ", target: ", targetThetaCurrent)
  } # End if
  # Trace
  trace <- matrix(ncol=length(currentTheta)+1, nrow=numberOfIterations, 0)
  colnames(trace) <- c(theta.estimated.names, "log.density")
  # Acceptance rate
  acceptance.rate <- 0
  # Scaling factor for covmat size
  scaling.sd  <- 1
  # Scaling multiplier
  scaling.multiplier <- 1
  # Empirical covariance matrix (0 everywhere initially)
  covmat.empirical <- covmat.proposal
  covmat.empirical[,] <- 0
  # Empirical mean vector
  theta.mean <- currentTheta
  # If printResultsEvery is null never print info
  if (is.null(printResultsEvery)) {
    printResultsEvery <- numberOfIterations + 1
  } # End if
  start_iteration_time <- Sys.time()
  for (ithIteration in seq_len(numberOfIterations)) {
    # Adaptive step
    if (!is.null(startSizeAdaptation) && ithIteration >= startSizeAdaptation &&
        (is.null(startShapeAdaptation) || acceptance.rate*ithIteration < startShapeAdaptation)) {
      if (!adapting.size) {
        message("\n---> Start adapting size of covariance matrix")
        adapting.size <- TRUE
      } # End if
      # Adapt size of covmat until we get enough accepted jumps
      scaling.multiplier <- exp(slowAdaptations^(ithIteration-startSizeAdaptation) * (acceptance.rate - 0.234))
      scaling.sd <- scaling.sd * scaling.multiplier
      scaling.sd <- min(c(scaling.sd,maxScalingStDv))
      # Only scale if it doesn't reduce the covariance matrix to 0
      covmat.proposal.new <- scaling.sd^2*covmat.proposal.init
      if (!(any(diag(covmat.proposal.new)[theta.estimated.names] <
                .Machine$double.eps))) {
        covmat.proposal <- covmat.proposal.new
      } # End if
    } # End if
    else if (!is.null(startShapeAdaptation) &&
               acceptance.rate*ithIteration >= startShapeAdaptation &&
               (adapting.shape == 0 || is.null(stopShapeAdaptation) ||
                ithIteration < adapting.shape + stopShapeAdaptation)) {
      if (!adapting.shape) {
        message("\n---> Start adapting shape of covariance matrix")
        # Flush.console()
        adapting.shape <- ithIteration
      } # End if
      # Adapt shape of covmat using optimal scaling factor for multivariate target distributions
      scaling.sd <- 2.38/sqrt(length(theta.estimated.names))
      covmat.proposal <- scaling.sd^2 * covmat.empirical
    } # End elif 
    else if (adapting.shape > 0) {
      message("\n---> Stop adapting shape of covariance matrix")
      adapting.shape <- -1
    } # End elif
    # Print info
    if (ithIteration %% ceiling(printResultsEvery) == 0) {
      message(Sys.time(), ", Iteration: ",ithIteration,"/", numberOfIterations,
              ", acceptance rate: ",
              sprintf("%.3f",acceptance.rate), appendLF=FALSE)
      if (!is.null(startSizeAdaptation) || !is.null(startShapeAdaptation)) {
        message(", scaling.sd: ", sprintf("%.3f", scaling.sd),
                ", scaling.multiplier: ", sprintf("%.3f", scaling.multiplier),
                appendLF=FALSE)
      } # End nested if
      message(", state: ",(printNamedVector(currentTheta)))
      message(", logdensity: ", targetThetaCurrent)
    } # End if
    # Propose another parameter set
    if (any(diag(covmat.proposal)[theta.estimated.names] <
            .Machine$double.eps)) {
      print(covmat.proposal[theta.estimated.names,theta.estimated.names])
      stop("non-positive definite covmat",call.=FALSE)
    } # End if
    if (length(theta.estimated.names) > 0) {
      theta.propose[theta.estimated.names] <-
        as.vector(rtmvnorm(1,
                           mean =
                             currentTheta[theta.estimated.names],
                           sigma =
                             covmat.proposal[theta.estimated.names,theta.estimated.names],
                           lower =
                             lower.proposal[theta.estimated.names],
                           upper = upper.proposal[theta.estimated.names]))
    } # End if
    # Evaluate posterior of proposed parameter
    target.theta.propose <- target(theta.propose)
    # If return value is a vector, set log.density and trace
    if (!is.finite(target.theta.propose)) {
      # If posterior is 0 then do not compute anything else and don't accept
      logAcceptance <- -Inf
    } # End if
    else{
      # Compute acceptance
      logAcceptance <- target.theta.propose - targetThetaCurrent
      logAcceptance <- logAcceptance +
        dtmvnorm(x = currentTheta[theta.estimated.names],
                 mean =
                   theta.propose[theta.estimated.names],
                 sigma =
                   covmat.proposal[theta.estimated.names,
                                   theta.estimated.names],
                 lower =
                   lower.proposal[theta.estimated.names],
                 upper =
                   upper.proposal[theta.estimated.names],
                 log = TRUE)
      logAcceptance <- logAcceptance -
        dtmvnorm(x = theta.propose[theta.estimated.names],
                 mean = currentTheta[theta.estimated.names],
                 sigma =
                   covmat.proposal[theta.estimated.names,
                                   theta.estimated.names],
                 lower =
                   lower.proposal[theta.estimated.names],
                 upper =
                   upper.proposal[theta.estimated.names],
                 log = TRUE)
    } # End else
    if (verbose) {
      message("Propose: ", theta.propose[theta.estimated.names],
              ", target: ", target.theta.propose,
              ", acc prob: ", exp(logAcceptance), ", ",
              appendLF = FALSE)
    } # End if
    if (is.accepted <- (log(runif (1)) < logAcceptance)) {
      # accept proposed parameter set
      currentTheta <- theta.propose
      targetThetaCurrent <- target.theta.propose
      if (verbose) {
        message("accepted")
      } # End nested if
    } # End if
    else if (verbose) {
      message("rejected")
    } # End elif
    trace[ithIteration, ] <- c(currentTheta, targetThetaCurrent)
    # Update acceptance rate
    if (ithIteration == 1) {
      acceptance.rate <- is.accepted
    } # End if
    else {
      acceptance.rate <- acceptance.rate +
        (is.accepted - acceptance.rate) / ithIteration
    } # End else
    # Update empirical covariance matrix
    if (adapting.shape >= 0) {
      tmp <- updateCovmat(covmat.empirical, theta.mean,
                          currentTheta, ithIteration)
      covmat.empirical <- tmp$covmat
      theta.mean <- tmp$theta.mean
    } # End if
  } # End for ithIteration
  return(list(trace = trace,
              acceptance.rate = acceptance.rate,
              covmat.empirical = covmat.empirical))
} # End multivariateMCMCMH()
# Run it
mcmcSIRdata3 <- multivariateMCMCMH(target = logPosteriorDensity_SIRdata3,
                       initialTheta = c(R0 = 2, D = 1),
                       proposalStDv = c(0.01, 0.1),
                       numberOfIterations = 1000)
traceD3<-mcmcSIRdata3$trace
head(traceD3)
# R0        D log.density
# [1,] 2.000000 1.000000  -14179.738
# [2,] 2.000000 1.000000  -14179.738
# [3,] 1.996765 1.115162  -11088.636
# [4,] 1.992155 1.199804   -9212.302
# [5,] 1.992155 1.199804   -9212.302
# [6,] 1.992155 1.199804   -9212.302

# Use coda; convert to recognizable format
library('coda')
mcmcTraceD3 <- mcmc(traceD3)
class(mcmcTraceD3)
# [1] "mcmc"
View(mcmcTraceD3)
# Summary statistics
summary(mcmcTraceD3)
# Iterations = 1:1000
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 1000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
# naive se: adjusts for sample size
# time-series se: corrects naive se for autocorrelations
#   
#   Mean        SD Naive SE Time-series SE
# R0             1.57    0.1458  0.00461        0.09452
# D              1.60    0.2154  0.00681        0.11698
# log.density -400.41 1056.4395 33.40755      188.19369
# 
# 2. Quantiles for each variable:
#   
#   2.5%      25%      50%      75%    97.5%
# R0              1.477    1.494    1.498    1.531    1.975
# D               1.463    1.489    1.504    1.541    2.185
# log.density -2298.474 -157.689 -142.285 -141.462 -141.049

# Calculate acceptance rate
acceptanceRate <- 1 - rejectionRate(mcmcTraceD3)
acceptanceRate
# R0           D log.density 
# 0.1321321   0.1321321   0.1321321 

# Calculate effective sample size
# Estimate number of independent samples generated by MCMC
effectiveSize(mcmcTraceD3)
# R0           D log.density 
# 2.378953    3.389006   31.512230 

# MIXING
# Ensure MCMC sampler explores parameter space efficiently 
# Dont want sampler to accept or reject too many

# Trace Plots
# Goal: Avoid flats (where chain stays in same state for a long time)
plot(mcmcTraceD3)
# Here, we see burn-in ~ 1000 iterations, then poor mixing
# R0: Sampler never moves past 2
# D (duration infection): Sampler never moves past 3.5
# Possible we need to start w/higher initial values of each

# Autocorrelations: Check for convergence
# lag-k: correlation b/t every sample & the sample k steps prior
# Expect: smaller as k increases (i.e., independent samples)
# If higher for high values of k, may indicate correlation b/t samples &/or high mixing
autocorr.plot(mcmcTraceD3)
# Here, we see autocorrelation decreases with increasing k

# Adaptive MCMC
# Automate process of adjusting proposal distributions
# Steps:
# Initialize: startSizeAdaptation, startShapeAdaptation and slowAdaptations
# Run MCMC & monitor acceptance rate
# After reaching the number of iterations (startSizeAdaptation), start adapting size of the steps
# Here, we scale the proposal distribution to smaller steps if large acceptance or larger steps if small acceptance
# When startShapeAdaptation proposals are accepted,
# Get empirical covariance matrix
# Adapt shape of proposal 
# Slow down adaptation as slowAdaptations approaches 1.

# Run Adaptive MCMC
traceAMCMC <- multivariateMCMCMH(target = logPosteriorDensity_SIRdata3,
                initialTheta = c(R0 = 2, D = 1),
                proposalStDv = c(1, 0.5),
                numberOfIterations = 5000,
                startSizeAdaptation = 100,
                startShapeAdaptation = 500,
                slowAdaptations=0.999,
                limits = list(lower = c(R0 = 0, D = 0)))
trace2AMCMC<-traceAMCMC$trace
head(trace2AMCMC)
# R0        D log.density
# [1,] 2.000000 1.000000 -14179.7378
# [2,] 2.000000 1.000000 -14179.7378
# [3,] 1.654035 2.125152   -824.8449
# [4,] 1.654035 2.125152   -824.8449
# [5,] 1.654035 2.125152   -824.8449
# [6,] 1.654035 2.125152   -824.8449
mcmcTrace2AMCMC<-mcmc(trace2AMCMC)
summary(mcmcTrace2AMCMC)
# Iterations = 1:5000
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 5000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean        SD  Naive SE Time-series SE
# R0             1.492   0.01838 0.0002600      0.0009477
# D              1.489   0.03757 0.0005313      0.0021791
# log.density -149.465 282.58874 3.9964083      7.2773390
# 
# 2. Quantiles for each variable:
#   
#   2.5%      25%      50%      75%    97.5%
# R0             1.474    1.487    1.490    1.496    1.507
# D              1.447    1.477    1.488    1.497    1.525
# log.density -148.978 -142.763 -141.616 -141.293 -140.985
plot(mcmcTrace2AMCMC)
# Compare acceptance rate at the beginning of adaptation vs
# shape rtoward the end
