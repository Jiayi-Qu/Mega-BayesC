library(data.table)
library(rrBLUP)
library(tidyr)
library(lme4qtl)
library(lme4)
library(MegaLMM)


source('Estimate_gcor_prediction.R')
set_MegaLMM_nthreads(RcppParallel::defaultNumThreads()-1)

foldid = as.numeric(commandArgs(t=T)[1])
if(is.na(foldid)) foldid = 1
set.seed(foldid)

source('data_prep_Krause.R')


# foldid = as.numeric(commandArgs(t=T)[1])
# if(is.na(foldid)) foldid = 1
# set.seed(foldid)

runID = sprintf('MegaLMM_Krause_BayesC_L0.9_Mest_nofixfactor_%d',foldid)

  run_parameters = MegaLMM_control(
    which_sampler = list(Y=4,F=4),
    run_sampler_times = 1,
    drop0_tol = 1e-10,
    scale_Y = F,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
    h2_divisions = 2, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
    h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
    burn = 2000,  # number of burn in samples before saving posterior samples
    K = 100 # number of factors
  )

  Y = cbind(data$yNA,HTP_wide[match(data$GID,HTP_wide$GID),-1])
  Y = scale(Y)

  priors = MegaLMM_priors(
    #tot_Y_var = list(V = 0.5,   nu = 10),      # Prior variance of trait residuals after accounting for fixed effects and factors
    tot_F_var = list(V = 18/20, nu = 10000),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
    tot_Y_var = list(V = .5*10/11,nu = 11),
    Lambda_prior = list(
      sampler = sample_Lambda_prec_BayesC, # function that implements the horseshoe-based Lambda prior described in Runcie et al 2020. See code to see requirements for this function.
      Lambda_df = 4,
      Lambda_scale = 1,
      delta_1 = list(shape = 20, rate = 1/2),
      delta_2 = list(shape = 3, rate = 1),
      fixed_pi = 0.9,
      delta_iterations_factor = 100   # parameter that affects mixing of the MCMC sampler. This value is generally fine.
    ),
    B2_prior = list(
      sampler = sample_B2_prec_BayesC,
      fixed_pi = NULL,
      B2_F_df = 4,
      B2_F_scale = 1,
      B2_R_df = 4,
      B2_R_scale = 1
    ),
    h2_priors_resids_fun = function(h2s,n) ifelse(h2s == 0,1,0),
    h2_priors_factors_fun = function(h2s,n) ifelse(h2s == 0,1,0)
    )

  MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                                ~ 1 + (1|GID),
                                extra_regressions = list(X = X,resids = F,factors=T),
                                data=data,         # the data.frame with information for constructing the model matrices
                                Lambda_fixed = matrix(c(1,rep(0,ncol(Y)-1)),nr=1),
                                run_parameters=run_parameters,
                                run_ID = runID
  )
  maps = make_Missing_data_map(MegaLMM_state,2,verbose=T)
  MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map)
  MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
  MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
  MegaLMM_state_base = MegaLMM_state
  MegaLMM_state = initialize_MegaLMM(MegaLMM_state)
  MegaLMM_state = clear_Posterior(MegaLMM_state)

  
  MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F','tot_F_prec','F_h2','tot_Eta_prec','resid_h2','B2_F','Lambda_pi','B2_F_pi')#,'Eta_mean')
  MegaLMM_state = clear_Posterior(MegaLMM_state)

  # set_MegaLMM_nthreads(7)
  n_iter = 100;  # how many samples to collect at once?
  
  # system(sprintf('rm %s/U_pred_samples.csv',runID))
  for(i  in 1:100) {
    print(sprintf('Run %d',i))
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)

    MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
    print(MegaLMM_state) # print status of current chain
    plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf

    # set of commands to run during burn-in period to help chain converge
    if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i <= 20) {
      MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
      MegaLMM_state = clear_Posterior(MegaLMM_state)
      print(MegaLMM_state$run_parameters$burn)
    }
  }

  rm(MegaLMM_state)
  gc()

 
  MegaLMM_state = MegaLMM_state_base
  MegaLMM_state$Posterior = readRDS(sprintf('%s/Posterior/Posterior_base.rds',runID))

  MegaLMM_state$Posterior$Lambda = load_posterior_param(MegaLMM_state,'Lambda')
  MegaLMM_state$Posterior$B2_F = load_posterior_param(MegaLMM_state,'B2_F')
  U = get_posterior_mean(MegaLMM_state,X2_F %*% B2_F %*% Lambda[,1])
  
  MegaLMM_state$Posterior$resid_h2 = load_posterior_param(MegaLMM_state,'resid_h2')
  MegaLMM_state$Posterior$F_h2 = load_posterior_param(MegaLMM_state,'F_h2')
  MegaLMM_state$Posterior$tot_Eta_prec = load_posterior_param(MegaLMM_state,'tot_Eta_prec')
  MegaLMM_state$Posterior$tot_F_prec = load_posterior_param(MegaLMM_state,'tot_F_prec')
  MegaLMM_state$Posterior$F = load_posterior_param(MegaLMM_state,'F')
  
  write.csv(get_posterior_mean(MegaLMM_state, Lambda), sprintf('%s/Posterior/Posterior_Lambda.csv',runID))
  write.csv(get_posterior_mean(MegaLMM_state, B2_F), sprintf('%s/Posterior/Posterior_B2_F.csv',runID))
  write.csv(get_posterior_mean(MegaLMM_state, resid_h2), sprintf('%s/Posterior/Posterior_resid_h2.csv',runID))
  write.csv(get_posterior_mean(MegaLMM_state, F_h2), sprintf('%s/Posterior/Posterior_F_h2.csv',runID))
  write.csv(get_posterior_mean(MegaLMM_state, tot_Eta_prec), sprintf('%s/Posterior/Posterior_tot_Eta_prec.csv',runID))
  write.csv(get_posterior_mean(MegaLMM_state, tot_F_prec), sprintf('%s/Posterior/Posterior_tot_F_prec.csv',runID))
  write.csv(get_posterior_mean(MegaLMM_state, F), sprintf('%s/Posterior/Posterior_F.csv',runID))

  print("csv saved.")
  results = data.frame(Method = 'MegaLMM_BayesC_L0.9_Mest_nofixfactor',
                       pearson = cor(data$BLUP[nas],U[nas,1])/sqrt(h2_BLUE),
                       g_cor = estimate_gcor(data.frame(ID=data$GID[nas],obs = data$BLUP[nas],pred = U[nas,1]),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']])
  results$fold = foldid
  write.csv(results,file = sprintf('%s/results_MegaLMM_BayesC_L0.9_Mest_fold_%d.csv',results_dir,foldid))
  
  
  
  
 
  
  


