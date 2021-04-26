function void = main_particle_MCMC_wilker(void)

clear all; close all; clc;

outfile = 'pMCMC_results_coinfection_Wilker_N100_nu100';

% set up MCMC parameters
MCMC_params.n_particles = 200;   % number of particles to use in particle filtering algorithm
MCMC_params.iterations = 50000;  % number of MCMC steps
MCMC_params.write_steps = 100;   % frequency at which MCMC chain is written out

data.time_days = [0 1 3 5];
data.f_time_days_ferret13 = [0.044 0.057 0.320 0.720];
data.f_time_days_ferret15 = [0.044 0.188 0.449 0.683];
data.f_time_days_ferret17 = [0.044 0.353 0.632 0.809];
data.f_time_days_ferret21 = [0.044 0.202 0.721 0.590];

data.t_generation = 8;      % 8 hrs = 1 generation

% convert to generations, and only use day 1 and day 3 data:
data.generation = (data.time_days(1:3))*24/data.t_generation;
data.f_generation(1,:) = data.f_time_days_ferret13(1:3);
data.f_generation(2,:) = data.f_time_days_ferret15(1:3);
data.f_generation(3,:) = data.f_time_days_ferret17(1:3);
data.f_generation(4,:) = data.f_time_days_ferret21(1:3);
data.n_ferrets = 4;

MCMC_params.nu_noise = 100;
MCMC_params.Nvirions = 100;
MCMC_params.threshold_f_curr = 0.001; % if simulations yield a fixed value (f = 0 or 1), beta distribution yields -Inf, causing a problem for the likelihood calculation
% this threshold value changes fixed values of 0 to the threshold value and fixed values of 1 to (1-threshold value)

% parameters to estimate in pMCMC, and initial guesses
MCMC_params.init_f0 = 0.08*ones(4,1); % stock value, for each ferret
MCMC_params.init_logMOI = log(3);   
MCMC_params.init_sigmam = log(2);   

% prior on log(MOI)
MCMC_params.logMOI_prior_mu = log(4);
MCMC_params.logMOI_prior_sigma = 0.4;

% variances for proposing new theta
% estimate f0, log(MOI), and log(fitness) (that is, sigma_m) to ensure
% positivity for MOI and fitness
MCMC_params.f0_variance = 0.0001;
MCMC_params.log_MOI_variance = 0.02;
MCMC_params.sigmam_variance = 0.002;

run_pMCMC_wilker(data, MCMC_params, outfile);
