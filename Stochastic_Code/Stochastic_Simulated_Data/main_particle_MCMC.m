function void = main_particle_MCMC(void)

clear all; close all; clc;

outfile = 'pMCMC_results_coinfection_sim_N100_nu100';

% set up MCMC parameters
MCMC_params.n_particles = 200;   % number of particles to use in particle filtering algorithm
MCMC_params.iterations = 50000;  % number of MCMC steps
MCMC_params.write_steps = 100;  % frequency at which MCMC chain is written out

data.generation = [5 10 15 20 25];
data.f_generation = [0.094 0.270 0.447, 0.546, 0.699];

MCMC_params.nu_noise = 100;
MCMC_params.Nvirions = 100;
MCMC_params.threshold_f_curr = 0.001; % if the variant allele becomes lost or fixed in a simulation (f = 0 or 1, respectively), beta distribution yields -Inf, causing a problem for the likelihood calculation
% this threshold value changes fixed values of 0 to the threshold value and fixed values of 1 to (1-threshold value)

% parameters to estimate in pMCMC, and initial guesses; true values are f0 = 0.1, MOI = 2, fitness (exp^sigma_m) = 1.5
% estimate log(MOI) and log(fitness) to accomodate positivity constraint
MCMC_params.init_f0 = 0.3;          % true: 0.1
MCMC_params.init_logMOI = log(3);   % true: log(2)
MCMC_params.init_sigmam = log(2);   % true: log(1.5)

% prior on log(MOI)
MCMC_params.logMOI_prior_mu = log(2);
MCMC_params.logMOI_prior_sigma = 0.5;

% variances for proposing new set of parameters (theta)
MCMC_params.f0_variance = 0.0001;
MCMC_params.log_MOI_variance = 0.02;
MCMC_params.sigmam_variance = 0.002;

run_pMCMC(data, MCMC_params, outfile);
