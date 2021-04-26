function void = main_plot_pMCMC(void)

clear all; close all; clc;

infile = 'pMCMC_results_coinfection_sim_N100_nu100'

load(infile);

burnin = 5000;
MCMC_samplingfrequency = 100;
num_reconStateVars_plot = 10;
max_sample = 50000;
%max_sample = length(MCMC_results.theta(:,1));

% MCMC iterations following burn-in:
locs = burnin:max_sample; %length(MCMC_results.theta(:,1));
locs_sample = burnin:MCMC_samplingfrequency:max_sample; %length(MCMC_results.theta(:,1));
index_locs = randperm(length(locs_sample), num_reconStateVars_plot);
locs_sample_recon = locs_sample(index_locs);

% RED AS TRUE; Black as inferred (credible intervals, etc.; recon state vars)! blue = MCMC data
figure(5);
subplot(2,2,1);
plot(0:max(data.generation), MCMC_results.recon_f_state_var(locs_sample_recon,:), 'k'); hold on;
plot(data.generation, data.f_generation, 'r.',  'MarkerSize',20); axis([-1 26 0 1]); % plot(0, 0.1, 'mo');  
%simulated_data = [0.1 0.13 0.1 0.14 0.15 0.19 0.19 0.24 0.25 0.28 0.41 0.37 0.37 0.42 0.48 0.48 0.47 0.58 0.68 0.79 0.82 0.87 0.85 0.87 0.93 0.91];
simulated_data = [0.1 0.07 0.1 0.14 0.13 0.13 0.12 0.17 0.18 0.17 0.25 0.3 0.31 0.38 0.39 0.44 0.45 0.48 0.54 0.57 0.64 0.71 0.71 0.69 0.72 0.69];

length(0:max(data.generation))
length(simulated_data)
plot(0:max(data.generation), simulated_data, 'r', 'LineWidth',2);
xlabel('Generation'); ylabel('Variant frequency');
axis([0 27 0 1]);

init_freq_ordered = sort(MCMC_results.theta(locs_sample,1));
n_samples = length(init_freq_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_freq_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_freq_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_freq_ordered(median_loc);

subplot(2,2,2); histogram(MCMC_results.theta(locs_sample,1), 'Normalization', 'probability');
axis([-0.050 0.4 0 0.3]); 
hold on; y = axis; plot([0.1 0.1], [0 y(4)], 'r', 'LineWidth',2); 
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('Initial frequency'); ylabel('Proportion');

init_MOI_ordered = sort(exp(MCMC_results.theta(locs_sample,2)));
n_samples = length(init_MOI_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_MOI_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_MOI_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_MOI_ordered(median_loc);

subplot(2,2,3); histogram(exp(MCMC_results.theta(locs_sample,2)), 'Normalization', 'probability'); axis([-0.050 7 0 0.3]); 
hold on; y = axis; plot([2 2], [0 y(4)], 'r', 'LineWidth',2); 
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('MOI'); ylabel('Proportion'); thisaxis = axis; axis([0 7 thisaxis(3) thisaxis(4)]); 

init_fitness_ordered = sort(exp(MCMC_results.theta(locs_sample,3)));
n_samples = length(init_fitness_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_fitness_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_fitness_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_fitness_ordered(median_loc);

subplot(2,2,4); histogram(exp(MCMC_results.theta(locs_sample,3)), 'Normalization', 'probability'); axis([0.75 3 0 0.3]); 
hold on; y = axis; plot([1.5 1.5], [0 y(4)], 'r', 'LineWidth',2); 
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('Variant fitness'); ylabel('Proportion');

% Supplemental figure 8 (trace plot):
figure(8);
subplot(1,3,1); plot(MCMC_results.theta(1:max_sample,1)); hold on; plot(0.1*ones(size(MCMC_results.theta(1:max_sample,1))), 'r', 'LineWidth',2); xlabel('MCMC iteration'); ylabel('Initial frequency'); y = axis; plot([burnin burnin], [y(3) y(4)], 'k--');
subplot(1,3,2); semilogy(exp(MCMC_results.theta(1:max_sample,2))); hold on; semilogy(2*ones(size(MCMC_results.theta(1:max_sample,1))), 'r', 'LineWidth',2);  xlabel('MCMC iteration'); ylabel('MOI'); y = axis; plot([burnin burnin], [y(3) y(4)], 'k--');
subplot(1,3,3); semilogy(exp(MCMC_results.theta(1:max_sample,3))); hold on; semilogy(1.5*ones(size(MCMC_results.theta(1:max_sample,1))), 'r', 'LineWidth',2);  xlabel('MCMC iteration'); ylabel('Variant fitness'); y = axis; plot([burnin burnin], [y(3) y(4)], 'k--');
