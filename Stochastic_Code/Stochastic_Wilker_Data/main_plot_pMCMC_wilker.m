function void = main_plot_pMCMC_wilker(void)

clear all; close all; clc;

infile = 'pMCMC_results_coinfection_Wilker_N100_nu100';

load(infile);

burnin = 5000;
MCMC_samplingfrequency = 100;
num_reconStateVars_plot = 10;
max_sample = 50000;

% MCMC iterations following burn-in:
locs = burnin:max_sample; %length(MCMC_results.theta(:,1));
locs_sample = burnin:MCMC_samplingfrequency:max_sample; %length(MCMC_results.theta(:,1));
index_locs = randperm(length(locs_sample), num_reconStateVars_plot);
locs_sample_recon = locs_sample(index_locs);

% Supplemental figure 9 (trace plot):
figure(9);
subplot(2,3,1); plot(MCMC_results.theta(1:max_sample,1)); hold on; xlabel('MCMC iteration'); ylabel('Initial frequency'); y = axis; plot([burnin burnin], [y(3) y(4)], 'k--');
subplot(2,3,2); plot(MCMC_results.theta(1:max_sample,2)); hold on; xlabel('MCMC iteration'); ylabel('Initial frequency'); y = axis; plot([burnin burnin], [y(3) y(4)], 'k--');
subplot(2,3,3); plot(MCMC_results.theta(1:max_sample,3)); hold on; xlabel('MCMC iteration'); ylabel('Initial frequency'); y = axis; plot([burnin burnin], [y(3) y(4)], 'k--');
subplot(2,3,4); plot(MCMC_results.theta(1:max_sample,4)); hold on; xlabel('MCMC iteration'); ylabel('Initial frequency'); y = axis; plot([burnin burnin], [y(3) y(4)], 'k--');
subplot(2,3,5); semilogy(exp(MCMC_results.theta(1:max_sample,5))); hold on; xlabel('MCMC iteration'); ylabel('MOI'); y = axis; plot([burnin burnin], [y(3) y(4)], 'k--');
subplot(2,3,6); semilogy(exp(MCMC_results.theta(1:max_sample,6))); hold on; xlabel('MCMC iteration'); ylabel('Variant fitness'); y = axis; plot([burnin burnin], [y(3) y(4)], 'k--');
%subplot(2,2,4); plot(MCMC_results.logL); ylabel('log-likelihood'); y = axis; hold on; plot([burnin burnin], [y(3) y(4)], 'k--');

% Supplemental figure 10 (joint density plot):
figure(10);
plot(exp(MCMC_results.theta(locs_sample,(data.n_ferrets + 2))), exp(MCMC_results.theta(locs_sample,(data.n_ferrets + 1))), 'b.'); ylabel('MOI'); xlabel('Variant fitness'); 
ax = gca; ax.YDir = 'reverse'; 
axis([0 20 0 10]);

% Supplemental figure 11 (initial frequency posteriors):
figure(11);

init_freq_ordered = sort(MCMC_results.theta(locs_sample,1));
n_samples = length(init_freq_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_freq_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_freq_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_freq_ordered(median_loc);

subplot(2,2,1); histogram(MCMC_results.theta(locs_sample,1), 'Normalization', 'probability'); hold on; axis([0 0.15 0 0.25]);
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('Initial frequency (ferret 13)'); ylabel('Proportion');

init_freq_ordered = sort(MCMC_results.theta(locs_sample,2));
n_samples = length(init_freq_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_freq_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_freq_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_freq_ordered(median_loc);

subplot(2,2,2); histogram(MCMC_results.theta(locs_sample,2), 'Normalization', 'probability'); hold on; axis([0 0.15 0 0.25]); 
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('Initial frequency (ferret 15)'); ylabel('Proportion');

init_freq_ordered = sort(MCMC_results.theta(locs_sample,3));
n_samples = length(init_freq_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_freq_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_freq_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_freq_ordered(median_loc);

subplot(2,2,3); histogram(MCMC_results.theta(locs_sample,3), 'Normalization', 'probability'); hold on; axis([0 0.15 0 0.25]);
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('Initial frequency (ferret 17)'); ylabel('Proportion');

init_freq_ordered = sort(MCMC_results.theta(locs_sample,4));
n_samples = length(init_freq_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_freq_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_freq_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_freq_ordered(median_loc);

subplot(2,2,4); histogram(MCMC_results.theta(locs_sample,4), 'Normalization', 'probability'); hold on; axis([0 0.15 0 0.25]); 
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('Initial frequency (ferret 21)'); ylabel('Proportion');

% main manuscript figure
figure(6);
init_MOI_ordered = sort(exp(MCMC_results.theta(locs_sample,5)));
n_samples = length(init_MOI_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_MOI_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_MOI_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_MOI_ordered(median_loc);
vals_MOI = [CI_low med_val CI_high]

subplot(1,2,1); histogram(exp(MCMC_results.theta(locs_sample,5)), 'Normalization', 'probability'); hold on; axis([0 10 0 0.20]); %y = axis;
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('MOI'); ylabel('Proportion'); thisaxis = axis; axis([0 10 thisaxis(3) thisaxis(4)]); 

init_fitness_ordered = sort(exp(MCMC_results.theta(locs_sample,6)));
n_samples = length(init_fitness_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_fitness_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_fitness_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_fitness_ordered(median_loc);

subplot(1,2,2); histogram(exp(MCMC_results.theta(locs_sample,6)), 'Normalization', 'probability'); hold on; axis([0 10 0 0.25]); %y = axis;
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('Variant fitness'); ylabel('Proportion');
axis([0 20 y(3) y(4)])
vals_fitness = [CI_low med_val CI_high]

clear all;

infile = 'pMCMC_results_coinfection_Wilker_N40_nu100';

load(infile);

burnin = 5000;
MCMC_samplingfrequency = 100;
num_reconStateVars_plot = 10;
max_sample = 50000;
%max_sample = length(MCMC_results.theta(:,1));

figure(12);

% MCMC iterations following burn-in:
locs = burnin:max_sample; 
locs_sample = burnin:MCMC_samplingfrequency:max_sample;
index_locs = randperm(length(locs_sample), num_reconStateVars_plot);
locs_sample_recon = locs_sample(index_locs);

init_MOI_ordered = sort(exp(MCMC_results.theta(locs_sample,5)));
n_samples = length(init_MOI_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_MOI_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_MOI_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_MOI_ordered(median_loc);
vals_MOI = [CI_low med_val CI_high]

subplot(1,2,1); histogram(exp(MCMC_results.theta(locs_sample,5)), 'Normalization', 'probability'); hold on; y = axis;
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('MOI'); ylabel('Proportion'); thisaxis = axis; axis([0 10 thisaxis(3) thisaxis(4)]); 

init_fitness_ordered = sort(exp(MCMC_results.theta(locs_sample,6)));
n_samples = length(init_fitness_ordered);
CI_low_loc = ceil(0.025*n_samples); CI_low = init_fitness_ordered(CI_low_loc);
CI_high_loc = floor(0.975*n_samples); CI_high = init_fitness_ordered(CI_high_loc);
median_loc = round(0.5*n_samples); med_val = init_fitness_ordered(median_loc);

subplot(1,2,2); histogram(exp(MCMC_results.theta(locs_sample,6)), 'Normalization', 'probability'); hold on; y = axis;
y = axis; plot([CI_low CI_low], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([CI_high CI_high], [0 y(4)], 'k--', 'LineWidth',2);
y = axis; plot([med_val med_val], [0 y(4)], 'k', 'LineWidth',2);
xlabel('Variant fitness'); ylabel('Proportion');
axis([0 20 y(3) y(4)])
vals_fitness = [CI_low med_val CI_high]

return;

figure; 
for this_f = 1:data.n_ferrets
    subplot(4,1,this_f); plot(0:max(data.generation), MCMC_results.recon_f_state_var(this_f).this_ferret(locs_sample_recon,:), 'r')
    hold on; plot([0 data.generation], [0.044 data.f_generation(this_f,:)], 'k'); plot([0 data.generation], [0.044 data.f_generation(this_f,:)], 'ko'); plot(0, 0.044, 'bo');  axis([-1 10 0 1]);
    xlabel('generation'); ylabel('mutant frequency');
end