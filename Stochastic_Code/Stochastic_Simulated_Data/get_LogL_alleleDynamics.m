function [theta_logL, recon_f_state_var] = get_LogL_alleleDynamics(theta, data, MCMC_params, C)

init_freq = theta(1);
MOI = exp(theta(2));
fitness = exp(theta(3));

wmean = [];

% initialize particles:
for i = 1:MCMC_params.n_particles
    particles(i).gen_list = 0;
    particles(i).f_list = init_freq;
    particles(i).w_list = [];
    particles(i).f_curr = init_freq;
end

% do Sequential Monte Carlo (SMC):
for g = 1:max(data.generation)
    
    % simulate one generation
    for i = 1:MCMC_params.n_particles
        particles(i).gen_list = [particles(i).gen_list g];
        particles(i) = simulate_one_generation(MCMC_params, particles(i), fitness, C);
    end
    
    % do selection on particles every time there is an observed data point:
    if ismember(g, data.generation)
        loc = find(data.generation == g);
        data_f = data.f_generation(loc);
        for i = 1:MCMC_params.n_particles
            % if particle's variant frequency falls below minimum frequency
            % or above (1-minimum frequency), set to threshold value:
            if particles(i).f_curr < MCMC_params.threshold_f_curr
                particles(i).f_curr = MCMC_params.threshold_f_curr;
            end
            if particles(i).f_curr > (1-MCMC_params.threshold_f_curr)
                particles(i).f_curr = 1-MCMC_params.threshold_f_curr;
            end
            % calculate particle weights:
            A = MCMC_params.nu_noise*particles(i).f_curr;
            B = MCMC_params.nu_noise*(1-particles(i).f_curr);
            this_w = betapdf(data_f,A,B);
            particles(i).w_list = [particles(i).w_list this_w];
            w_vector(i,1) = this_w;
        end
        % resample particles:
        k_vector = fast_resample(w_vector, MCMC_params);
        for i = 1:MCMC_params.n_particles
            newparticles(i) = particles(k_vector(i));
        end
        particles = newparticles;
        % calculate the mean particle weight from this timepoint, and tag
        % on to a vector that keeps this info around:
        wmean = [wmean mean(w_vector)];
    end
end

% this is the marginal (log) likelihood:
theta_logL = sum(log(wmean));

% sample a random particle at the very end and return it as the reconstructed state variable: 
rand_val = randi(MCMC_params.n_particles);
recon_f_state_var = particles(rand_val).f_list;
