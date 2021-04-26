function [theta_logL, recon_f_state_var] = get_LogL_alleleDynamics_wilker(theta, data, MCMC_params, C)

init_freqs = theta(1:data.n_ferrets);
MOI = exp(theta((data.n_ferrets + 1)));
fitness = exp(theta((data.n_ferrets + 2)));

wmean = [];
for i = 1:MCMC_params.n_particles
    particles(i).gen_list = 0;
    particles(i).f_list = init_freqs';
    particles(i).f_curr = init_freqs';
    particles(i).w_list = [];
end

for g = 0:max(data.generation)

    if g ~= 0
        % simulate one passage
        for i = 1:MCMC_params.n_particles
            particles(i).gen_list = [particles(i).gen_list g];
            particles(i) = simulate_one_generation_wilker(MCMC_params, particles(i), fitness, C, data.n_ferrets);
        end
    end
    
    if ismember(g, data.generation)
        loc = find(data.generation == g);
        for i = 1:MCMC_params.n_particles
            for ferret = 1:data.n_ferrets
                data_f = data.f_generation(ferret, loc);
                % probability of data given model
                if particles(i).f_curr(ferret,1) < MCMC_params.threshold_f_curr
                    particles(i).f_curr(ferret,1) = MCMC_params.threshold_f_curr;
                end
                if particles(i).f_curr(ferret,1) > (1-MCMC_params.threshold_f_curr)
                    particles(i).f_curr(ferret,1) = 1-MCMC_params.threshold_f_curr;
                end
                A = MCMC_params.nu_noise*particles(i).f_curr(ferret,1);
                B = MCMC_params.nu_noise*(1-particles(i).f_curr(ferret,1));
                this_w(ferret,1) = betapdf(data_f,A,B);
            end
            particles(i).w_list = [particles(i).w_list prod(this_w)];
            w_vector(i,1) = prod(this_w);
        end
        k_vector = fast_resample(w_vector, MCMC_params);
        wmean = [wmean mean(w_vector)];
        for i = 1:MCMC_params.n_particles
            newparticles(i) = particles(k_vector(i));
        end
        particles = newparticles;
    end
    
end

theta_logL = sum(log(wmean));

rand_val = randi(MCMC_params.n_particles);
recon_f_state_var = particles(rand_val).f_list;
