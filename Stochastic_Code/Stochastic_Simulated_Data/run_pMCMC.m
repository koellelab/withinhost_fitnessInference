function void = run_pMCMC(data, MCMC_params, outfile)

if (min(data.f_generation) <= 0) || (max(data.f_generation) >= 1)
    error('code cannot handle variants that are lost or fixed, or carry values that are negative or above one.')
end

% theta holds [init freq, log(MOI), and log(fitness) = sigma_m]
curr_theta = [MCMC_params.init_f0 MCMC_params.init_logMOI MCMC_params.init_sigmam];
% no covariances implemented:
MCMC_params.theta_covar = [MCMC_params.f0_variance 0 0; 0 MCMC_params.log_MOI_variance 0; 0 0 MCMC_params.sigmam_variance];

if 0
    % checking what the proposals would look like:
    vals = mvnrnd(curr_theta, MCMC_params.theta_covar,1000);
    subplot(2,2,1); hist(vals(:,1)); xlabel('f0 proposal'); ylabel('frequency');
    subplot(2,2,2); hist(exp(vals(:,2))); xlabel('MOI proposal'); ylabel('frequency');
    subplot(2,2,3); hist(exp(vals(:,3))); xlabel('fitness proposal'); ylabel('frequency');
    x = -4:0.001:4; y = normpdf(x, MCMC_params.logMOI_prior_mu, MCMC_params.logMOI_prior_sigma)
    subplot(2,2,4); plot(exp(x), y); thisaxis = axis; axis([0 6 thisaxis(3) thisaxis(4)]);
    xlabel('MOI (prior)'); ylabel('density');
    return;
end

C = round(MCMC_params.Nvirions/exp(MCMC_params.init_logMOI));
[curr_logL, curr_recon_f_state_var] = get_LogL_alleleDynamics(curr_theta, data, MCMC_params, C);
MCMC_results.MCMC_iteration = 1; MCMC_results.theta = curr_theta;
MCMC_results.logL = curr_logL; MCMC_results.recon_f_state_var = curr_recon_f_state_var;
curr_priorVal = lognpdf(exp(curr_theta(2)), MCMC_params.logMOI_prior_mu, MCMC_params.logMOI_prior_sigma);
MCMC_results.recon_f_state_var(1,:) = curr_recon_f_state_var;
    
for m = 2:MCMC_params.iterations

    % propose another set of parameters: mvnrnd takes in as arguments means and covariance matrix
    while 1
        proposal_theta = mvnrnd(curr_theta, MCMC_params.theta_covar);
        % check to make sure that f0 falls in range (0, 1)
        if (proposal_theta(1) > 0) && (proposal_theta(1) < 1)
            break;
        end
        display('proposed f0 needs to be positive-- looping');
    end
           
    % evaluate loglikelihood of proposed theta:
    C = round(MCMC_params.Nvirions/exp(proposal_theta(2)));
    [proposal_logL, proposed_recon_f_state_var] = get_LogL_alleleDynamics(proposal_theta, data, MCMC_params, C);
    
    % in the absence of a prior, accept the proposed theta if its likelihood if greater than that of the current theta's likelihood, or accept proposed theta with some probability if lower than current theta's likelihood:
    %prob_accept = exp(proposal_logL - curr_logL);  % this is with no prior
        
    % with prior on MOI:
    proposal_priorVal = lognpdf(exp(proposal_theta(2)), MCMC_params.logMOI_prior_mu, MCMC_params.logMOI_prior_sigma);
    p_theta_prop = exp(proposal_logL)*proposal_priorVal;
    p_theta_curr = exp(curr_logL)*curr_priorVal;
    prob_accept = min(1, (p_theta_prop/p_theta_curr));
    
    if prob_accept > rand
        curr_logL = proposal_logL;
        curr_theta = proposal_theta;
        curr_recon_f_state_var = proposed_recon_f_state_var;
        curr_priorVal = proposal_priorVal;
    end
    
    MCMC_results.MCMC_iteration = [MCMC_results.MCMC_iteration; m];
    MCMC_results.logL = [MCMC_results.logL; curr_logL];
    MCMC_results.theta = [MCMC_results.theta; curr_theta];
    MCMC_results.recon_f_state_var = [MCMC_results.recon_f_state_var; curr_recon_f_state_var];
        
    if mod(m, MCMC_params.write_steps) == 0
        m
        save(outfile, 'data', 'MCMC_params', 'MCMC_results');
    end
end
