function void = run_pMCMC_wilker(data, MCMC_params, outfile)

if (min(min(data.f_generation)) == 0) || (max(max(data.f_generation)) == 1)
    error('code cannot handle mutants that are lost or fixed')
end

% theta holds [init freq, log(MOI), and sigmam]
curr_theta = [MCMC_params.init_f0' MCMC_params.init_logMOI MCMC_params.init_sigmam];
MCMC_params.theta_covar = zeros(length(curr_theta));

for i = 1:data.n_ferrets
    MCMC_params.theta_covar(i,i) = MCMC_params.f0_variance;
end
MCMC_params.theta_covar((data.n_ferrets + 1), (data.n_ferrets + 1)) = MCMC_params.log_MOI_variance;
MCMC_params.theta_covar((data.n_ferrets + 2), (data.n_ferrets + 2)) = MCMC_params.sigmam_variance;

if 0
    % checking what the proposals would look like:
    vals = mvnrnd(curr_theta, MCMC_params.theta_covar,1000);
    subplot(2,2,1); hist(vals(:,1)); xlabel('f0 proposal'); ylabel('frequency');
    subplot(2,2,2); hist(exp(vals(:,(data.n_ferrets + 1)))); xlabel('MOI proposal'); ylabel('frequency');
    subplot(2,2,3); hist(exp(vals(:,(data.n_ferrets + 2)))); xlabel('fitness proposal'); ylabel('frequency');
    x = -4:0.001:4; y = normpdf(x, MCMC_params.logMOI_prior_mu, MCMC_params.logMOI_prior_sigma)
    subplot(2,2,4); 
    figure; plot(exp(x), y); thisaxis = axis; axis([0 10 thisaxis(3) thisaxis(4)]);
    xlabel('MOI (prior)'); ylabel('density');
    return;
end

C = round(MCMC_params.Nvirions/exp(MCMC_params.init_logMOI));
[curr_logL, curr_recon_f_state_var] = get_LogL_alleleDynamics_wilker(curr_theta, data, MCMC_params, C);
MCMC_results.MCMC_iteration = 1; MCMC_results.theta = curr_theta;
MCMC_results.logL = curr_logL; 
for i = 1:data.n_ferrets
    MCMC_results.recon_f_state_var(i).this_ferret(1,:) = curr_recon_f_state_var(i,:);
end
curr_priorVal_logMOI = lognpdf(exp(curr_theta(data.n_ferrets + 1)), MCMC_params.logMOI_prior_mu, MCMC_params.logMOI_prior_sigma);
curr_priorVal_f0 = 1;
curr_priorVal = curr_priorVal_f0*curr_priorVal_logMOI;
    
for m = 2:MCMC_params.iterations
    
    % propose another set of parameters: mvnrnd takes in as arguments means and covariance matrix
    while 1
        proposal_theta = mvnrnd(curr_theta, MCMC_params.theta_covar);
        if min(proposal_theta(1:data.n_ferrets)) > 0
            break;
        end
        display('proposed f0 for each ferret needs to be positive-- looping');
    end
           
    % evaluate loglikelihood of proposed theta:
    C = round(MCMC_params.Nvirions/exp(proposal_theta(data.n_ferrets + 1)));
    [proposal_logL, proposed_recon_f_state_var] = get_LogL_alleleDynamics_wilker(proposal_theta, data, MCMC_params, C);
    
    proposal_priorVal_logMOI = lognpdf(exp(proposal_theta(data.n_ferrets + 1)), MCMC_params.logMOI_prior_mu, MCMC_params.logMOI_prior_sigma);

    proposal_priorVal_f0 = 1;
    proposal_priorVal = proposal_priorVal_f0*proposal_priorVal_logMOI;
    
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
    for this_f = 1:data.n_ferrets
        MCMC_results.recon_f_state_var(this_f).this_ferret = [MCMC_results.recon_f_state_var(this_f).this_ferret; curr_recon_f_state_var(this_f,:)];
    end
        
    if mod(m, MCMC_params.write_steps) == 0
        m
        curr_logL
        save(outfile, 'data', 'MCMC_params', 'MCMC_results');
    end
    
end
