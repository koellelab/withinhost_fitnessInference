function particle = simulate_one_generation_wilker(MCMC_params, particle, fitness_m, C, n_ferrets)

% manuscript text needs to be made consistent with this (eqns 8 and 9) 

for ferret = 1:n_ferrets
    
    nMutantVirions = round(MCMC_params.Nvirions*particle.f_curr(ferret));
    nWTVirions = MCMC_params.Nvirions - nMutantVirions;
    
    p = (1/C)*ones(C, 1);

    r_mutant = mnrnd(nMutantVirions, p);
    r_wt = mnrnd(nWTVirions, p);
 
    sum_f_m = 0;
    sum_f_wt = 0;
    for i = 1:C
        F = GetF(r_mutant(i), r_wt(i), fitness_m);
        sum_f_m = sum_f_m + r_mutant(i)*F;
        sum_f_wt = sum_f_wt + r_wt(i)*F;
    end
    mean_fitness_mutant = sum_f_m/nMutantVirions;
    mean_fitness_wt = sum_f_wt/nWTVirions;
    f_next = (particle.f_curr(ferret,1)*mean_fitness_mutant)/(particle.f_curr(ferret,1)*mean_fitness_mutant + (1-particle.f_curr(ferret,1))*mean_fitness_wt);
    f_next = binornd(MCMC_params.Nvirions, f_next)/MCMC_params.Nvirions;
    
    particle.f_curr(ferret,1) = f_next;
    if nMutantVirions == 0
        particle.f_curr(ferret,1) = 0;
    end
    if nMutantVirions == MCMC_params.Nvirions
        particle.f_curr(ferret,:) = 1;
    end

end

particle.f_list = [particle.f_list particle.f_curr];
