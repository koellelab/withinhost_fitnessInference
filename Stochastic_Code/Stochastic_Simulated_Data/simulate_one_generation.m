function particle = simulate_one_generation(MCMC_params, particle, fitness_m, C)

% calculate discrete number of variant (mutant) virions, and discrete
% number of wild-type virions:
nMutantVirions = round(MCMC_params.Nvirions*particle.f_curr);   % equation 7
nWTVirions = MCMC_params.Nvirions - nMutantVirions;             % equation 8

% short-circuit rest of code if variant is either lost or fixed:
if nMutantVirions == 0
    particle.f_list = [particle.f_list 0];
    particle.f_curr = 0;
    return;
end
if nMutantVirions == MCMC_params.Nvirions
    particle.f_list = [particle.f_list 1];
    particle.f_curr = 1;
    return;
end

% vector of 1/C values for multinomial draw:
p = (1/C)*ones(C, 1);
% drawing from multinomial distribution:
r_mutant = mnrnd(nMutantVirions, p);
r_wt = mnrnd(nWTVirions, p);
 
% calculating mean fitness of variant virus and mean fitness of wild-type virus:
sum_f_m = 0;
sum_f_wt = 0;
for i = 1:C
    F = GetF(r_mutant(i), r_wt(i), fitness_m);
    sum_f_m = sum_f_m + r_mutant(i)*F;
    sum_f_wt = sum_f_wt + r_wt(i)*F;
end
mean_fitness_mutant = sum_f_m/nMutantVirions;   % equation 9
mean_fitness_wt = sum_f_wt/nWTVirions;          % equation 10

% equation 11:
f_next_determ = (particle.f_curr*mean_fitness_mutant)/(particle.f_curr*mean_fitness_mutant + (1-particle.f_curr)*mean_fitness_wt);
f_next_realized = binornd(MCMC_params.Nvirions, f_next_determ)/MCMC_params.Nvirions;

particle.f_list = [particle.f_list f_next_realized];
particle.f_curr = f_next_realized;
