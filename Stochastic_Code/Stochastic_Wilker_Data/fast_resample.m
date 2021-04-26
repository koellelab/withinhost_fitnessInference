function k_vector = fast_resample(w, MCMC_params)

w_cdf = cumsum(w)/sum(w);
if max(w_cdf) ~= 1
    k_vector = 1:length(w); % don't resample
    return;
end

k_vector(1, 1:MCMC_params.n_particles) = NaN;
for i = 1:MCMC_params.n_particles
    locs = find(w_cdf > rand);  
    k_vector(i) = locs(1);
end
