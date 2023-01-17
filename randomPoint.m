function p = randomPoint(mu_bound, SIGMA_bound, D)
%
%     p = randomPoint(mu_bound, SIGMA_bound, D)
%
% This function returns a random multivariate normal as a struct with 'mu'
% and 'SIGMA' fields. The mu value is uniformly distributed in the
% hypercube [-mu_bound, mu_bound]^D, and the covariance matrix has
% eigenvalues distributed lognormally between [1/SIGMA_bound, SIGMA_bound]
% and with a random orientation

% Random mean, simply uniformly distributed in hypercube
p.mu = -mu_bound + rand(D,1) * 2 * mu_bound;

% Generate random eigenvalues, loguniformly distributed
lambda = exp( -log(SIGMA_bound) + rand(D,1) * 2 * log(SIGMA_bound) ); 

% Generate a random orthogonal matrix to decide the orientation
R = randomOrthogonalMatrix(D);

% Use this together with eigenvalues to define the covariance matrix
p.SIGMA = R * diag(lambda) * R';

end

