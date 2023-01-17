function u = randomVector(magnitude,D)
%
%     u = randomVector(magnitude,D)
%
% This function creates a random vector in the D-variate normal space,
% outputting a struct u with fields 'mu' and 'SIGMA'. The vector has the
% requested magnitude.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select a random amount to weight the mu vs SIGMA magnitude
r = rand;
mu_mag = sqrt(r);
SIGMA_mag = sqrt( 2 * (1 - r) );

% Pick a random direction for the mu and SIGMA components
r = randn(D,1);
mu_dir = r / norm(r);
R = randn(D);
S = 0.5 * (R + R');
SIGMA_dir = S / norm(S, 'fro');

% The direction is then the selected weighting of each direction, further
% scaled by the magnitude
u.mu = mu_mag * mu_dir * magnitude;
u.SIGMA = SIGMA_mag * SIGMA_dir * magnitude;

