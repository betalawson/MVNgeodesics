function sKL = symKL( p1, p2 )
%
%     sKL = symKL( p1, p2 )
%
% This function calculates the symmeterised KL divergence,
%
%     sKL = 1/2 KL( p1, p2 ) + 1/2 KL( p2, p1 ),
%
% between multivariate normal distributions, expressed as two structs p1 
% and p2 that each has a 'mu' and 'SIGMA' component

% Find the inverse of the covariance matrices
d = length(p1.mu);
invSIGMA1 = p1.SIGMA \ eye(d);
invSIGMA2 = p2.SIGMA \ eye(d);

% Calculate the symmeterised KL
sKL = 0.25 * ( trace( invSIGMA1 * p2.SIGMA ) + trace( invSIGMA2 * p1.SIGMA ) + ( p1.mu - p2.mu )' * (invSIGMA1 + invSIGMA2) * ( p1.mu - p2.mu ) - 2*d);

end