function v = randomUnitVector(d, SIGMA)
% This function generates a random unit vector from the tangent space for
% MVNs. The user may optionally provide the inverse covariance (base
% point), otherwise it will be assumed to be the origin

% Assume the origin if no base point provided, otherwise, calculate inverse
if nargin < 2
    SIGMA = eye(d);
    invSIGMA = eye(d);
else
    invSIGMA = SIGMA \ eye(d);
end

% Generate d values distributed according to the base point covariance for
% the mu vector (these will be 'sphere' in the sense of the metric, which
% has term dmu^T Sigma^-1 dmu
vmu = mvnrnd(zeros(d,1),SIGMA)';

% Prepare a set of lower-triangular elements
tri_ind = [];
loc = 0;
for k = 1:d
    tri_ind(loc+1:loc+1+d-k) = d*(k-1)+(k-1)+(1:1+d-k);
    loc = loc + 1+d-k;
end

% Generate (d+1)d/2 normally-distributed values for SIGMA-component
vSIGMA = zeros(d);
vSIGMA(tri_ind) = randn(d*(d+1)/2,1);
% Copy the off-diagonal elements across to achieve symmetry
vSIGMA = vSIGMA + vSIGMA' - diag(diag(vSIGMA));
% Double the variance along the diagonal to counteract the value of 1/2 in 
% the metric (off-diagonals contribute twice,so they counteract it alraeady)
vSIGMA(logical(eye(d))) = vSIGMA(logical(eye(d))) * sqrt(2);
% Affine transform this so that it's sphere in the sense of the metric,
% consider || SIGMA0^(-1/2) dSIGMA SIGMA0^(-1/2) ||_F^2
vSIGMA = SIGMA^(1/2) * vSIGMA * SIGMA^(1/2);

% Calculate the length of the resultant vector
v = struct('mu',vmu,'SIGMA',vSIGMA);
L = sqrt(innerProduct(v,v,invSIGMA));
% Normalise
v.mu = v.mu / L;
v.SIGMA = v.SIGMA / L;

end