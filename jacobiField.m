function J = jacobiField(G, dv, t)
%
%     J = jacobiField(G, dv, t)
%
% This function numerically approximates the Jacobi field, J(t), associated
% with the input geodesic G and an update to its velocity/tangent dv.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the norm of the update vector to select an update length
h = 1e-6 / sqrt( innerProduct(dv, dv, eye(size(dv.SIGMA)) ) );

% Create geodesic with small shift from G in direction dv
Gh = G;
Gh.v = struct('mu', G.v.mu + h*dv.mu, 'SIGMA', G.v.SIGMA + h*dv.SIGMA);

% Fire the two geodesics
p = fireGeodesic(G,t);
ph = fireGeodesic(Gh,t);

% Finite difference the two endpoints
J.mu = (ph.mu - p.mu) / h;
J.SIGMA = (ph.SIGMA - p.SIGMA) / h;
