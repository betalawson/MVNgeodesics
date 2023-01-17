function value = innerProduct(u,v,invSIGMA)
%
%     value = innerProduct(u, v, invSIGMA)
%     value = innerProduct(u, v)
%
% This function calculates the inner product between two vectors u and
% v in the space of multivariate normal distributions, with each provided 
% as a struct defining its mu and SIGMA components. 
%
% The inner product is the one that induces the Riemannian metric on this 
% statistical manifold. Given the warping of space, the calculation depends
% on the precision  matrix where the angle is being calculated. If no input
% precision invSIGMA is provided, the origin invSIGMA = eye(D) is used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate inner product by sticking u, v and the precision into the metric
value = u.mu' * invSIGMA * v.mu + 1/2 * trace( invSIGMA * u.SIGMA * invSIGMA * v.SIGMA );

end

