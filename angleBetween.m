function value = angleBetween(u, v, invSIGMA)
%
%     value = angleBetween(u, v, invSIGMA)
%     value = angleBetween(u, v)
%
% This function calculates the angle (in radians) between two vectors u and
% v in the space of multivariate normal distributions, with each provided 
% as a struct defining its mu and SIGMA components.
%
% The angle is calculated using the normalised inner product, with the
% inner product being the one that induces the Riemannian metric on this
% statistical manifold. Given the warping of space, the calculation depends
% on the precision matrix where the angle is being calculated. If no input
% precision invSIGMA is provided, the origin invSIGMA = eye(D) is used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume calculation takes place at the origin if no point provided
if nargin < 3
    invSIGMA = eye(size(v.SIGMA));
end

% Find the angle using cos theta = <u,v> / sqrt( <u,u> <v,v> )
value = acos( innerProduct(u,v,invSIGMA) / sqrt( innerProduct(u,u,invSIGMA) * innerProduct(v,v,invSIGMA) ) );

end

