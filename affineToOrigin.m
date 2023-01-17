function [O,pt,P,r] = affineToOrigin(p1,p2)
% 
%     [O,pt,P,r] = affineToOrigin(p1,p2)
%
% This function as inputs two multivariate normals, p1 and p2, specified
% as structs with a 'mu' and 'SIGMA' field. An affine transformation is
% selected such that p1 becomes the origin, O, and p2 becomes a new
% transformed point, pt. The parameters (P,r) of the transform FROM the
% origin TO the target point are also returned.
%
% Note: There is an infinite number of ways to select the transform P -
% any arbitrary rotation of a decomposition P P^T = SIGMA is valid. This
% function uses the matrix square root, so as to return a symmetric result.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find a decomposition P P^T = SIGMA1
P = p1.SIGMA^(1/2);

% The shift in mean is simply the shift that gets back to mu1
r = p1.mu;

% The second point shifts according to the transform in reverse
invP = p1.SIGMA^(-1/2);    % Assuming this is better than a backslash for the inverse
pt.mu = invP * (p2.mu - p1.mu);
pt.SIGMA = invP * p2.SIGMA * invP';

% We know we transform to the origin, so just set it manually
O.mu = zeros(length(p1.mu),1);
O.SIGMA = eye(length(p1.mu));

end

