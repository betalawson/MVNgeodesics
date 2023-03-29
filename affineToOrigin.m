function [O,varargout] = affineToOrigin(p1,p2)
% 
%     [O,P,r] = affineToOrigin(p1)
%     [O,pt,P,r] = affineToOrigin(p1,p2)
%
% This function takes as input at least one multivariate normal, p1,
% specified as a struct with a 'mu' field and a 'SIGMA' field. An affine
% transformation is selected such that p1 becomes the origin, O, and the
% parameters (P,r) of the affine transform FROM the origin BACK TO the
% target point are also returned.
%
% If a second point, p2, is provided, the function also returns (as its
% second output) the second point's new location after application of the
% same affine transform.
%
% Note: There is an infinite number of ways to select the transform P -
% any arbitrary rotation of a decomposition P P^T = SIGMA is valid. This
% function uses the matrix square root, so as to return a symmetric result.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find a decomposition P P^T = SIGMA1
P = p1.SIGMA^(1/2);

% The shift in mean is simply the shift that gets back to mu1
r = p1.mu;

% We know we transform to the origin, so just set it manually
O.mu = zeros(length(p1.mu),1);
O.SIGMA = eye(length(p1.mu));

% If a second point was provided, transform it also
if nargin > 1
    
    % Find the inverse transform to be applied to p2
    invP = p1.SIGMA^(-1/2);
    
    % Apply the inverse transform to p2
    pt.mu = invP * (p2.mu - p1.mu);
    pt.SIGMA = invP * p2.SIGMA * invP';
    
    % Gather outputs
    varargout{1} = pt;
    varargout{2} = P;
    varargout{3} = r;
    
% Otherwise, just gather the outputs
else
   
    varargout{1} = P;
    varargout{2} = r;
    
end

end