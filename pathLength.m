function L = pathLength(path, useSYMKL)
% This function takes a path through MVN space as input, and returns the
% information geometric distance along the MVN statistical manifold. This
% distance is estimated by adding the distance associated with each
% increment in the path - as such, a path consisting of a large number of
% points should be provided for reasonable accuracy

% Read out the number of points in the path
N = length(path);

% Assume use of metric-based approximations if not specified
if nargin < 2
    useSYMKL = false;
end

% Loop over each element in the provided path, and add its contribution to
% the length using the requested approximation
L = 0;
for k = 1:N-1
    
    % Symmetrised KL as approximation - closed form for MVNs
    if useSYMKL
        
        % Read out the current and next path points
        mu1 = path{k}.mu;
        SIGMA1 = path{k}.SIGMA;
        mu2 = path{k+1}.mu;
        SIGMA2 = path{k+1}.SIGMA;
        % Use these to calculate
        L = L + sqrt(0.5 * ( trace(SIGMA2 \ SIGMA1) + trace(SIGMA1 \ SIGMA2) + (mu2-mu1)' * ( SIGMA2 \ (mu2-mu1) + SIGMA1 \ (mu2-mu1) ) - 2*length(mu1) ) );
        
    % Metric-based approximation - use differential form but now we have
    % \Delta \theta  instead of d\theta
    else
       
        % Store the tangent vector object
        dtheta.mu = path{k+1}.mu - path{k}.mu;
        dtheta.SIGMA = path{k+1}.SIGMA - path{k}.SIGMA;
        % Find inverse and evaluate inner product at start point
        invSIGMA = path{k}.SIGMA \ eye(size(path{k}.SIGMA));
        L = L + sqrt( innerProduct(dtheta, dtheta, invSIGMA) );
        
    end
    
end