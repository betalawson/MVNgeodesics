function L = pathLength(path)
% This function takes a path through MVN space as input, and returns the
% information geometric distance along the MVN statistical manifold. This
% distance is estimated by adding the distance associated with each
% increment in the path - as such, a path consisting of a large number of
% points should be provided for reasonable accuracy

% Read out the number of points in the path
N = length(path);

% Loop over each element in the provided path, and add its contribution to
% the length
L = 0;
for k = 1:N-1
    
    % Read out the increments in the variables 
    dtheta.mu = path{k+1}.mu - path{k}.mu;
    dtheta.SIGMA = path{k+1}.SIGMA - path{k}.SIGMA;
        
    % Add on the corresponding length increment
    invSIGMA = inv(path{k}.SIGMA);
    L = L + sqrt( innerProduct(dtheta, dtheta, invSIGMA) );
    
end