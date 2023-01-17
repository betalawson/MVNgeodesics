function [pts, D] = validatePoints(pts)
%
% This is an internal function that validates the cell array of input
% points in MVN space, pts, consists of valid MVNs and of matching
% dimension. The pts are output, with all mean vectors forced to be column
% vectors. The dimension of the inputs is also output.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If one point, process it individually
if ~iscell(pts)
    
    [pts, D] = validator(pts);
    
% Otherwise, loop through each point in the cell array, get its dimension, 
% then check that they match
else
    
    D = zeros(length(pts),1);
    for k = 1:length(pts)
        [pts{k}, D(k)] = validator(pts{k});
    end
    
    if ~all( D - D(1) == 0 )
        error('Ensure all points are of the same dimension');
    else
        D = D(1);
    end
    
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [pt, d] = validator(pt)
% This subfunction validates the input (mu, SIGMA) point as a valid point
% in the space of multivariate normals, and returns its dimension if it is 
% valid

% Validate point given as a valid struct with (mu, SIGMA) data
if ~isfield(pt,'mu') || ~isfield(pt,'SIGMA')
    error('Ensure points are specified as structs with a ''mu'' field and a ''SIGMA'' field');
else
    
    % Validate mu data is a vector
    if min(size(pt.mu)) ~= 1
        error('Ensure all points have ''mu'' field equal to a vector');
    else
        
        % Force mean vector to be a column vector
        pt.mu = pt.mu(:);
        d = length(pt.mu);
        
        % Validate dimension of covariance
        if ~isequal( size(pt.SIGMA), [d d] )
            error('Ensure all points have ''SIGMA'' field equal to a matrix of dimension equal to their ''mu'' vector');
        else
            
            % Make covariance matrix is symmetric and warn if this is a big
            % change
            if norm(pt.SIGMA - pt.SIGMA') > 1e-10
                warning('A point''s covariance matrix is assymetric, only its symmetric component will be used');
            end
            pt.SIGMA = 0.5 * (pt.SIGMA + pt.SIGMA');
            
            % Make sure the eigenvalues are all positive
            [V, LAM] = eig(pt.SIGMA);
            if any(LAM < 0)
                warning('A point''s covariance matrix is not positive definite, all negative eigenvalues will be zeroed');
                LAM(LAM < 0) = 0;
                pt.SIGMA = V * (LAM / V);
            end
            
        end 
    end 
end               
            