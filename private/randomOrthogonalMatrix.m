function V = randomOrthogonalMatrix(input)
% This internal function generates a random orthogonal matrix of the 
% requested dimension. "input" may either be the dimension desired (if a 
% scalar), or a matrix of directions that should be directions in the
% resulting matrix (they will be forced orthogonal if not, in order across
% columns).

% Read out dimension and number of pre-provided vectors
if length(input) == 1
    D = input;
    P = 0;
else
    [D,P] = size(input);
end

% Initialise matrix
V = zeros(D);

for j = 1:D
    
    % If a direction not pre-provided, create a random direction
    if j > P
        r = randn(D,1);    
    else
        if norm(input(:,j)) > 1e-15
            r = input(:,j);
        else
            r = randn(D,1);
        end
    end
    
    % Normalise to simplify following formula
    r = r / norm(r);
    
    % Subtract out any components of this vector in existing directions
    for k = 1:j-1
        r = r - (r' * V(:,k)) * V(:,k);
    end
    
    % Re-normalise and add to matrix
    V(:,j) = r / norm(r);
    
end

% Ensure the determinant of the transformation is determinant one
if det(V) < 0
    V(:,end) = -V(:,end);
end