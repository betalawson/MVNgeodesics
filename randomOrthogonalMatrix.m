function R = randomOrthogonalMatrix(D, V)
%
%     R = randomOrthogonalMatrix( D )
%     R = randomOrthogonalMatrix( D, V )
%
% This function generates a random orthogonal matrix of the requested
% dimension, D. If a vector or matrix, V, is also supplied, then its
% columns will be used as directions in the final created matrix (if
% provided directions are not orthogonal, they will be forced orthogonal
% via Gram-Schmidt, in left-to-right order across columns).

% If no vector/matrix of directions provided, this is zero pre-specified
% directions
if nargin < 2
    P = 0;
else
    P = size(V,2);
end

% Initialise matrix
R = zeros(D);

% Create vector columns of R one by one via Gram-Schmidt
for j = 1:D
    
    % If a direction not pre-provided, create a random direction
    if j > P
        r = randn(D,1);    
    else
        if norm(V(:,j)) > 1e-15
            r = V(:,j);
        else
            r = randn(D,1);
        end
    end
    
    % Normalise to simplify following formula
    r = r / norm(r);
    
    % Subtract out any components of this vector in existing directions
    for k = 1:j-1
        r = r - (r' * R(:,k)) * R(:,k);
    end
    
    % Re-normalise and add to matrix
    R(:,j) = r / norm(r);
    
end

% Ensure the determinant of the transformation is determinant one
if det(R) < 0
    R(:,end) = -R(:,end);
end