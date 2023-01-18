function A = numericalProtection(A)
%
% This internal function simply ensures that the matrix A is symmetric, and
% any tiny elements that are likely rounding errors are deleted

% Specify tolerance
tol = 1e-12;

% Symmetrise
A = 0.5 * (A + A');
% Remove any elements very close to the identity
dA = A - eye(size(A));
dA( abs(dA) < tol ) = 0;
A = eye(size(A)) + dA;