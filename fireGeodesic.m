function p = fireGeodesic(G,t)
%
%     p = fireGeodesic(G)
%     p = fireGeodesic(G, t)
%
% This function fires the input geodesic G out to time t. The geodesic
% should be provided as a struct with components P, r and v. The velocity
% v is the velocity in "origin space", and should be specified as a struct
% itself with a 'mu' field and a 'SIGMA' field.
%
% If no input 't' is provided, the geodesic is fired out to t=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    t = 1;
end

% Grab out the initial velocity components
x = G.v.mu;
B = G.v.SIGMA;

% Read out the dimension
d = length(x);

% Formulate the matrix that may be exponentiated to find geodesic paths
A = [ [    -B    ,  x  , zeros(d)  ];
      [    x'    ,  0  ,    -x'    ];
      [ zeros(d) , -x  ,     B     ] ];

% Exponentiate the matrix
LAMBDA = expm( A * t );

% Read out the geodesic (in temrs of canonical variables)
DELTA = LAMBDA(1:d, 1:d);
delta = LAMBDA(1:d, d+1);
% Convert back to true components
SIGMA = numericalProtection( DELTA \ eye(d) );

% Store the point - still in origin space
p.mu = SIGMA * delta;
p.SIGMA = SIGMA;

% Transform back to the original space as defined by the geodesic object
p = affineTransform( p, G.P, G.r );
    
end
