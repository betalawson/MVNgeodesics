function u = parallelTransport(u, G, t_range)
%
%     u = parallelTransport(u, G)
%     u = parallelTransport(u, G, 'backward')
%     u = parallelTransport(u, G, t_range)
%
% This function carries out the process of parallel transporting the input
% vector 'u' along the input geodesic G. The user may optionally provide
% the range of t values along G where the transport should take place.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                 INTERNAL PARAMETERS                 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define number of steps per unit length to use for numerical integration
steps_per_length = 250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Default to transport from beginning of geodesic to end (t = 1) if no
% third input argument is provided
if nargin < 3
    t_range = [0, 1];
elseif strcmp( t_range, 'backward' )
    t_range = [1, 0];
end

% Calculate the length being integrated along
L = sqrt(innerProduct(G.v, G.v, eye(size(G.P)))) * abs(t_range(2) - t_range(1));

% Calculate the number of steps based on this and the input steps/length
N = ceil(L * steps_per_length);

% Find timestep
dt = ( t_range(2) - t_range(1) ) / N;

% Prepare the right hand side of the parallel transport ODEs
RHSfun = @(t, u) parallelTransportRHS(t, u, G);

% Initialise time
t = t_range(1);

% Loop timesteps
for k = 1:N

    % First RK4 term
    f1 = RHSfun( t, u );
    
    % Second RK4 term
    u1.mu = u.mu + dt/2*f1.mu;
    u1.SIGMA = u.SIGMA + dt/2*f1.SIGMA;
    f2 = RHSfun( t + dt/2, u1);
        
    % Third RK4 term
    u2.mu = u.mu + dt/2*f2.mu;
    u2.SIGMA = u.SIGMA + dt/2*f2.SIGMA;
    f3 = RHSfun( t + dt/2, u2 );
    
    % Fourth RK4 term
    u3.mu = u.mu + dt * f3.mu;
    u3.SIGMA = u.SIGMA + dt * f3.SIGMA;
    f4 = RHSfun( t + dt, u3 );
    
    % RK4 Update
    t = t + dt;
    u.mu = u.mu + dt/6 * ( f1.mu + 2*f2.mu + 2*f3.mu + f4.mu );
    u.SIGMA = u.SIGMA + dt/6 * ( f1.SIGMA + 2*f2.SIGMA + 2*f3.SIGMA + f4.SIGMA );
        
end


function du = parallelTransportRHS(t, u, G )
% Subfunction defining the RHS of the parallel transport ODEs for a given
% time, vector and geodesic. These equations take the form
%       du/dt = f(t,u,G)

% Fire the geodesic out to the current timepoint
p = fireGeodesic(G,t);
  
% Find the velocity of the geodesic at this timepoint (used in ODE)
dmu_dt = p.SIGMA * G.v.mu;
dSIGMA_dt = numericalProtection( p.SIGMA * ( G.v.SIGMA - G.v.mu * p.mu' ) );

% Find the inverse of the covariance
invSIGMA = p.SIGMA \ eye(size(p.SIGMA));

% Evaluate the ODE RHS
du.mu = 1/2 * (dSIGMA_dt * invSIGMA * u.mu + u.SIGMA * invSIGMA * dmu_dt);
du.SIGMA = 1/2 * (dSIGMA_dt * invSIGMA * u.SIGMA + u.SIGMA * invSIGMA * dSIGMA_dt - dmu_dt * u.mu' - u.mu * dmu_dt');