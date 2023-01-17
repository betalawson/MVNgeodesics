function options = OPTIONSdefaults()
% This function defines the default options used by the shooting methods
% for finding geodesic connections between multivariate normals. These
% options can be manually overriden in calls to those functions, or the
% default behaviours may be adjusted here.

%%% GENERAL PARAMETERS

options.symKL_tol = 1e-8;          % Tolerance to reach before terminating
options.visualise = false;         % Specify whether to visualise or not
options.make_animation = false;    % Specify whether to save an animation (only used if 'visualise' is true)
options.verbose = false;           % Specify whether to provide text output to user

%%% INITIAL GEODESIC

options.Ginit = [];                % This should not be changed. However, you may provide a geodesic with which to initialise the onePointShooting method by setting this value in the input options

%%% SINGLE SHOOTING PARAMETERS

% Method used for approximate geodesics that inform velocity corrections.
% Approaches included are:
%   'straight' - Velocity assuming Euclidean space (also 'euclid')
%   'eigen' - Eigendecomposition-based velocity
%   'taylor' - Taylor series approach assuming small velocity in mu space
%   'transport' - Parallel transport based approach
options.approx_method = 'euclid';

% Maximum norm of the velocity vector before its contribution to the
% velocity update is downweighted
options.max_vnorm = 0.5;

% Looping options
options.max_iters = 5000;          % Maximum number of iterations to try



%%% MULTI SHOOTING PARAMETERS

% Looping options
options.max_multi_iters = 500;     % Maximum number of iterations to try

% Number of points to use (note: will be forced odd for multi point and to
% next value of 2^m + 1 for adaptive multi point)
options.N_points = 128;

% Path between p1 and p2 used for point initialisation. Options are:
%   'euclid' - Direct interpolation of (mu, SIGMA) between p1 and p2
%   'geometric' - Direct interpolation of (delta, DELTA) between p1 and p2
%   'moment' - Interpolation of the first and second moments of MVNs
%   'hybrid' - Average of the geometric and moment paths in (mu,SIGMA)
options.initial_path = 'euclid';

%%% ADAPTIVE MULTI SHOOTING PARAMETERS

% Required relative improvement (reduction) in length per path refinement
% in the multi-point iteration 
options.req_L_improvement = 1e-2;   % (Require at least a 1% reduction)

% Minimum number of iterations to spend at each number of points
options.min_Npts_iters = 5;

end

