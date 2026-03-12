function FIGURE6_covarianceEigenvalues
% This function generates a series of geodesics in bivariate MVN space,
% where the eigenvalues of the target covariance have fixed values (but
% target means and the orientation of target covariance eigenvectors are
% allowed to vary at random). All geodesics begin at the origin, henec
% representing the case where the problem has already been translated to
% the origin. The generated geodesics are used to plot the eigenvalues over
% their length, demonstrating the lower bound on eigenvalues placed by the 
% log-linear move from the origin (log lambda = 0) to the target.

% Specify the number of geodesics to generate and plot
N_geodesics = 50;

% Specify the number of plotting points along the t = [0,1] curve
N_plotpts = 251;

% Specify the bounds in mu-space in which target means can be placed
mu_min = -8;
mu_max = 8;

% Specify the fixed eigenvalues for the target covariance matrix
target_eigvals = [exp(1); exp(-1)];

% Set the random seed for reproducibility
rng(7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure the target eigenvalues are in ascending order (matches plotting
% below, which also sorts the same way)
target_eigvals = sort(target_eigvals);

% Prepare the origin in 2D MVN space
O = struct('mu',zeros(2,1),'SIGMA',eye(2));

% Split the t = [0,1] segment into the requested number of points
t_vals = linspace(0,1,N_plotpts);

% Generate the requested number of targets (done outside of parallel loop
% for easier reproducibility)
mu_stars = cell(1,N_geodesics);
SIGMA_stars = cell(1,N_geodesics);
for k = 1:N_geodesics
    
    % Generate a random target mean (uniformly) from the bounds
    mu_stars{k} = mu_min * [1;1]  + (mu_max - mu_min) * rand(2,1);
    
    % Generate a random target covariance using below sub-function
    SIGMA_stars{k} = specifiedEigvalCovariance(target_eigvals);
    
end

% Now find the geodesics connecting the origin to all targets in parallel
% for speed
Gs = cell(1,N_geodesics);
parfor k = 1:N_geodesics
    
    % Find the geodesic to the point
    Gs{k} = onePointShooting(O, struct('mu',mu_stars{k},'SIGMA',SIGMA_stars{k}) );
    
end

% Prepare figure
figure('units','Normalized','OuterPosition',[0.1 0.1 0.3 0.7]); hold on;

% Loop over each generated geodesic, and add to the plot two lines, one for
% each of its eigenvalues over the course of its length
for k = 1:N_geodesics
   
    % Grab out current geodesic
    G = Gs{k};
    
    % Calculate eigenvalues along the geodesic
    eig_vec = zeros(2,N_plotpts);
    for n = 1:N_plotpts
       
        % Fire geodesic to this point
        P = fireGeodesic(G,t_vals(n));
        
        % Calculate eigenvalues of the current covariance and sort them
        eig_vals = sort(eig(P.SIGMA));
        
        % Add these to the plot of eigenvalues along geodesic
        eig_vec(:,n) = eig_vals;
        
    end
    
    % Add the eigenvalues to the plot
    plot(t_vals, log(eig_vec(1,:)), 'LineWidth', 1.5, 'Color', [1.0, 0.4, 0.4]);
    plot(t_vals, log(eig_vec(2,:)), 'LineWidth', 1.5, 'Color', [0.4, 0.4, 1.0]);
    
end

% Plot the straight-line log-lambda connections from origin (0) to targets
plot(t_vals, log(target_eigvals(1))*t_vals, '--', 'LineWidth', 3.5, 'Color', [0.6, 0.0, 0.0]);
plot(t_vals, log(target_eigvals(2))*t_vals, '--', 'LineWidth', 3.5, 'Color', [0.0, 0.0, 0.6]);

% Clean up the plot
axis square;
xlim([0 1]);
ylim([-1.25 2.75]);
xlabel('Geodesic Parameter, t');
ylabel('log \lambda');
set(gca,'FontSize',24);
set(gca,'LineWidth',2.5);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SIGMA = specifiedEigvalCovariance(eigvals)
% This sub-function is used to generate a random covariance matrix that has
% the eigenvalues provided in the first input argument (as a vector)

% Read the dimension
d = length(eigvals);
% Generate a random orientation matrix for these eigenvalues
V = randomOrthogonalMatrix(d);
% Generate a covariance matrix with these eigenvalues
SIGMA = V * diag(eigvals) * V';