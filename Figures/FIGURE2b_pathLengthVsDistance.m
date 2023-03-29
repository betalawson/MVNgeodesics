function FIGURE2b_pathLengthVsDistance(regenerate)
% This function generates a large number of end points for geodesics
% (beginning at the origin), and produces a box plot comparing the lengths
% of various non-geodesic paths

% Specify number of samples
N_samples = 10000;

% Define the dimension
D = 2;

% Set the bounds for mean space
mu_bound = 10;
% Set the bounds for covariance eigenvalues (multiplicative factor)
SIGMA_bound = 10;

% Define plot colours
plot_colors = [
    [1.0, 0.4, 0.4];
    [0.4, 0.4, 1.0];
    [0.9, 0.5, 0.0];
    ];

% Define path names (for legend)
path_names = {'Euclidean', 'Annealing', 'Moment-Averaging', 'Wasserstein'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume, if no flag provided, not to regenerate data
if nargin < 1
    regenerate = false;
end

% Regenerate data if requested/needed
if ~exist('data_pathPerformance.mat','file') || regenerate
    
    % Move back into the main folder
    cd('..');
    
    
    % Define the origin in this dimension
    O.mu = zeros(D,1);
    O.SIGMA = eye(D);
    
    % Loop over requested number of samples
    parfor k = 1:N_samples
        
        % Create a random endpoint
        pt = randomPoint(mu_bound, SIGMA_bound, D);
        
        % Find the length of the geodesic to this point
        G = onePointShooting(O,pt);
        Lgeo = G.L;
        
        % Now find the non-geodesic paths to this location
        Peuclid = closedFormPath(O, pt, 5001, 'euclid');
        Panneal = closedFormPath(O, pt, 5001, 'anneal');
        Pmoment = closedFormPath(O, pt, 5001, 'moment');
        Pwasserstein = closedFormPath(O, pt, 5001, 'transport');
        
        % Find the lengths of these paths
        Leuclid = pathLength(Peuclid);
        Lanneal = pathLength(Panneal);
        Lmoment = pathLength(Pmoment);
        Lwasserstein = pathLength(Pwasserstein);
        
        % Add to the data matrix the log % error
        rel_err(k,:) = [Leuclid - Lgeo, Lanneal - Lgeo, Lmoment - Lgeo, Lwasserstein - Lgeo] / Lgeo;
        
    end
    
    % Remove any values that show up as negative relative error, as these
    % correspond to the numerical length calculation incorrectly asserting a
    % better-than-geodesic length
    rel_err( rel_err < 0 ) = NaN;
    
    
    % Move back to original folder
    cd('Figures');
    
    % Save the data
    save('data_pathPerformance.mat','rel_err');
    
else
    
    % Load in the data
    load('data_pathPerformance.mat','rel_err');
    
end

% Specify which columns of the data to plot (and order of plotting)
use_columns = [2, 3, 4];
% use_columns = 1:size(rel_err,2)

% Change the shape of the error object into one compatible with the boxplot
% function used, also take the logarithm because it's a ratio
c = 0;
for k = use_columns
    c = c + 1;
    log_rel_err(1,c,:) = log10( rel_err(:,k) );
end

% Initialise figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);

% Create base boxplot
h = boxplot2( log_rel_err, [1] );

% Color the different path types
for k = 1:size(log_rel_err,2)
    structfun(@(x) set( x(k,:), 'color', plot_colors(k,:), 'markeredgecolor', plot_colors(k,:), 'linewidth', 2.5), h);
end

% Grab out locations of each bit of data to set Xticks
for k = 1:size(log_rel_err,2)
    x_locs(k) = h.lwhis(k).XData(1);
end

% Make these the XTicks and label them
set(gca,'XTick',x_locs);
set(gca,'XTickLabel',path_names(use_columns));

% Set x-axis limits
x_diff = x_locs(2) - x_locs(1);
xlim([x_locs(1) - 0.5*x_diff, x_locs(end) + 0.5*x_diff]);

% Other plot cleanup
ylabel('log_1_0 (L - d_F)/d_F','FontSize',24);
set(gca,'Fontsize', 24);
legend('off');

% Print out the median error of each error category
for k = 1:size(log_rel_err,2)
   
    % Read out the median line for this class
    err = 10^( h.med(k).YData(1) );
    
    % Print out the error
    fprintf('Median error for %s is %g%%\n', path_names{use_columns(k)}, err*100);
    
end