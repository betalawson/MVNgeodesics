function FIGURE3_univariatePaths()
% This function visualises the different paths between an example pair of 
% univariate normals that are generated using popular path-generating 
% approaches in the literature, or that are optimal according to different
% distance metrics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SPECIFICATIONS

% Specify the two distributions to connect
P = struct('mu', 0, 'SIGMA', 0.25);
Q = struct('mu', 2.5, 'SIGMA', 0.05);

% Other parameters
vpts = 201;             % Number of evaluation points for plotting
Npts = 11;              % Number of path points to show on plot
Npts_calc = 50001;      % Number of path points to use in length calculation
epsilon1 = 1;           % Extent of entropy regularisation (first subfig)
epsilon2 = 5;           % Extent of entropy regularisation (second subfig)

% Plotting color
plot_color = [0.5, 0.5, 1.0];
           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare the figure
figure('units','Normalized','OuterPosition',[0.26, 0.11, 0.44, 0.73]);
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% Prepare a vector of evaluation points
minval = min([P.mu - 3*sqrt(P.SIGMA); Q.mu - 3*sqrt(Q.SIGMA)]);
maxval = max([P.mu + 3*sqrt(P.SIGMA); Q.mu + 3*sqrt(Q.SIGMA)]);
% Use the more distant limit as the actual limit-setter for symmetry
mu_mean = 0.5 * (P.mu + Q.mu);
L = 2.5*max([maxval - mu_mean, mu_mean - minval]);
minval = mu_mean - L/2;
maxval = mu_mean + L/2;
% Set up the vector of evaluation points (for pdf plotting)
xv = linspace(minval, maxval, vpts);

% Generate the paths
Wpath = closedFormPath(P, Q, Npts, 'W2');
entWpath1 = closedFormPath(P, Q, Npts, 'EW2', epsilon1);
entWpath2 = closedFormPath(P, Q, Npts, 'EW2', epsilon2);
Apath = closedFormPath(P, Q, Npts, 'annealing');
Mpath = closedFormPath(P, Q, Npts, 'moment');

% Generate Fisher-Rao path using shooting for ease
G = onePointShooting(P,Q);
FRpath = traceGeodesic(G,[0 1],Npts);

% Generate high-resolution versions of the paths
Wpath_calc = closedFormPath(P, Q, Npts_calc, 'W2');
entWpath1_calc = closedFormPath(P, Q, Npts_calc, 'EW2', epsilon1);
entWpath2_calc = closedFormPath(P, Q, Npts_calc, 'EW2', epsilon2);
Apath_calc = closedFormPath(P, Q, Npts_calc, 'annealing');
Mpath_calc = closedFormPath(P, Q, Npts_calc, 'moment');

% Use these to approximate the distance under Fisher metric along them
dW = pathLength(Wpath_calc);
deW1 = pathLength(entWpath1_calc);
deW2 = pathLength(entWpath2_calc);
dA = pathLength(Apath_calc);
dM = pathLength(Mpath_calc);
% Length of Fisher geodesic is stored in the geodesic object already
dFR = G.L;

% Put all the paths in one array for plotting ease
plot_paths = {Wpath, Apath, Mpath, entWpath1, entWpath2, FRpath};
path_lengths = [dW, dA, dM, deW1, deW2, dFR];
path_names = {'$d_{W}$-geodesic','Annealing', 'Moment-averaged', ['$d_{\epsilon W}$-geodesic ($\epsilon = ', num2str(epsilon1), '$)'], ['$d_{\epsilon W}$-geodesic ($\epsilon = ', num2str(epsilon2), '$)'], '$d_{F\!R}$-geodesic'};

% Loop over each path in the list, adding each to the plot
for p = 1:length(plot_paths)

    % Prepare plot location
    nexttile; hold on;
    
    % Plot the points along this path
    for k = 2:Npts-1
        pkv = normpdf(xv, plot_paths{p}{k}.mu, sqrt(plot_paths{p}{k}.SIGMA));
        plot(xv, pkv, 'LineWidth', 2, 'Color', plot_color);
    end
    
    % Plot the base Gaussians last so they are on top
    pv = normpdf(xv,P.mu,sqrt(P.SIGMA));
    plot( xv, pv, 'k', 'LineWidth', 3);
    qv = normpdf(xv,Q.mu,sqrt(Q.SIGMA));
    plot( xv, qv, 'k', 'LineWidth', 3);
    
    % Clean up plot
    axis square;
    axis off;
    title(path_names{p},'FontSize', 20, 'Interpreter', 'LaTeX');
    
    % Add the distance as an xlabel
    text(mean(xv), -0.20, ['$L = ',num2str(path_lengths(p),'%.4f'),'$'],'FontSize',22, 'HorizontalAlignment', 'Center', 'Interpreter', 'LaTeX');
    
end  



end
