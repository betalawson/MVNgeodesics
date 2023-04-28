function FIGURE4a_paths2D()
% This function shows a couple of example geodesics in the space of
% bivariate normals, as well as visualising the corresponding paths formed
% using a Euclidean connection and an annealing connection

% Initialise the figure
figure('units','normalized','OuterPosition',[0 0 1 1]);
hold on;

% Specify the origin for bivariate MVNs
O.mu = [0;0];
O.SIGMA = eye(2);

% Specify a few points to find geodesics to
T_mus = [ [2;3.5], [2; 0], [-2; 0] ];
T_SIGMAS = {[5, 0; 0, 0.5], [0.5, 0; 0, 5], [0.8, -0.6; -0.6, 0.8] };

% Specify the text offset
dy = 0.225;

% Specify the text locations
text_x = [2.35, 2.2, -2.1];
text_y = [3.3, 0.3, 1];

% Specify whether to plot all curves, or just example geodesics
show_all_curves = false;

% Return to the original folder
cd('..');

% Loop over each, plotting the geodesic path along with the other paths
for k = 1:size(T_mus,2)
    
    % Grab out the mu and SIGMA here
    T.mu = T_mus(:,k);
    T.SIGMA = T_SIGMAS{k};
    
    % Find the geodesic path
    G = onePointShooting(O, T);
        
    % Find the three paths
    if show_all_curves
        Peuclid = closedFormPath(O, T, 4001, 'euclid');
        Panneal = closedFormPath(O, T, 4001, 'anneal');
        Pmoment = closedFormPath(O, T, 4001, 'moment');
        Pwasserstein = closedFormPath(O, T, 4001, 'transport');
    
        % Find the lengths of these paths (numerically)
        Lanneal = pathLength(Panneal);
        Lmoment = pathLength(Pmoment);
        Lwasserstein = pathLength(Pwasserstein);
    end
       
    % Plot the geodesic
    plotGeodesic(G, 'pathColor', [0, 0, 0], 'ellipseColor', [0 0 0], 'ellipseWidth', 3, 'pathWidth', 4.5, 'covFrequency', 1000, 'Npts', 4001);

    % Plot the three paths if requested
    if show_all_curves
        
        % Curve plotting
        plotPath(Panneal, 'pathColor', [1, 0.4, 0.4], 'ellipseColor', [1 0.4 0.4], 'ellipseWidth', 3, 'pathWidth', 3, 'covFrequency', 1000);
        plotPath(Pmoment, 'pathColor', [0.4, 0.4, 1], 'ellipseColor', [0.4 0.4 1], 'ellipseWidth', 3, 'pathWidth', 3, 'covFrequency', 1000);
        plotPath(Pwasserstein, 'pathStyle','--', 'pathColor', [0.9, 0.5, 0.0], 'ellipseColor', [0.9 0.5 0.0], 'ellipseWidth', 3, 'pathWidth', 3, 'covFrequency', 1000);
        
        % Append the lengths of each path to the plot
        text( text_x(k), text_y(k) - 0*dy, ['L = ',num2str(Lanneal,4)], 'Color', [1 0.4 0.4], 'FontSize', 24);
        text( text_x(k), text_y(k) - 1*dy, ['L = ',num2str(Lmoment,4)], 'Color', [0.4 0.4 1], 'FontSize', 24);
        text( text_x(k), text_y(k) - 2*dy, ['L = ',num2str(Lwasserstein,4)], 'Color', [0.9 0.5 0], 'FontSize', 24);
        
    end
        
    % Append just the geodesic's length
    text( text_x(k), text_y(k) - 3*dy, ['L = ',num2str(G.L,4)], 'Color', [0 0 0], 'FontSize', 24);
        
end

% Return to the Figures folder
cd('Figures');

% Append labels and legend if plotting multiple curves
xlabel('\mu_1','FontSize',24);
ylabel('\mu_2','FontSize',24);
if show_all_curves
    legend({'Annealing Path', '', '', '', '', '', '', '', '', 'Moment-Averaged Path', '', '', '', '', '', '', '', '', 'Wasserstein Path', '', '', '', '', '', '', '', '', 'Geodesic Path'}, 'FontSize', 24, 'Location','NorthWest');
end
set(gca, 'FontSize', 24);

% Fix axis limits and thickness
xlim([-2.25 2.75]);
ylim([-1 4]);
axis square;
set(gca,'LineWidth',2);

