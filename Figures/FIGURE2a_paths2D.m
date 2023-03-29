function FIGURE2a_paths2D()
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
T_SIGMAS = {[5, 0; 0, 0.5], [0.5, 0; 0, 5], [0.4, -0.375; -0.375, 0.4] };

% Specify the text offset
text_offset = [0.1, -0.1];

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
    Peuclid = closedFormPath(O, T, 5001, 'euclid');
    Panneal = closedFormPath(O, T, 5001, 'anneal');
    Pmoment = closedFormPath(O, T, 5001, 'moment');
    Pwasserstein = closedFormPath(O, T, 5001, 'transport');
    
    % Plot the three paths
    plotPath(Panneal, 'pathColor', [1, 0.4, 0.4], 'ellipseColor', [1 0.4 0.4], 'ellipseWidth', 2, 'pathWidth', 2.25, 'covFrequency', 1000);
    plotPath(Pmoment, 'pathColor', [0.4, 0.4, 1], 'ellipseColor', [0.4 0.4 1], 'ellipseWidth', 2, 'pathWidth', 2.25, 'covFrequency', 1000);
    plotPath(Pwasserstein, 'pathStyle','--', 'pathColor', [0.9, 0.5, 0.0], 'ellipseColor', [0.9 0.5 0.0], 'ellipseWidth', 2, 'pathWidth', 2.25, 'covFrequency', 1000);
    plotGeodesic(G, 'pathColor', [0, 0, 0], 'ellipseColor', [0 0 0], 'ellipseWidth', 3, 'pathWidth', 4.5, 'covFrequency', 1000, 'Npts', 5001);

    % Calculate the lengths of each path
    Lgeo = sqrt( innerProduct(G.v, G.v, eye(2) ) );
    Lanneal = pathLength(Panneal);
    Lmoment = pathLength(Pmoment);
    Lwasserstein = pathLength(Pwasserstein);
    
    % Trace out the geodesic halfway to get the point where to put the
    % label
    Pgeo_mid = fireGeodesic(G,0.5);
    
    % Append the lengths of each path at their middle point
    text( Panneal{2501}.mu(1) + text_offset(1), Panneal{2501}.mu(2) + text_offset(2), ['L = ',num2str(Lanneal)], 'Color', [1 0.4 0.4], 'FontSize', 24);
    text( Pmoment{2501}.mu(1) + text_offset(1), Pmoment{2501}.mu(2) + text_offset(2), ['L = ',num2str(Lmoment)], 'Color', [0.4 0.4 1], 'Fontsize', 24);
    text( Pwasserstein{2501}.mu(1) + text_offset(1), Pwasserstein{2501}.mu(2) + text_offset(2), ['L = ',num2str(Lwasserstein)], 'Color', [0.9 0.5 0], 'Fontsize', 24);
    text( Pgeo_mid.mu(1) + text_offset(1), Pgeo_mid.mu(2) + text_offset(2), ['L = ',num2str(Lgeo)], 'Color', [0 0 0], 'FontSize', 24);        
    
end

% Return to the Figures folder
cd('Figures');

% Append labels and legend
xlabel('\mu_1','FontSize',24);
ylabel('\mu_2','FontSize',24);
legend({'Annealing Path', '', '', '', '', '', '', '', '', 'Moment-Averaged Path', '', '', '', '', '', '', '', '', 'Wasserstein Path', '', '', '', '', '', '', '', '', 'Geodesic Path'}, 'FontSize', 24, 'Location','NorthWest');
set(gca, 'FontSize', 24);

% Fix axis limits and thickness
xlim([-2.25 2.75]);
ylim([-1 4]);
axis square;
set(gca,'LineWidth',2);

