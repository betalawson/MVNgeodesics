function FIGURE5a_exampleGeodesics
%
% This function plots a series of geodesics in the space of bivariate MVNs,
% together with the geodesics obtained by approximately selecting the
% initial velocity. All geodesics begin at the origin, and the user
% provides a list of endpoints.

% List the points to target
targets{1} = struct('mu', [-0.4;0.75],'SIGMA', [1, 0.6; 0.6, 0.6]);
targets{2} = struct('mu', [1;-0.5],'SIGMA', [2.5,0;0,0.5]);
targets{3} = struct('mu', [-1, -0.5],'SIGMA', [1.2, 0; 0, 0.8]);
%targets{4} = struct('mu', ,'SIGMA', );

% Create another point that deliberately has a straight line solution
targets{4} = struct('mu', [1; 1], 'SIGMA', [1.45, -0.75; -0.75, 1.45]);


% Define the colours
plot_colors = [ 
                0.40, 0.40, 0.40;        % Grey - Euclidean
                1.00, 0.40, 0.40;        % Red - Taylor
                0.40, 0.40, 1.00;        % Blue - Eigen
              ];  
          
% List method names
methods = {'euclid','taylor','eigen'};

% Corresponding method names for display on legend
legend_names = {'Euclidean', 'Taylor', 'Eigen'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise figure
figure('units', 'normalized', 'OuterPosition', [0 0 1 1]);
hold on;

% Read out dimension
D = length(targets{1}.mu);

% Create origin in this dimension
O = struct('mu',zeros(D,1),'SIGMA',eye(D));

% Jump back to main folder
cd('..');

% Append geodesic to legend names
legend_names{length(legend_names)+1} = 'Correct Geodesic';

% Loop over the provided target points
for k = 1:length(targets)
    
    % Loop over the methods
    for m = 1:length(methods)
        
        % Formulate the approximate geodesics
        Gapprox = approximateGeodesic(O, targets{k}, methods{m});
    
        % Plot this geodesic on the figure
        plotGeodesic(Gapprox, 'pointSize', 1, 'ellipseColor', plot_colors(m,:), 'pathColor', plot_colors(m,:), 'endColor', plot_colors(m,:), 'covFrequency', 50, 'Npts', 151, 'ellipseWidth',1.5, 'pathWidth', 2.5);
        
        % Find the geodesic's end point separately to store it
        endpt = fireGeodesic(Gapprox,1);
        endpt_mu(m,1:D) = endpt.mu;
        
    end    

    % Find the true geodesic to this point
    G = onePointShooting(O, targets{k});
    
    % Plot the true geodesic in black
    plotGeodesic(G, 'pathStyle', ':', 'pathColor',[0 0 0], 'pathWidth', 4, 'pointSize', 40, 'ellipseColor', [0 0 0], 'endColor', [0 0 0], 'covFrequency', 50, 'Npts', 151, 'ellipseWidth', 1.5);
    
    % Plot the approximate geodesic end points again on top
    for m = 1:length(methods)
        plot(endpt_mu(m,1), endpt_mu(m,2), '.', 'MarkerSize', 25, 'MarkerFaceColor', plot_colors(m,:), 'MarkerEdgeColor', plot_colors(m,:));
    end
    
    % Also plot the origin on top to make sure it's black
    plot(0, 0, '.', 'MarkerSize', 40, 'MarkerEdgeColor', [0 0 0]);
    
end

% Set up legend text (counts how many lines plotted to avoid labelling
% them)
N_ellipse = floor(151/50);
c = 0;
for m = 1:length(methods)+1
    
    % Blank labels for ellipses
    for k = 1:N_ellipse+1 
        c = c + 1;
        legend_text{c} = '';
    end
    
    % Method label for paths
    c = c + 1;
    legend_text{c} = legend_names{m};
    
    % Blank labels for start/end points
    c = c + 1;
    legend_text{c} = '';
    c = c + 1;
    legend_text{c} = '';
    
end
    
% Append legend
legend_obj = legend(legend_text,'FontSize',24);

% Append axis labels and increase fontsize
xlabel('\mu_1');
ylabel('\mu_2');
set(gca,'FontSize',24);

% Adjust axis limits
xlim([-1.25,2.05]);
ylim([-0.7,1.25]);

% Make lines thicker
set(gca,'LineWidth',2);
set(legend_obj,'LineWidth',2);

% Return to original folder
cd('Figures');