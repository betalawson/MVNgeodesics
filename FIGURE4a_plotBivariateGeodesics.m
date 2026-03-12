function FIGURE4a_plotBivariateGeodesics
% This function plots a set of example geodesics through the space of
% bivariate normals

%%% SPECIFICATIONS

% Define the fixed length of all geodesics to be plotted
fix_length = true;
fixed_length = 2.5;

% Define the list of target points
T{1} = struct('mu', 4 * [cos(-pi/3); sin(-pi/3)], 'SIGMA',[cos(-pi/3),-sin(-pi/3);sin(-pi/3),cos(-pi/3)] * diag([1;4]) * [cos(pi/3),-sin(pi/3);sin(pi/3),cos(pi/3)]);
T{2} = struct('mu', [3.5; -0.5], 'SIGMA', [1.1, 1.0; 1.0, 1.5]);
T{3} = struct('mu', [2; 3.5], 'SIGMA', eye(2));
T{4} = struct('mu', [-3; 3], 'SIGMA', diag([0.25, 12]));
T{5} = struct('mu', [-1.5; -1.5], 'SIGMA', [1, 0; 0, 0.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the origin
O = struct('mu',zeros(2,1), 'SIGMA',eye(2));

% Prepare figure
figure('units','normalized','OuterPosition',[0 0 1 1]); hold on;

% Loop over each target point, and find the geodesic to it first
Gs = cell(1,length(T));
for k = 1:length(T)
    
    % Grab current target and find the geodesic to it
    Gs{k} = onePointShooting(O,T{k});
    % Normalise this geodesic to have the set length (if requested)
    if fix_length
        rescale_factor = fixed_length / Gs{k}.L;
        Gs{k}.v = struct('mu', Gs{k}.v.mu * rescale_factor, 'SIGMA', Gs{k}.v.SIGMA * rescale_factor);
    end
    
end

% Loop over each target point, drawing the covariance ellipses first
for k = 1:length(T)
       
    % Plot this geodesic
    plotGeodesic(Gs{k},'Npts',141,'covFrequency',20,'scale',0.2, 'pathColor', [0 0 0], 'ellipseColor', [0.5 0.5 0.5]);
    
end

% Now do another loop to plot only the paths themselves, and the MVNs
% (ensuring these are drawn over the covariance ellipses)
for k = 1:length(T)    

    % Plot this geodesic (path only)
    plotGeodesic(Gs{k},'Npts',141,'covFrequency',1e6,'scale',0.2, 'pathColor', [0 0 0], 'ellipseColor', [0.5 0.5 0.5], 'ellipseWidth', 1.75);
    % Plot the endpoint MVN
    plotMVN(fireGeodesic(Gs{k},1), 'ellipseWidth', 5, 'ellipseColor', [0 0 0], 'scale', 0.2, 'pointColor', [0 0 0], 'pointSize', 35);
    
end

% Plot the origin at the end
plotMVN(O, 'ellipseWidth', 4, 'ellipseColor', [0 0 0], 'scale', 0.2, 'pointColor', [0 0 0], 'pointSize', 35);
    
% Clean up the plot
axis([-3.5 3.5 -3.5 3.5]);
axis square;
xlabel('$\mu_1$','Interpreter','LaTeX');
ylabel('$\mu_2$','Interpreter','LaTeX');
set(gca,'FontSize', 24);
