function plotUnivariateGeodesics()
% This function plots many geodesics through the space of univariate
% normals, each starting from the origin of the space

% Specify the number of geodesics to plot
N_geos = 20;

% Define the origin
O = struct('mu',0,'SIGMA',1);

% Create a new figure
figure; hold on;

% Loop to plot the requested number of geodesics
for k = 1:N_geos
    
    % Place a random point to find a geodesic to
    P = struct('mu',5*rand,'SIGMA',exp(-2 + 4*rand));
    
    % Find the geodesic
    G = onePointShooting(O,P);
    
    % Plot the geodesic
    plotGeodesic(G);    
    
end

end