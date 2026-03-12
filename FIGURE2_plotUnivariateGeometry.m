function FIGURE2_plotUnivariateGeometry
% This function plots many geodesics through the space of univariate
% normals, each starting from the origin of the space

% Specify whether to compare these against Euclidean geodesics or not
compare_wasserstein = true;
show_tangent = false;

% Specify the number of geodesics to plot with thick lines
N_coords_thick = 5;
% Specify the number of geodesics to plot with thin lines
N_lines_thin = 201;
% Number of points for plotting Wasserstein geodesics
N_Wpts = 101;

% Select the (i,j) values for which to draw in the tangent vector
%   (These can be set to -1 if don't want to plot)
tangent_i = 3;
tangent_j = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the origin
O = struct('mu',0,'SIGMA',1);

% Create vectors of targets (thick line)
mean_max = 4;
target_means = linspace(-mean_max, mean_max, N_coords_thick);
var_max = 5;
target_vars = linspace(1, var_max, N_coords_thick);

% Use the maximum co-ordinate difference to set the scaling for thin lines
% (they should extend outside the plot)
scale = 1.25*max([mean_max-0; var_max-1]);

% Create a vector of angles for thin lines
thetas = linspace(0,pi,N_lines_thin);

% Create a new figure
if compare_wasserstein
    figure;
else
    figure; hold on;
end

% Loop to plot the requested number of geodesics
t = linspace(0, scale, N_Wpts);
for i = 1:N_lines_thin
    
    % Special processing for the Euclidean plot if requested
    if compare_wasserstein
        % Move to Euclidean subplot
        subplot(1,2,1); hold on;
        % Find the velocity corresponding to the current angle (just use
        % angles in standard deviation space for ease)
        vmu = cos(thetas(i));
        vsigma = sqrt(1 - vmu^2);
        % Fire the Wasserstein geodesic (straight line in std. dev. space)
        mu_W = t * vmu;
        sigma2_W = (1 + t * vsigma).^2;
        % Also get negative variance lines - set negatives to zero before squaring
        sigma2_Wminus = (1 - t * vsigma);
        sigma2_Wminus(sigma2_Wminus < 0) = 0;
        sigma2_Wminus = sigma2_Wminus.^2;
        % Plot these curves
        plot(mu_W, sigma2_W, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5);
        plot(mu_W, sigma2_Wminus, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5);
        % Move back to Fisher-Rao subplot
        subplot(1,2,2); hold on;
    end
    
    % find the velocity corresponding to the current angle
    vmu = cos(thetas(i));
    vsigma2 = sqrt(2)*sqrt(1 - vmu^2);
    % Plot the geodesic with this velocity, scaled
    G = struct('P',1,'r',0,'v',struct('mu',vmu*scale,'SIGMA',vsigma2*scale));
    plotGeodesic(G, 'PathColor', [0.75 0.75 0.75], 'pathWidth', 1.5);
    % Also plot the equivalent geodesic but moving downwards in variance
    G = struct('P',1,'r',0,'v',struct('mu',vmu*scale,'SIGMA',-vsigma2*scale));
    plotGeodesic(G, 'PathColor', [0.75 0.75 0.75], 'pathWidth', 1.5);
    
end


% Loop to plot the requested number of geodesics
t = linspace(0, 1, N_Wpts);
for i = 1:N_coords_thick
    for j = 1:N_coords_thick
    
        % Special processing for the Euclidean plot if requested
        if compare_wasserstein
            % Move to Euclidean subplot
            subplot(1,2,1); hold on;
            % Plot a Wasserstein geodesic
            mu_W = t * target_means(i);
            sigma2_W = (1 + t * (sqrt(target_vars(j)) - 1)).^2;
            plot(mu_W, sigma2_W, 'Color', [0 0 0], 'LineWidth', 4.5);
            % Plot the target point
            plot(target_means(i),target_vars(j),'k.','MarkerSize',40);
            % Move back to Fisher-Rao subplot
            subplot(1,2,2); hold on;
        end
        
        % Place grid point and find the geodesic to it
        P = struct('mu',target_means(i),'SIGMA',target_vars(j));
        G = onePointShooting(O,P);
    
        % Plot the geodesic and target point
        plotGeodesic(G, 'PathColor', [0 0 0], 'pathWidth', 4.5);
        plot(target_means(i),target_vars(j),'k.','MarkerSize',40);
   
    end
end

% If plotting the tangent lines, do so at the end so they overwrite
% existing objects
if show_tangent   
    
    % Special processing for the Euclidean plot if requested
    if compare_wasserstein
        % Move to Euclidean subplot
        subplot(1,2,1); hold on;
        % Read out the target point for which to show a tangent
        tangent_x = target_xs(tangent_i);
        tangent_y = target_ys(tangent_j);
        % Plot the Euclidean tangent (straight line) and re-colour point
        plot([0, tangent_x], [0, tangent_y], 'Color', [1.0, 0.3, 0.3], 'LineWidth', 8);
        plot(tangent_x,tangent_y,'.','MarkerEdgeColor',[0.8 0.0 0.0],'MarkerSize',65);
        % Move back to Fisher-Rao subplot
        subplot(1,2,2); hold on;
    end
    
    % Read out the target point for which to show a tangent
    tangent_mean = target_means(tangent_i);
    tangent_var = target_vars(tangent_j);
    % Prepare target point
    P = struct('mu',tangent_mean,'SIGMA',tangent_var);
    % Re-find the geodesic for this point
    G = onePointShooting(O,P,struct('symKL_tol',1e-12));
    % Plot its tangent vector using geodesic velocity, also re-colour point
    plot([0, G.v.mu], [1, 1+G.v.SIGMA], 'Color', [1.0, 0.3, 0.3], 'LineWidth', 8);
    plot(tangent_mean,tangent_var,'.','MarkerEdgeColor',[0.8 0.0 0.0],'MarkerSize',65);
    
end

% Add an origin to all plots
if compare_wasserstein
    subplot(1,2,1);
    plot(0,1,'w.','MarkerSize',50);
    subplot(1,2,2);
end
plot(0,1,'w.','MarkerSize',50);

% Plot clean up
if compare_wasserstein
    subplot(1,2,1); hold on;
    axis([-mean_max, mean_max, 1-0.1*var_max, var_max*1.1]);   
    xlabel('\mu');
    ylabel('\sigma^2');
    set(gca,'FontSize',24,'LineWidth',2);
    title('2-Wasserstein Geometry');
    subplot(1,2,2); hold on;
end
xlabel('\mu');
ylabel('\sigma^2');
set(gca,'FontSize',24,'LineWidth',2);
title('Fisher Geometry');
axis([-mean_max, mean_max, 1-0.1*var_max, var_max*1.1]);
% Ensure axis lines are top of the plot
if compare_wasserstein
    subplot(1,2,1);
    set(gca,'Layer','top');
    subplot(1,2,2);
end
set(gca,'Layer','top');
end