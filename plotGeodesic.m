function plotGeodesic(varargin)
%
%     plotGeodesic(G, options)
%     plotGeodesic(G, range, options)
%     plotGeodesic(ax, ...)
%
% This function plots the provided geodesic object 'G'. If no range is
% specified, the geodesic will be plotted over the range [0,1]. The user
% may optionally specify an axis to plot on, otherwise the current axis
% will be used (as per MATLAB's standard plot function).
%
% Options should be provided as name-value pairs, and include:
%
%             Npts - Number of points to plot along the given range
%        pathColor - Color of the line through mu-space on the plot
%        pathWidth - Width of the line through mu-space on the plot
%     ellipseColor - Color of the covariance ellipses on the plot
%     ellipseWidth - Width of the covariance ellipses on the plot
%     covFrequency - How many points between each covariance ellipse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INPUT HANDLING

% First store all inputs just to separate from 'varargin'
inputs = varargin;

% Check if the first input was an axis object
if isAxes( inputs{1} )
    % Use provided axes
    ax = inputs{1};
    % Remove axes from the input set
    inputs = inputs(2:end);
    
else
    
    % Use current axes (or create a new figure if one is not present)
    ax = gca;
    
end

% Next input is always the geodesic
G = inputs{1};
inputs = inputs(2:end);

% Check if the next input is a range, if not assume the range is [0,1]
range_given = false;
if ~isempty(inputs) > 0
    if ismatrix(inputs{1}) && length(inputs{1}) == 2
        range_given = true;
    end
end

if range_given
    range = inputs{1};
    inputs = inputs(2:end);
else
    range = [0, 1];
end



%%% ADDITIONAL INPUTS

% Create input parsing object
extra_options = inputParser();

% Add parameters to this object
addParameter(extra_options, 'lineSpec', '-');
addParameter(extra_options, 'Npts', 51);
addParameter(extra_options, 'ellipseColor', [0.4 0.4 0.4] );
addParameter(extra_options, 'ellipseWidth', 0.85) ;
addParameter(extra_options, 'pathColor', [0 0 1] );
addParameter(extra_options, 'pathWidth', 2 );
addParameter(extra_options, 'pathStyle', '-' );
addParameter(extra_options, 'scale', 0.15);
addParameter(extra_options, 'pointSize', 20);
addParameter(extra_options, 'startColor', [1, 0, 0]);
addParameter(extra_options, 'endColor', [0 0 0]);
addParameter(extra_options, 'covFrequency', 10);

% Read out the parameter informaiton provided
parse(extra_options, inputs{:});
plot_options = extra_options.Results;


%%% GEODESIC EVALUATION

% Use geodesic tracing function to get a set of (mu, SIGMA) values along the path
path = traceGeodesic( G, range, plot_options.Npts );

%%% PLOTTING

% Hold the axes
hold(ax, 'on');

% Split according to the dimension of the geodesic
switch length(path{1}.mu)
    
    
    % In one dimension, plot mu against sigma squared on a 2D plot
    case 1
        
        % Run along the geodesic, storing the mean and variance values
        mu_path = zeros(length(path),1);
        for k = 1:length(path)
            
            % Store the mean value
            mu_path(k) = path{k}.mu';
            % Store the covariance value
            SIGMA_path(k) = path{k}.SIGMA;
            
        end
        
        % Plot mu versus SIGMA
        plot(ax, mu_path, SIGMA_path, 'color', plot_options.pathColor, 'LineWidth', plot_options.pathWidth, 'LineStyle', plot_options.pathStyle );
        
        
        
        
    % In two dimensions, plot the movement of the mean on a 2D plot, with
    % ellipses indicating covariance
    case 2
        
        % Run along the geodesic, storing the mean path and plotting covariances
        mu_path = zeros(length(path),2);
        for k = 1:length(path)
            
            % Plot the covariance ellipse if frequency is hit
            if mod(k-1,plot_options.covFrequency) == 0 && plot_options.covFrequency >= 0
                [x,y] = makeCovEllipse( path{k}.mu(1), path{k}.mu(2), path{k}.SIGMA, plot_options.scale );
                plot(ax, x, y, 'Color', plot_options.ellipseColor, 'LineWidth', plot_options.ellipseWidth);
            end
            
            % Store the mean vector
            mu_path(k,:) = path{k}.mu';
            
        end
        
        % Plot the mu path
        plot(ax, mu_path(:,1), mu_path(:,2), plot_options.lineSpec, 'color', plot_options.pathColor, 'LineWidth', plot_options.pathWidth, 'LineStyle', plot_options.pathStyle  );
        % Plot the starting and ending points
        plot(ax, mu_path(1,1), mu_path(1,2), '.', 'MarkerSize', plot_options.pointSize, 'MarkerEdgeColor', plot_options.startColor);
        plot(ax, mu_path(end,1), mu_path(end,2), '.', 'MarkerSize', plot_options.pointSize, 'MarkerEdgeColor', plot_options.endColor);
        
        
        
        
    % In three dimensinos, plot the movement of the mean on a 3D plot, with
    % ellipsoids indicating covariance
    case 3
        
        % Run along the geodesic, storing the mean path and plotting covariances
        mu_path = zeros(length(path),3);
        for k = 1:length(path)
            
            % Plot the covariance ellipse if frequency is hit
            if mod(k-1,plot_options.covFrequency) == 0 && plot_options.covFrequency >= 0
                
                % Find eigendecomposition of the covariance matrix
                [V, LAM] = eig(path{k}.SIGMA);
                
                % Convert the variance amounts into "standard deviation"
                % amounts, scaled by the plotting scale
                STD = sqrt(diag(LAM)) * plot_options.scale;
                
                % Plot the ellipsoid in non-rotated, non-translated space
                [x,y,z] = ellipsoid( 0,0,0, STD(1), STD(2), STD(3) );
                matsize = size(x);
                X = [x(:),y(:),z(:)];
                
                % Transform the ellipsoid using rotation and translation
                X = (V * X' + path{k}.mu)';
                x = X(:,1);
                y = X(:,2);
                z = X(:,3);
                x = reshape(x, matsize);
                y = reshape(y, matsize);
                z = reshape(z, matsize);
                
                % Plot this as a surface
                surf(x,y,z, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'FaceColor', plot_options.ellipseColor);
                
            end
            
            % Store the mean vector
            mu_path(k,:) = path{k}.mu';
            
        end
        
        % Plot the mu path
        plot3(ax, mu_path(:,1), mu_path(:,2), mu_path(:,3), 'color', plot_options.pathColor, 'LineWidth', plot_options.pathWidth, 'LineStyle', plot_options.pathStyle  );
        % Plot the starting and ending points
        plot3(ax, mu_path(1,1), mu_path(1,2), mu_path(1,3), '.', 'MarkerSize', plot_options.pointSize, 'MarkerEdgeColor', plot_options.startColor);
        plot3(ax, mu_path(end,1), mu_path(end,2), mu_path(end,3), '.', 'MarkerSize', plot_options.pointSize, 'MarkerEdgeColor', plot_options.endColor);
        
        
        
        
    % Otherwise, do nothing except warning the user
    otherwise
        warning('Unable to visualise geodesics in dimensions greater than three!');
        
end


end



% A simple subfunction that checks if the given input is an axes object
function output = isAxes(input)
try
    output = strcmp(get(input, 'type'), 'axes');
catch
    output = false;
end
end