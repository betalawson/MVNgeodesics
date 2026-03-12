function plotMVN( varargin )
%
%     plotMVN(p)
%     plotMVN(p, options)
%     plotMVN(ax, ...)
%
% This function plots the given multivariate normal distribution, expressed
% as a struct 'p' with fields 'mu' and 'SIGMA'. The user
% may optionally specify an axis to plot on, otherwise the current axis
% will be used (as per MATLAB's standard plot function). The plotting
% routine depends on the dimension, and a warning is output for dimensions
% greater than 3, for which no visualisation routine has been implemented.
% Note also that valid name-value pairs for plotting will depend on the
% dimension of the normal being plotted.

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

% Next input is always the MVN to be plotted
p = inputs{1};
inputs = inputs(2:end);


%%% ADDITIONAL INPUTS

% Create input parsing object
extra_options = inputParser();

% Add parameters to this object
addParameter(extra_options, 'pathStyle', '-');
addParameter(extra_options, 'ellipseColor', [0.4 0.4 0.4] );
addParameter(extra_options, 'ellipseWidth', 0.85) ;
addParameter(extra_options, 'ellipsoidFaceOpacity', 0.2);
addParameter(extra_options, 'ellipsoidEdgeOpacity', 0.2);
addParameter(extra_options, 'scale', 0.15);
addParameter(extra_options, 'pointSize', 20);
addParameter(extra_options, 'pointColor', [1, 0, 0]);

% Read out the parameter informaiton provided
parse(extra_options, inputs{:});
plot_options = extra_options.Results;



% Plotting depends on dimension
switch length(p.mu)
    
    % 1-D: A point in (mu, sigma^2) space
    case 1
        
        plot(ax,p.mu, p.SIGMA, '.', 'MarkerSize', plot_options.point_size, 'MarkerEdgeColor', plot_options.pointColor);
        
    % 2-D: A point in 2-D mu-space with covariance indicated by an ellipse
    case 2
        
        % Plot the mean point
        plot(ax,p.mu(1), p.mu(2), '.', 'MarkerSize', plot_options.pointSize, 'MarkerEdgeColor', plot_options.pointColor);
        % Make and plot the covariance matrix as an ellipse
        [x,y] = makeCovEllipse( p.mu(1), p.mu(2), p.SIGMA, plot_options.scale );
        plot(ax,x, y, 'LineWidth', plot_options.ellipseWidth, 'Color', plot_options.ellipseColor, 'LineStyle', plot_options.pathStyle );
        
    % 3-D: A point in 3-D mu-space with covariance indicated by ellipsoid
    case 3
        
        % Plot the mean point
        plot3(ax,p.mu(1), p.mu(2), p.mu(3), '.', 'MarkerSize', plot_options.pointSize, 'MarkerEdgeColor', plot_options.pointColor);
        % Find eigendecomposition of the covariance matrix
        [V, LAM] = eig(p.SIGMA);
        
        % Convert the variance amounts into "standard deviation"
        % amounts, scaled by the plotting scale
        STD = sqrt(diag(LAM)) * plot_options.scale;
        
        % Plot the ellipsoid in non-rotated, non-translated space
        [x,y,z] = ellipsoid( 0,0,0, STD(1), STD(2), STD(3) );
        matsize = size(x);
        X = [x(:),y(:),z(:)];
        
        % Transform the ellipsoid using rotation and translation
        X = (V * X' + p.mu)';
        x = X(:,1);
        y = X(:,2);
        z = X(:,3);
        x = reshape(x, matsize);
        y = reshape(y, matsize);
        z = reshape(z, matsize);
        
        % Plot this as a surface
        surf(ax,x,y,z, 'FaceAlpha', plot_options.ellipsoidFaceOpacity, 'EdgeAlpha', plot_options.ellipsoidEdgeOpacity, 'FaceColor', plot_options.ellipseColor);
        
    otherwise
        warning('Unable to plot multivariate normals in dimensions greater than three!');
        
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