function plotMVN( p )
%
%     plotMVN( p )
%
% This function plots the given multivariate normal distribution, expressed
% as a struct with 'mu' component and 'SIGMA' component. The plotting
% routine depends on the dimension, and a warning is output for dimensions
% greater than 3, for which no visualisation routine has been implemented.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale = 0.15;           % Controls the size of ellipses for 2-D and 3-D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting depends on dimension
switch length(p.mu)
    
    % 1-D: A point in (mu, sigma^2) space
    case 1
        
        plot(p.mu, p.SIGMA, 'k.', 'MarkerSize', 20);
        
    % 2-D: A point in 2-D mu-space with covariance indicated by an ellipse
    case 2
        
        % Plot the mean point
        plot(p.mu(1), p.mu(2), 'k.', 'MarkerSize', 20);
        % Make and plot the covariance matrix as an ellipse
        [x,y] = makeCovEllipse( p.mu(1), p.mu(2), p.SIGMA, scale );
        plot(x, y, 'k', 'LineWidth', 2);
        
    % 3-D: A point in 3-D mu-space with covariance indicated by ellipsoid
    case 3
        
        % Plot the mean point
        plot3(p.mu(1), p.mu(2), p.mu(3), 'k.', 'MarkerSize', 20);
        % Find eigendecomposition of the covariance matrix
        [V, LAM] = eig(p.SIGMA);
        
        % Convert the variance amounts into "standard deviation"
        % amounts, scaled by the plotting scale
        STD = sqrt(diag(LAM)) * scale;
        
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
        surf(x,y,z, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2, 'FaceColor', [0.5 0.5 1]);
        
    otherwise
        warning('Unable to plot multivariate normals in dimensions greater than three!');
        
end

