function [x, y] = makeCovEllipse( xc, yc, SIGMA, scale )
%
% This internal function creates a series of (x,y) pairs of points that
% sketch out an ellipse with centre (xc, yc), orientation and axis lengths
% defined by the input positive-definite covariance matrix SIGMA, and with
% axis lengths also optionally scaled by an additional factor

% Default to no scaling if no scale factor provided
if nargin < 4
    scale = 1;
end

% Create equally spaced points around a circle
t = linspace(0,2*pi,5000);

% Find eigendecomposition of the covariance matrix
[V, LAM] = eig(SIGMA);
LAM = diag(LAM);

% Check for a valid eigendecomposition
if all( LAM >= 0 )

    % Define ellipse in polar co-ordinates
    r = scale ./ sqrt( cos(t - atan2(V(2,1),V(1,1))).^2 / LAM(1)^2 + sin(t - atan2(V(2,1),V(1,1))).^2 / LAM(2)^2 );

    % Map to (x,y) points with requested centre
    x = xc + r .* cos(t);
    y = yc + r .* sin(t);
    
else
    
    warning('Unable to represent this covariance as an ellipse, will not appear on plot');
    x = [];
    y = [];
    
end