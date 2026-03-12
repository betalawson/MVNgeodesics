function FIGURE5_geodesicsPlane
% This function shows a set of geodesics on the space of trivariate
% normals, highlighting the spaces of points that can or cannot be reached
% by a straight-line geodesic when one eigenvalue of the target covariance
% is repeated.


%%% SPECIFICATIONS

% Random seed for reproducibility
rng(7);

% Starting location
P0 = struct('mu', [3;-1;2], 'SIGMA', [4.0, 0.4, -0.4; 0.4, 4.2, 1.0; -0.4, 1.0, 2.0]);

% Target covariance (constructed to have a repeated eigenvalue)
alpha = pi/8;
beta = -pi/9;
V_star = [ cos(alpha)*cos(beta), -sin(alpha),  cos(alpha)*sin(beta);
           sin(alpha)*cos(beta),  cos(alpha),  sin(alpha)*sin(beta);
                     -sin(beta),           0,  cos(beta)             ];
SIGMA_star = V_star * diag([5.5; 5.5; 0.5]) * V_star';

% Specify the number of points to plot in the plane, and off the plane/line
Npts = 5;

% Specify min/max Euclidean lengths of vectors in mu-space to plane points
L_min = 0.75;
L_max = 1;
% Specify min/max shift in plane perpendicular direction for "off" points
delta_min = 0.5;
delta_max = 1;

% Specify the size of the plane to plot
plane_size = 1;


%%% PLANE PREPARATION

% Generate vectors that move along the plane, defined by the two
% repeated-eigenvalue eigenvectors, but transformed to lie in actual space
V_plane = P0.SIGMA^(1/2) * V_star(:,1:2);

% Generate corners (in plane parameter space) for a slice of the plane
c1 = plane_size*[-1, -1, 1, 1];
c2 = plane_size*[-1, 1, 1, -1];

% Convert corners to lie in actual space
C = V_plane * [c1;c2] + P0.mu;

% Find the perpendicular vector to the plane (visually perpendicular)
V_perp = cross(V_plane(:,1), V_plane(:,2));
V_perp = V_perp / norm(V_perp);

%%% TARGET POINT GENERATION (TARGETS IN PLANE)

% Create angles for the requested number of plane points around the circle
t = linspace(0,2*pi,Npts+1);
t = t(1:Npts);

% Convert these into vectors (in plane parameter space) of random length
M_plane = ( L_min + (L_max - L_min) * rand(1,Npts) ) .* [cos(t); sin(t)];

% Convert target covariance from origin-transformed space to actual space
P1SIGMA = P0.SIGMA^(1/2) * SIGMA_star * P0.SIGMA^(1/2);

% Convert from plane parameter space to actual target points
P1s_plane = cell(1,Npts);
for k = 1:Npts
    P1s_plane{k} = struct('mu', P0.mu + V_plane*M_plane(:,k), 'SIGMA', P1SIGMA);
end


%%% TARGET POINT GENERATION (ALIGNED WITH REMAINING EIGENVECTOR)

% Target mean is a step along remaining eigenvector (shifted)
P1_line = struct('mu', 1.2*P0.SIGMA^(1/2) * V_star(:,3) + P0.mu, 'SIGMA', P1SIGMA);


%%% TARGET POINT GENERATION (NOT FALLING ON THE PLANE OR LINE)

% Take the existing target points and shift them in perpendicular direction
% (allowing for better visualisation)
P1s_off = cell(1,Npts);
for k = 1:Npts
    P1s_off{k} = struct('mu', P1s_plane{k}.mu + 2 * (delta_min + (delta_max - delta_min) * rand) * V_perp, 'SIGMA', P1SIGMA);
end


%%% GEODESIC FINDING TO ALL POINTS

% The set of targets on the plane, and the points above them
for k = 1:Npts
    Gs_plane{k} = onePointShooting(P0, P1s_plane{k});
    Gs_off{k} = onePointShooting(P0, P1s_off{k});
end

% The line target
G_line = onePointShooting(P0, P1_line);


%%% PLOTTING

% Initialise figure
figure('units', 'normalized', 'OuterPosition', [0 0 1 1]);
hold on;

% Plot the plane as the patch between these points
patch( C(1,:), C(2,:), C(3,:), [0.5 0.5 0.5], 'faceAlpha', 0.85 );

% Plot the geodesics on the plane
for k = 1:Npts
    plotGeodesic(Gs_plane{k}, 'ellipseColor', [0.75 0.75 0.5], 'pathColor', [1.0 0.9 0.2], 'Npts', 101, 'covFrequency', 100, 'ellipsoidEdgeOpacity', 0.75, ...
                              'ellipsoidFaceOpacity', 0.75, 'ellipsoidEdgewidth', 0.75, 'ellipsoidPoints', 15, 'pathWidth', 4);
end

% Plot a dotted line along the third eigenvector
line_lim = [-0.75 1.75];
line_ends = line_lim .* (P0.SIGMA^(1/2) * V_star(:,3)) + P0.mu;
plot3([line_ends(1,1), P0.mu(1)], [line_ends(2,1), P0.mu(2)], [line_ends(3,1), P0.mu(3)], 'k--', 'LineWidth', 3.5);   % Line bottom to "origin"
plot3([line_ends(1,2), P1_line.mu(1)], [line_ends(2,2), P1_line.mu(2)], [line_ends(3,2), P1_line.mu(3)], 'k--', 'LineWidth', 3.5); % Geodesic end to line top

% Plot the geodesic on the line
plotGeodesic(G_line, 'ellipseColor', [0.75 0.75 0.5], 'pathColor', [1.0 0.9 0.2], 'Npts', 101, 'covFrequency', 100, 'ellipsoidEdgewidth', 0.75, 'ellipsoidEdgeOpacity', 0.75, ...
                     'ellipsoidFaceOpacity', 0.75, 'ellipsoidPoints', 15, 'pathWidth', 4);

% Plot the geodesics off the plane
for k = 1:Npts
    plotGeodesic(Gs_off{k}, 'ellipseColor', [0.5 0.5 0.75], 'pathColor', [0.2 0.2 0.8], 'Npts', 101, 'covFrequency', 100, 'ellipsoidEdgewidth', 0.75, 'ellipsoidEdgeOpacity', 0.75, ...
                            'ellipsoidFaceOpacity', 0.75, 'ellipsoidPoints', 15, 'pathWidth', 4);
end

% Plot the line over this straight-line geodesic just to improve the graphical appearance of the figure
plot3([P0.mu(1) P1_line.mu(1)], [P0.mu(2) P1_line.mu(2)], [P0.mu(3) P1_line.mu(3)], 'LineWidth', 4, 'Color', [1.0 0.9 0.2]);

% Turn off axes and set viewing angle
axis equal;
axis off;
xlim([-1 7]);
ylim([-5.25 3.25]);
zlim([0.5 4.25]);
set(gca,'view',[-51.9923, -3.5311]);

end