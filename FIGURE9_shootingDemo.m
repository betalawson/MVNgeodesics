function FIGURE9_shootingDemo
% This function creates a figure that demonstrates the concept of shooting
% and using parallel transport to form velocity corrections out of
% "residual" vectors

%%% SPECIFICATIONS

% Start and end MVN for the example
P0 = struct('mu',[1;1],'SIGMA',diag([5;0.5]));
P1 = struct('mu',[3;3],'SIGMA',diag([0.5;5]));

% Specify axis limits to contain the figure
x_limits = [0.4, 3.25];
y_limits = [0.75, 3.6];

% Scaling for showing vector mean components 
vel_scale = 0.5;
vec_width = 2.75;     % Line width

% Number of arrows used to demonstrate parallel transport
N_arrows = 16;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise figure
figure('units', 'normalized', 'OuterPosition', [0.2625, 0.0389, 0.4969, 0.9204]);
hold on;

% Clean up and label the figure
set(gca, 'FontSize', 24);
xlim(x_limits);
ylim(y_limits);
xlabel('\mu_1');
ylabel('\mu_2');
axis square;
set(gca,'LineWidth',2);

% Transform problem to the origin so parallel transport equations work
O = struct('mu',zeros(2,1),'SIGMA',eye(2));
T = struct('mu',P0.SIGMA^(-1/2)*(P1.mu-P0.mu),'SIGMA',P0.SIGMA^(-1/2)*P1.SIGMA*P0.SIGMA^(-1/2));

% Generate and plot the approximate geodesic connecting these points
G = approximateGeodesic(O,T,'eigen');

% Generate the endpoint reached by this geodesic
P_G = fireGeodesic(G,1);
invSIGMA = P_G.SIGMA \ eye(size(P_G.SIGMA));

% Find approximate geodesic connecting geodesic endpoint and target point
G_c = approximateGeodesic(P_G, T, 'eigen');

% Find the euclidean connection from p to pt too
R = struct('mu', T.mu - P_G.mu, 'SIGMA', T.SIGMA - P_G.SIGMA);

% Grab connecting velocity (in origin-transformed space)
w = affineTransform( G_c.v, G_c.P, zeros(2,1) );

% Parallel transport connecting geodesic's velocity back to start point
dv = parallelTransport(w,G,[1 0]);

% Find Jacobi field associated with this shift to geodesic's endpoint
J = jacobiField(G,dv,1);
        
% Project the residual onto the Jacobi field to select a step length
s = innerProduct(R,J,invSIGMA) / innerProduct(J,J,invSIGMA);

% Apply velocity update with this step length
G_up = G;
G_up.v = struct('mu',G.v.mu + s*dv.mu,'SIGMA',G.v.SIGMA + s*dv.SIGMA);

% Create the "base" geodesic in the original space and plot it
Gbase.P = P0.SIGMA^(1/2)*G.P;
Gbase.r = P0.mu;
Gbase.v = G.v;        % Geodesic objects always describe velocity at origin (so no transform needed)
plotGeodesic(Gbase, 'pathColor', [1 0.3 0.3], 'ellipseColor', [1 0.6 0.6], 'pathWidth', 5, 'endColor', [1 0.3 0.3], 'ellipseWidth', 2, 'Npts', 1001, 'covFrequency', 200, 'endColor', [1 0.3 0.3], 'pointSize', 30, 'scale', 0.1);

% Create the updated geodesic in the original space and plot it
G_up.P = P0.SIGMA^(1/2)*G.P;
G_up.r = P0.mu;
plotGeodesic(G_up, 'pathColor', [0.3 0.6 1], 'ellipseColor', [0.6 0.8 1], 'pathWidth', 5, 'endColor', [0.3 0.6 1], 'ellipseWidth', 2, 'Npts', 1001, 'covFrequency', 200, 'endColor', [0.3 0.6 1], 'pointSize', 30, 'scale', 0.1);

% Create the connecting geodesic in the original space and plot it
G_c.P = P0.SIGMA^(1/2)*G_c.P;
G_c.r = P0.mu + P0.SIGMA^(1/2)*G_c.r;
plotGeodesic(G_c, 'pathColor',[0.4 0.4 0.4], 'ellipseColor', [0.4 0.4 0.4], 'pathWidth', 3, 'endColor', [0.4 0.4 0.4], 'ellipseWidth', 2, 'Npts', 1001, 'covFrequency', -1, 'endColor', [0 0 1], 'pointSize', 0.001, 'scale', 0.1);

% Use a loop to plot the parallel transported vector at regular intervals
t_vec = linspace(1,0,N_arrows);
for k = 1:N_arrows
    
    % Find where the reduced geodesic ends up (in origin-transformed space)
    Tk = fireGeodesic(G,t_vec(k));
    % Perform partial backwards parallel transport along original geodesic
    wk = parallelTransport(w,G,[1,t_vec(k)]);
    
    % Convert both objects back into the original domain
    %  (Only doing mean component because that's all we plot)
    wk.mu = P0.SIGMA^(1/2)*wk.mu;           % Tangent transforms P mu
    Tk.mu = P0.SIGMA^(1/2)*Tk.mu + P0.mu;   % Position transforms P mu + r
    
    % Plot the vectors as lines plus triangle arrowheads
    plotArrow(Tk.mu(1), Tk.mu(2), vel_scale*wk.mu(1), vel_scale*wk.mu(2), vec_width, [0.4, 0.4, 0.4], true);
        
end

% Plot the mean components of the initial velocity vectors too
v1_pos = P0.mu;
v1_vel = vel_scale * Gbase.P * Gbase.v.mu;
v2_pos = P0.mu;
v2_vel = vel_scale * G_up.P * G_up.v.mu;
plotArrow(v1_pos(1), v1_pos(2), v1_vel(1), v1_vel(2), vec_width, [1 0.3 0.3], true);
plotArrow(v2_pos(1), v2_pos(2), v2_vel(1), v2_vel(2), vec_width, [0.3 0.6 1.0], true);
% Plot the correction vector over here too
plotArrow(v1_pos(1)+v1_vel(1), v1_pos(2)+v1_vel(2), vel_scale*wk.mu(1), vel_scale*wk.mu(2), vec_width, [0.4 0.4 0.4], true);
plotArrow(v1_pos(1)+v1_vel(1), v1_pos(2)+v1_vel(2), vel_scale*s*wk.mu(1), vel_scale*s*wk.mu(2), vec_width, [0.4 0.4 0.4], false);

% Plot the targets
plotMVN(P0, 'ellipseColor',[0 0 0], 'ellipseWidth', 4, 'pointColor', [0 0 0], 'pointSize', 35, 'scale', 0.1);
plotMVN(P1, 'ellipseColor',[0 0 0], 'ellipseWidth', 4, 'pointColor', [0 0 0], 'pointSize', 35, 'scale', 0.1);

% Plot clean-up


end

function plotArrow(x,y,u,v,line_width,arrow_clr, dashed)

% Assume no dashing on the line unelss specifically asked for it
if nargin < 7
    dashed = false;
end

L = 0.05;

% Calculate angle orientation of provided vector
theta = atan2(v,u);

% Plot the line part of the vector - dashed if requested
if dashed
    plot(x + [0, u - L/2*cos(theta)], y + [0, v - L/2*sin(theta)], ':', 'LineWidth', line_width, 'Color', arrow_clr);
else
    plot(x + [0, u - L/2*cos(theta)], y + [0, v - L/2*sin(theta)], '-', 'LineWidth', line_width, 'Color', arrow_clr);
end

% Generate three points for the arrowhead
p(1,:) = [x + u, y + v];
p(2,:) = [x + u - L * cos(theta-pi/8), y + v - L * sin(theta-pi/6)];
p(3,:) = [x + u - L * cos(theta+pi/8), y + v - L * sin(theta+pi/6)];

patch(p(:,1), p(:,2), arrow_clr, 'EdgeColor', 'none');

end