function FIGURE4b_plotTrivariateGeodesics
% This function plots a set of example geodesics through the space of
% bivariate normals

% Define the fixed length of all geodesics to be plotted
fix_length = true;
fixed_length = 2.0;

% Define the list of initial velocities
v{1} = struct('mu', [0.2; -0.2; -0.1], 'SIGMA', 1.05*[-0.3, 0.2, -0.2; 0.2, -0.3, 0.3; -0.2, 0.3, 1.2]);
v{2} = struct('mu', [0.2; -0.4; 0], 'SIGMA', 1.5*[0.1, -0.2, 0.5; -0.2, 0.5, -0.1; 0.5, -0.1, 0.2]);
v{3} = struct('mu', [-0.3; 0.2; 0.1], 'SIGMA', [0.2, 0.1, 0.1; 0.1, 0.3, 0.1; 0.1, 0.1, -0.3]);
v{4} = struct('mu', [0; 0.6; 0.3], 'SIGMA', [0.4, -0.2, -0.2; -0.2, -0.6, 0.2; -0.2, 0.2, 0.1]);
% Generate a velocity for which the matrices uu^T and U share eigenvectors
%   (this is a way to create a straight line geodesic)
u = [0.5;0.2;0.9]; V = randomOrthogonalMatrix(3,u); U = V * diag([0.2; 0.4; 0.2]) * V';
v{5} = struct('mu', u, 'SIGMA', U);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the origin
O = struct('mu',zeros(3,1), 'SIGMA',eye(3));

% Prepare figure
figure('units','normalized','OuterPosition',[0.2620    0.1722    0.5042    0.6722]); hold on;

% Loop over each velocity, and fire its corresponding geodesic
Gs = cell(1,length(v));
for k = 1:length(v)
    
    % Normalise length if that flag was set
    if fix_length
        L = innerProduct(v{k},v{k},eye(3));
        rescale_factor = fixed_length / L;
        Gs{k} = struct('P',eye(3),'r',zeros(3,1),'v',struct('mu', v{k}.mu * rescale_factor, 'SIGMA', v{k}.SIGMA * rescale_factor));
    end
    
end

% Loop over each target point, drawing the covariance ellipses first
for k = 1:length(v)
       
    % Plot this geodesic
    plotGeodesic(Gs{k},'Npts',14000,'covFrequency',2000,'scale',0.2, 'pathColor', [0 0 0], 'ellipseColor', [0.85 0.85 0.85],'ellipsoidFaceOpacity', 0.25, 'ellipsoidEdgeOpacity', 0.1);
    
end

% Now do another loop to plot only the paths themselves, and the MVNs
% (ensuring these are drawn over the covariance ellipses)
for k = 1:length(v)    

    % Plot this geodesic (path only)
    plotGeodesic(Gs{k},'Npts',141,'covFrequency',1e6,'scale',0.2, 'pathColor', [0 0 0], 'ellipseColor', [0.5 0.5 0.5], 'ellipseWidth', 1.75);
    % Plot the endpoint MVN
    plotMVN(fireGeodesic(Gs{k},1), 'ellipseColor', [0.45 0.45 0.45], 'scale', 0.2, 'pointColor', [0 0 0], 'pointSize', 35, 'ellipsoidFaceOpacity', 0.5, 'ellipsoidEdgeOpacity', 0.5);
    
end

% Plot the origin at the end
plotMVN(O, 'ellipseColor', [0.5 0.5 0.5], 'scale', 0.2, 'pointColor', [0 0 0], 'pointSize', 35, 'ellipsoidFaceOpacity', 1, 'ellipsoidEdgeOpacity', 1);

    
% Clean up the plot
xlim = 2;
ylim = 2;
zlim = 1.5;
axis equal;
axis([ xlim*[-1,1], ylim*[-1,1], zlim*[-1,1] ]);
xlabel('$\mu_1$','Interpreter','LaTeX');
ylabel('$\mu_2$','Interpreter','LaTeX');
zlabel('$\mu_3$','Interpreter','LaTeX');
set(gca,'FontSize', 24);
set(gca,'View', [-56.7716   10.3933]);
% Set up the ticks to add grid lines separate to the ticks on axes
tick_vals{1} = -xlim:0.5:xlim;
tick_vals{2} = -ylim:0.5:ylim;
tick_vals{3} = -zlim:0.5:zlim;
tick_labels = cell(3,1);
for d = 1:3
    tick_labels{d} = cell(1,length(tick_vals{d}));
    for k = 1:length(tick_vals{d})
        if abs(tick_vals{d}(k) - round(tick_vals{d}(k)) ) < 1e-8
            tick_labels{d,k} = num2str(tick_vals{d}(k));
        else
            tick_labels{d,k} = '';
        end
    end
end
xticks(tick_vals{1});
yticks(tick_vals{2});
zticks(tick_vals{3});
xticklabels(tick_labels(1,:));
yticklabels(tick_labels(2,:));
zticklabels(tick_labels(3,:));
xtickangle(0);
ytickangle(0);
ztickangle(0);
grid on;
% Thicken gridlines, and thicken axes more (separately)
set(gca,'LineWidth',1.5);
plot3([-2 2;],[-2 -2],[-1.5 -1.5],'k-', 'LineWidth', 2.5);
plot3([-2 -2;],[-2 2],[-1.5 -1.5],'k-', 'LineWidth', 2.5);
plot3([-2 -2;],[2 2],[-1.5 1.5],'k-', 'LineWidth', 2.5);
