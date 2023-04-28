function FIGURES3_pathRefinement
% This function visualises a step of path refinement for instructive
% purposes

% Return to the main directory
cd('..');

% Specify number of points to use in the example
N = 9;

% Define the problem to use as an example
P1 = struct('mu',[0;0],'SIGMA',[10,0;0,1]);
P2 = struct('mu',[5;5],'SIGMA',[1,0;0,10]);

% Form the Euclidean path between these
pts = closedFormPath(P1,P2,N,'euclid');
new_pts = pts;

% Join up the even index points with geodesics
for k = 3:2:N-1
    G{k} = onePointShooting(pts{k-1},pts{k+1});
    new_pts{k} = fireGeodesic(G{k},0.5);
end

% Join up the new odd-index points with geodesics
for k = 2:2:N-1
    G{k} = onePointShooting(new_pts{k-1},new_pts{k+1});
    new_pts{k} = fireGeodesic(G{k},0.5);
end

% Also find the true geodesic
trueG = onePointShooting(P1,P2);

% Prepare figure
figure; hold on;

% Plot the geodesics used
for k = 2:N-1
    if mod(k,2) == 0
        plotGeodesic(G{k},'Npts',101,'covFrequency',1001,'pathColor',[0 0 1],'startColor',[0 0 1],'endColor',[0 0 1],'ellipseColor',[0 0 1],'ellipseWidth',2.5,'scale',0.05);
    else
        plotGeodesic(G{k},'Npts',101,'covFrequency',1001,'pathColor',[1 0 0],'startColor',[1 0 0],'endColor',[1 0 0],'ellipseColor',[1 0 0],'ellipseWidth',2.5,'scale',0.05);
    end
end

% Plot the original path
plotPath(pts, 'pathColor',[0.5 0.5 0.5], 'covFrequency', 1, 'startColor',[0.5 0.5 0.5], 'endColor', [0.5 0.5 0.5], 'scale', 0.05);
for k = 1:length(pts)
    plotPath(pts(k), 'pathColor',[0.5 0.5 0.5], 'covFrequency', 1, 'startColor',[0.5 0.5 0.5], 'endColor', [0.5 0.5 0.5], 'ellipseWidth',2.5, 'ellipseColor',[0.5 0.5 0.5],'scale',0.05);
end

% Plot the new path
for k = 1:length(pts)
    plotPath(new_pts(k), 'pathColor',[0 0 1], 'covFrequency', 1, 'startColor',[0 0 1], 'endColor', [0 0 1], 'ellipseWidth',2.5, 'ellipseColor',[0 0 1],'scale',0.05);
end

% Plot the true geodesic
plotGeodesic(trueG,'Npts',101,'covFrequency',20,'pathStyle','--','pathColor',[0 0 0],'startColor',[0 0 0],'endColor',[0 0 0],'ellipseColor',[0 0 0],'ellipseWidth',2.5,'scale',0.05);


% Append labels and legend
xlabel('\mu_1','FontSize',24);
ylabel('\mu_2','FontSize',24);
set(gca, 'FontSize', 24);

% Fix axis limits and thickness
xlim([-0.6 5.1]);
ylim([-0.1 5.6]);
axis square;
set(gca,'LineWidth',2);

% Return to the original directory
cd('Figures');


end

