function FIGURE11_pathRefinement
% This function visualises a single step of the path refinement algorithm,
% and in the paper is run for N=5 and N=9 to also demonstrate how larger
% numbers of path points result in slower convergence of the algorithm

% Specify number of points to use in the example
N = 11;

% Define the problem to use as an example
P1 = struct('mu',[0;0],'SIGMA',[10,0;0,1]);
P2 = struct('mu',[5;5],'SIGMA',[1,0;0,10]);

% Specify the colours for plotting
init_clr = [0.5, 0.5, 0.5];         % Initial path
geo_clr = [0, 0, 0];                % True geodesic between the points
inter_clr = [1, 0, 0];              % Intermediate geodesics between even-index points
iter_clr = [0, 0, 1];               % Colour of new path after one iteration

% Specify axis limits for the plot
xmin = -0.6;
xmax = 5.25;
ymin = -0.25;
ymax = 5.6;

% Specify the text locations on figure (for length scores)
text_x = -0.1;
text_y = 5.3;
% Specify the text offset between length scores
dy = 0.4;

% Specify a custom scaling for covariances (makes figure more legible)
cov_scale = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Form the Euclidean path between these
pts = closedFormPath(P1,P2,N,'euclid');
new_pts = pts;

% Calculate the length of the original path (recreate it with many more
% points for numerical length calculation)
Epath = closedFormPath(P1,P2,10001,'euclid');
E_L = pathLength(Epath);


% Join up the even index points with geodesics
for k = 3:2:N-1
    G{k-1} = onePointShooting(pts{k-1},pts{k+1});
    new_pts{k} = fireGeodesic(G{k-1},0.5);
end

% Join up the new odd-index points with geodesics
for k = 2:2:N-1
    G{k-1} = onePointShooting(new_pts{k-1},new_pts{k+1});
    new_pts{k} = fireGeodesic(G{k-1},0.5);
end

% Count up the total length of the updated path (it is constructed from
% geodesic segments between the odd-index points so we can simply add their
% stored lengths)
PR_L = 0;
for k = 1:2:N-1
    PR_L = PR_L + G{k}.L;
end



% Also find the true geodesic and read out its length
trueG = onePointShooting(P1,P2);
G_L = trueG.L;

% Prepare figure
figure; hold on;

% Plot the geodesics used
for k = 1:N-2
    if mod(k,2) == 1
        plotGeodesic(G{k}, 'pathWidth', 3,'Npts',101,'covFrequency',50,'pathColor',iter_clr,'startColor',iter_clr,'endColor',iter_clr,'ellipseColor',iter_clr,'ellipseWidth',3,'scale',cov_scale);
    else
        plotGeodesic(G{k}, 'pathWidth', 3,'Npts',101,'covFrequency',50,'pathColor',inter_clr,'startColor',inter_clr,'endColor',inter_clr,'ellipseColor',inter_clr,'ellipseWidth',3,'scale',cov_scale);
    end
end

% Plot the original path
plotPath(pts, 'pathWidth', 3, 'pathColor',init_clr, 'covFrequency', 1, 'startColor',init_clr, 'endColor', init_clr, 'scale', cov_scale);
for k = 1:length(pts)
    plotPath(pts(k), 'pathWidth', 3, 'pathColor',init_clr, 'covFrequency', 1, 'startColor',init_clr, 'endColor', init_clr, 'ellipseWidth',3, 'ellipseColor',[0.5 0.5 0.5],'scale',cov_scale);
end

% Plot the new path
%for k = 1:length(pts)
%    plotPath(new_pts(k), 'pathWidth', 3, 'pathColor',iter_clr, 'covFrequency', 1, 'startColor',iter_clr, 'endColor', iter_clr, 'ellipseWidth',3, 'ellipseColor',[0 0 1],'scale',cov_scale);
%end

% Plot the true geodesic
plotGeodesic(trueG, 'pathWidth', 3,'Npts',1+100*N,'covFrequency',100,'pathStyle','--','pathColor',geo_clr,'startColor',geo_clr,'endColor',geo_clr,'ellipseColor',geo_clr,'ellipseWidth',3,'scale',cov_scale);


% Append labels and legend
xlabel('\mu_1','FontSize',24);
ylabel('\mu_2','FontSize',24);
set(gca, 'FontSize', 24);

% Fix axis limits and thickness
xlim([xmin xmax]);
ylim([ymin ymax]);
axis square;
set(gca,'LineWidth',2);

% Set figure shape
set(gcf,'Units','Normalized');
set(gcf,'OuterPosition',[0.35, 0.14, 0.43, 0.74]);
set(gcf,'Position',[0.35, 0.25, 0.42, 0.66]);

% Append length scores for the different paths
plot_lengths = [E_L, PR_L, G_L];
plot_clrs = [init_clr; iter_clr; geo_clr];
for p = 1:length(plot_lengths)
        
    % Text label
    text( text_x, text_y - (p-1)*dy, ['L = ',num2str(plot_lengths(p),4)], 'Color', plot_clrs(p,:), 'FontSize', 24);

end

