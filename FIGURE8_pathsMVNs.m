function FIGURE8_pathsMVNs(regenerate)
% This function creates the content of Figure 8 and Figure S1, showing some
% example paths in MVN space (Figure 8a), and then the performance of those
% different paths (in terms of information efficiency) across many 
% randomly-generated examples for two-dimensional MVNs (Figure 8b) and for
% ten-dimensional MVNs (Figure S1)
% 
% If performance data is not present at the requested filename, or the
% optional flag to regenerate data is set, then this data will be
% generated.


%%% Shared specifications

% Specify the list of paths to try
path_list = {'anneal', 'moment', 'wasserstein', 'entropy', 'projection'};
% Specify the names as they will appear on the legend (matching above list of paths, but with the geodesic path attached to the end)
legend_names = {'Annealing', 'Moment-averaged', 'd_{\epsilonW} (\epsilon = 0)', 'd_{\epsilonW} (\epsilon = 10)', 'Embedding'};


% Specify number of points along paths (for numerical length calculation)
N_pathpts = 10001;

% Specify the entropy parameter to use for entropy-reguarlised Wasserstein
epsilon = 10;

% Specify the colours of the different plots
plot_clrs = [1.0, 0.4, 0.4;
             0.4, 0.4, 1.0;
             0.7, 0.0, 0.7;
             1.0, 0.4, 1.0;
             0.8, 0.6, 0.2  ];
             
% Figure size (as fraction of window)
figPos = [0.1, 0.1, 0.8, 0.725];
% Axis size (as fraction of figure, leaves room for legend)
axPos = [0.075, 0.20, 0.4, 0.75];
% Separation between axes
ax_sep = 0.1125;


%%% Path example specification

% Specify the list of target points to find geodesics to
T_examples{1} = struct('mu',[4;3],'SIGMA',[1.3, -0.3;-0.3, 0.3]);

% Specify how many points along paths to place between covariance ellipses
N_pathpts_ellipse = 2000;

% Axis limits
ex_xlim = [-0.25 4.25];
ex_ylim = [-0.75 3.75];

% Specify the text locations on figure (provide as vectors if there are
% multiple examples)
text_x = -0.1;
text_y = 3.3;
% Specify the text offset
dy = 0.2;


%%% Performance testing specifications

% Specify the filename where the run data is stored (will be loaded from if
% it exists, and the regenerate flag was not manually set)
filename = 'DATA_pathTesting';

% Specify the random seed for reproducibility
rng(7);

% Specify the list of dimensions to calculate over
D_list = [2,10];

% Number of geodesics to include in the test
N_geo = 5000;

% Specify the list of d_{FR} values to test - these are the true geodesic
% distance (of randomly-generated test geodeiscs), and in general longer
% geodesics make it harder for approximations to be accurate
L_list = [1,2,3,4,5,6,7];

% For plotting, specify the y-axis limits (used to keep the focus of the
% boxplots on the central data and remove some plotting of extreme results)
ymin = -6;
ymax = 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% PERFORMANCE DATA LOADING/GENERATION

% Assume data is not to be regenerated unless specifically told so
if nargin < 1
    regenerate = false;
end

% Append .mat to the filename if it wasn't given already
if ~strcmp(filename(end-3:end),'.mat')
    filename = [filename, '.mat'];
end

% Check if data is present with specified filename, otherwise generate it
if regenerate || ~exist(filename, 'file')

    % Count number of paths and dimensions requested to calculate for
    N_paths = length(path_list);
    N_Ds = length(D_list);
    N_Ls = length(L_list);
    
    % Prepare storage object
    rel_lengths = cell(1,N_Ds);
    
    % Loop over each dimension requested
    for d = 1:length(D_list)
        
        % Store current dimension
        D_here = D_list(d);
        
        % Prepare origin in this dimension
        O = struct('mu',zeros(D_here,1),'SIGMA',eye(D_here));
       
        % Loop over each length requested
        dFR_err{d} = zeros(N_Ls, N_paths, N_geo);
        for l = 1:N_Ls
            
            % Read out the current length
            L = L_list(l);
            
            % Generate requested number of target points (outside of
            % parallel loop for easier reproducibility)
            Ts = cell(1,N_geo);
            for k = 1:N_geo
            
                % Generate a random tangent space vector (initial geodesic velocity)
                v = randomUnitVector(D_here);           
            
                % Create geodesic with this velocity
                G = struct('P', eye(D_here), 'r', zeros(D_here,1), 'v', v);
            
                % Fire this geodesic to current length to get target point
                Ts{k} = fireGeodesic(G, L);
                
            end
        
            % Repeat across all the target points, generating the paths to each
            rel_lengths_here = zeros(N_paths, N_geo);
            parfor k = 1:N_geo
    
                % Now calculate the lengths of the listed set of paths
                for p = 1:N_paths
                    % Generate path and find its length
                    path = closedFormPath(O, Ts{k}, N_pathpts, path_list{p}, epsilon);
                    pathL = pathLength(path);
                    % Determine its relative "error" (wasted length) and store
                    rel_lengths_here(p,k) = (pathL - L) / L;
                end
    
            end
            
            % Store path length data in the larger object
            rel_lengths{d}(l,:,:) = rel_lengths_here;
            
        end
    
    end
        
    % Store the results and the requested scenario parameters
    save(filename, 'rel_lengths', 'path_list', 'epsilon', 'L_list', 'D_list');   
    

% If data exists and not specifically requested, just load it in
else
    
    load(filename, 'rel_lengths', 'path_list', 'epsilon', 'L_list', 'D_list');
    N_Ds = length(D_list);
    [N_Ls, N_paths, N_geo] = size(rel_lengths{1});
    
end


%%% PLOT INITIALISATION

% Initialise the figure
first_figobj = figure('units','normalized','OuterPosition',[0.2526, 0.1593, 0.6427, 0.6833]);
% Prepare axes for showing the path example
axes('Position',axPos); 
hold on;


%%% PATH EXAMPLE GENERATION

% Define the origin for bivariate MVNs
O = struct('mu',zeros(2,1),'SIGMA',eye(2));

% Loop over each target point, plotting the geodesic path along with the other paths
for k = 1:length(T_examples)
    
    % Grab out current target point
    T = T_examples{k};
    
    % Find the geodesic path to the target
    G = onePointShooting(O, T);
        
    % Generate the requested closed-form paths to the target and calculate
    % their length according to the metric (numerically)
    cPaths = cell(1,N_paths);
    cPath_Ls = zeros(1,N_paths);
    for p = 1:N_paths
        cPaths{p} = closedFormPath(O, T, N_pathpts, path_list{p}, epsilon);
        cPath_Ls(p) = pathLength(cPaths{p});
    end

    % Plot the closed-form paths
    for p = 1:N_paths
        
        % Path itself
        plotPath(cPaths{p}, 'pathColor', plot_clrs(p,:), 'ellipseColor', plot_clrs(p,:), 'ellipseWidth', 2.5, 'pathWidth', 5, 'covFrequency', N_pathpts_ellipse, 'Npts', N_pathpts);
        % Text label
        text( text_x(k), text_y(k) - (p-1)*dy, ['L = ',num2str(cPath_Ls(p),4)], 'Color', plot_clrs(p,:), 'FontSize', 24);
        
    end
       
    % Plot the geodesic on top (use dashed line)
    plotGeodesic(G, 'pathStyle', ':', 'pathColor', [0, 0, 0], 'ellipseColor', [0 0 0], 'ellipseWidth', 2.5, 'pathWidth', 5, 'covFrequency', N_pathpts_ellipse, 'Npts', N_pathpts, 'StartColor', [0 0 0]);
    
    % Append just the geodesic's length
    text( text_x(k), text_y(k) - p*dy, ['L = ',num2str(G.L,4)], 'Color', [0 0 0], 'FontSize', 24, 'HorizontalAlignment', 'Left');
        
end

% Fix axis limits and thickness
xlabel('\mu_1','FontSize',24);
ylabel('\mu_2','FontSize',24);
axis equal;
xlim(ex_xlim); 
ylim(ex_ylim);
set(gca,'YTickMode','manual');
% Manually set ticks to be in spacings of 0.5
xticks( floor(2*ex_xlim(1))/2 : 0.5 : ceil(2*ex_xlim(2))/2 );
yticks( floor(2*ex_ylim(1))/2 : 0.5 : ceil(2*ex_ylim(2))/2 );
set(gca,'LineWidth',2, 'FontSize', 24);



%%% PLOT PERFORMANCE DATA

% Loop over each dimension, plotting a separate figure for each
for d = 1:N_Ds
    
    % For first plot, put it on a second axis on the first figure
    if d == 1
        axPos = axPos;    
        axPos(1) = axPos(1) + axPos(3) + ax_sep;
        axes('Position',axPos);
        hold on
    
    % Otherwise, just create a new figure and use that same axis position
    else
        figure('units','normalized','OuterPosition',[0.2526, 0.1593, 0.6427, 0.6833]);
        axes('Position',axPos); 
        hold on;
    end
        
    % Eliminate any numerical error data (lengths so close to geodesic that
    % numerical error registers them as shorter, which is impossible)
    rel_lengths{d}(rel_lengths{d} < 0) = NaN;
        
    % Create base boxplot
    h = boxplot2( log10(rel_lengths{d}), 'Whisker', Inf );
    
    % Color the different path types
    for k = 1:N_paths
        structfun(@(x) set( x(k,:), 'color', plot_clrs(k,:), 'markeredgecolor', plot_clrs(k,:), 'linewidth', 2.5), h);
    end
    
    % Plot clean up
    xticks(L_list);
    xlabel('$d_{FR}$','Interpreter','LaTeX');
    ylabel('$\log_{10} \Bigl[(L - d_{FR})/d_{FR}\Bigr]$','Interpreter','LaTeX');
    set(gca,'FontSize',24, 'LineWidth', 2);
    
    % Axis limits and ticks
    ylim([ymin - 0.02*(ymax-ymin), ymax + 0.02*(ymax-ymin)]);
    xlim([0.5,N_Ls+0.5]);
    xticks(1:N_Ls);
    for l = 1:N_Ls
        L_labels{l} = num2str(L_list(l));
    end
    xticklabels(L_labels);
    
    % Add vertical lines to separate categories
    for k = 1:N_Ls-1
        plot( k+0.5*ones(1,2), [ymin - 0.05*(ymax-ymin), ymax + 0.05*(ymax-ymin)], 'k--', 'LineWidth', 2);
    end
    
    % Add spacing after legend names for spacing
    for k = 1:N_paths-1
        legend_names{k} = [legend_names{k},'{   }'];
    end
    
    % Add legend to each figure
    leg_obj = legend(legend_names,'Location','SouthOutside','Orientation','Horizontal');
    
    % Set figure size after all else done
    set(gcf,'Position',figPos);
    set(gca,'Position',axPos);
   
end

% Return to first figure and add the a) and b) labels
figure(first_figobj);
annotation('textbox', [0.01, 1, 0, 0], 'string', '\bf a) \rm', 'FontSize', 20)
annotation('textbox', [0.52, 1, 0, 0], 'string', '\bf b) \rm', 'FontSize', 20)