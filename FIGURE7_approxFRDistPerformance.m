function FIGURE7_approxFRDistPerformance(regenerate)
% This function displays the performance of different types of
% approximations for the Fisher-Rao distance, in a box plot that shows
% results across a large number of runs. If data is not present at the
% requested filename, or the optional flag to regenerate data is set, then
% this data will be generated.

%%% Testing specifications

% Specify the filename where the run data is stored (will be loaded from if
% it exists, and the regenerate flag was not manually set)
filename = 'DATA_distTesting';

% Specify the random seed for reproducibility
rng(7);

% Specify the list of dimensions to calculate over
D_list = [2,10];

% Number of geodesics per length to include in the test
N_geo = 5000;

% Specify the list of d_{FR} values to test - these are the true geodesic
% distance (of randomly-generated test geodeiscs), and in general longer
% geodesics make it harder for approximations to be accurate
L_list = [1,2,3,4,5,6,7];

% Specify the list of approximations to try
approx_list = {'taylor', 'eigen', 'projection', 'path'};
% Specify the names as they will appear on the legend (matching above list of paths, but with the geodesic path attached to the end)
legend_names = {'Small-u', 'Axis-aligned', 'Embedded v', 'Embedded path'};
                                           
% Specify whether to plot log10( abs(error)) for legible magnitudes, or
% just raw error data (useful for checking overestimation/underestimation
plot_logs = true;

%%% Plotting Specification

% Specify the colours of the different plots
plot_clrs = [1.0, 0.4, 0.4;
             0.4, 0.4, 1.0;
             0.7, 0.0, 0.7;
             1.0, 0.4, 1.0;
             0.8, 0.6, 0.2  ];
         
% For plotting, specify the y-axis limits (used to keep the focus of the
% boxplots on the central data and remove some plotting of extreme results)
%   ADJUST THESE IF SWITCHING OFF OF 'true' FOR plot_logs
ymin = -6;
ymax = 0.3;
         
% Figure size (as fraction of window)
figPos = [0.1, 0.1, 0.8, 0.725];
% Axis size (as fraction of figure, leaves room for legend)
axPos = [0.1, 0.20, 0.4, 0.75];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume data is not to be regenerated unless specifically told so
if nargin < 1
    regenerate = false;
end

% Append .mat to the filename if it wasn't given already
if ~strcmp(filename(end-3:end),'.mat')
    filename = [filename, '.mat'];
end


%%% DATA LOADING/GENERATION

% Check if data is present with specified filename, otherwise generate it
if regenerate || ~exist(filename, 'file')

    % Count number of approximations, dimensions and lengths to test
    N_approx = length(approx_list);
    N_Ds = length(D_list);
    N_Ls = length(L_list);
    
    % Prepare storage object
    dFR_err = cell(1,N_Ds);
    
    % Loop over each dimension requested
    for d = 1:N_Ds
        
        % Store current dimension
        D_here = D_list(d);
        
        % Prepare origin in this dimension
        O = struct('mu',zeros(D_here,1),'SIGMA',eye(D_here));
       
        % Loop over each length requested
        dFR_err{d} = zeros(N_Ls, N_approx, N_geo);
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
            
            % Calculate approximations for the target points in parallel
            dFR_err_here = zeros(N_approx, N_geo);
            parfor k = 1:N_geo
                for a = 1:N_approx
                    Ga = approximateGeodesic(O,Ts{k},approx_list{a});
                    dFR_err_here(a,k) = (Ga.L - L)/L;
                end
            end
    
            % Store approximation error data in the larger object
            dFR_err{d}(l,:,:) = dFR_err_here;
            
        end
    
    end
        
    % Store the results and the requested scenario parameters
    save(filename, 'dFR_err', 'approx_list', 'L_list', 'D_list');   
    

% If data exists and not specifically requested, just load it in
else
    
    load(filename, 'dFR_err', 'approx_list', 'L_list', 'D_list');
    N_Ds = length(D_list);
    [N_Ls, N_approx, N_geo] = size(dFR_err{1});
    
end


%%% PLOTTING

% Initialise the figure
figure('units','normalized','OuterPosition',[0.2526, 0.1593, 0.6427, 0.6833]);

% Plot a figure for each separate dimension
for d = 1:N_Ds

    % Increment the horizontal position of each subsequent plot
    if d > 1
        axPos(1) = axPos(1) + axPos(3) + 0.075;
    end
    axes('Position',axPos); 
    hold on;
   
    % Create base boxplot
    if plot_logs
        h = boxplot2( log10(abs(dFR_err{d})), 'Whisker', Inf );
    else
        h = boxplot2( dFR_err{d}, 'Whisker', Inf );
    end
    
    % Color the different approximations
    for k = 1:N_approx
        structfun(@(x) set( x(k,:), 'color', plot_clrs(k,:), 'markeredgecolor', plot_clrs(k,:), 'linewidth', 2.5), h);
    end
    
    % Plot clean up
    xticks(L_list);
    xlabel('$d_{FR}$','Interpreter','LaTeX');
    if d == 1
        ylabel('$\log_{10} \Bigl[(\hat{d}_{FR} - d_{FR})/d_{FR}\Bigr]$','Interpreter','LaTeX');
    end
    set(gca,'FontSize',24, 'LineWidth', 2);
    title(['n = ',num2str(D_list(d))],'FontSize',24);
    
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
    for k = 1:N_approx-1
        legend_names{k} = [legend_names{k},'{   }'];
    end
    
    % Add legend on the rightmost figure
    if d == N_Ds
        leg_obj = legend(legend_names,'Location','SouthOutside','Orientation','Horizontal');
    end
    
    % Set figure and axis sizes after all else done
    set(gcf,'Position',figPos);
    set(gca,'Position',axPos);
    
    % Report worst performance of the two best-performing methods
    fprintf('\n');
    for l = 1:N_Ls
        fprintf('Approximations for dimension d = %g and points d_FR = %g away had at most:\n', D_list(d), L_list(l));
        for m = 1:N_approx
            fprintf('%g%% error for method ''%s''\n',max(abs(dFR_err{d}(l,m,:)))*100,approx_list{m});
        end
        fprintf('\n');
    end
        
end

