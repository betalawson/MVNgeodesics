function FIGURE10_shootingPerformance(regenerate)
% This function creates the figure that displays the performance of
% shooting, and how it is improved by the use of residual vectors that are
% more representative of the MVN manifold's structure. If data is not
% present with the specified filename, or the optional "regenerate" flag is
% set to true, data is regenerated using the specifications below


%%% TESTING SPECIFICATIONS

% Specify the filename where the run data is stored (will be loaded from if
% it exists, and the regenerate flag was not manually set)
filename = 'DATA_shootingResiduals';

% Specify the random seed for reproducibility
rng(7);

% Specify list of dimensions to perform the test over
D_list = [2, 10];

% Number of geodesics per length to include in the test
N_geo = 250;

% Shooting algorithm parameters
max_vnorm = 0.5;
symKL_tol = 1e-7;

% Specify the list of d_{FR} values to test - these are the true geodesic
% distance (of randomly-generated test geodesics), and in general longer
% geodesics make it harder for approximations to be accurate
L_list = [1,2,3,4,5,6,7];

% Specify the list of approximations to use as residual vectors (repeats
% appear because initialisation status also varies)
method_list = {'euclid', 'eigen', 'projection', 'eigen', 'projection'};
init_flag = [false, false, false, true, true];

% Corresponding names for legend
legend_names = {'Euclidean', 'Axis-aligned', 'Embedded v', 'Axis-aligned (Initialised)', 'Embedded v (Initialised)'};


%%% Plotting Specification
         
% Select whether to plot runtimes or iteration counts (y limits will need to be adjusted)
show_runtimes = true;

% Figure size (as fraction of window)
figPos = [0.1, 0.1, 0.8, 0.725];
% Axis size (as fraction of figure, leaves room for legend)
axPos = [0.1, 0.19, 0.4, 0.75];
% Separation horizontally between axes on same plot
ax_sep = 0.075;

% y-limits are specified manually to allow better display of information
% (for example, to place more focus on general performance versus rare
% extremes)
ymin = log10(0.01);
ymax = log10(100);

% List the methods to actually show results for
plot_methods = [1,2,4,3,5];

% Specify the colours of the different methods
method_clrs = [0.4, 0.4, 0.4;
               0.4, 0.4, 1.0;
               0.7, 0.0, 0.7;
               0.7, 0.7, 1.0
               1.0, 0.7, 1.0  ];              
              
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
    N_methods = length(method_list);
    N_Ds = length(D_list);
    N_Ls = length(L_list);
    
    % Prepare storage objects
    iter_counts = cell(1,N_Ds);
    runtimes = cell(1,N_Ds);

    % Loop over each dimension requested
    for d = 1:N_Ds
        
        % Store current dimension
        D_here = D_list(d);
        
        % Prepare origin in this dimension
        O = struct('mu',zeros(D_here,1),'SIGMA',eye(D_here));
       
        % Loop over each length requested
        iter_counts{d} = zeros(N_Ls, N_methods, N_geo);
        runtimes{d} = zeros(N_Ls, N_methods, N_geo);
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
            iter_count_here = zeros(N_methods, N_geo);
            runtime_here = zeros(N_methods, N_geo);
            parfor k = 1:N_geo
                              
                % Run each method (residual choice + initialisation)
                for m = 1:N_methods
                    
                    % Call the requested shooting routine
                    if init_flag(m)
                        [~,diags] = onePointShooting(O,Ts{k},struct('symKL_tol',symKL_tol,'max_vnorm',max_vnorm,'approx_method',method_list{m},'Ginit',approximateGeodesic(O,Ts{k},method_list{m})));
                    else
                        [~,diags] = onePointShooting(O,Ts{k},struct('symKL_tol',symKL_tol,'max_vnorm',max_vnorm,'approx_method',method_list{m}));
                    end
                    
                    % If the method converged, store performance results
                    if diags.converged
                        iter_count_here(m,k) = diags.iterations;
                        runtime_here(m,k) = diags.runtime;
                    % Otherwise, mark failure with a result of "-1"
                    else
                        iter_count_here(m,k) = -1;
                        runtime_here(m,k) = -1;
                    end
                    
                end
                
            end
            fprintf('\nTime to calculate all methods for dimension n = %g and length dFR = %g was:\n',D_here,L);
            fprintf('Time: %g seconds\nTime per iteration: %g seconds\n\n',sum(runtime_here,'all'), sum(runtime_here,'all') / sum(iter_count_here,'all'));         
            
            % Store iteration and runtime data in the larger objects
            iter_counts{d}(l,:,:) = iter_count_here;
            runtimes{d}(l,:,:) = runtime_here;
            
        end
        
    end
        
    % Store the iteration count data and requested scenario parameters
    save(filename, 'iter_counts', 'runtimes', 'method_list', 'L_list', 'D_list');
   
% If data exists and not specifically requested, just load it in
else
    
    load(filename, 'iter_counts', 'runtimes', 'method_list', 'L_list', 'D_list');
    N_Ds = length(D_list);
    [N_Ls, N_methods, N_geo] = size(iter_counts{1});
    
end
    
    
%%% PLOTTING

% Initialise the figure
figure('units','normalized','OuterPosition',[0.2526, 0.1593, 0.6427, 0.6833]);
    
% Reset axis position
this_axPos = axPos;
    
% Create a subplot for each dimension
for d = 1:N_Ds
        
    % Increment the horizontal position of each plot after the first
    if d > 1
        this_axPos(1) = this_axPos(1) + this_axPos(3) + ax_sep;
    end
    axes('Position',this_axPos); 
    hold on;
        
    % Plot dummy data for each method on the last boxplot, regardless
    % of whether they appear on it - used for legend purposes later
    if d == N_Ds
        for k = 1:length(plot_methods)
            leg_dummies(k) = plot([NaN NaN], [NaN NaN],'-','Color', method_clrs(plot_methods(k),:), 'LineWidth', 2);
        end
    end
        
    % Remove any -1's from the boxplot, but count them for reporting
    % separately, these represent cases where shooting failed to converge
    % to requested tolerance
    if show_runtimes
        plot_data = runtimes{d}(:,plot_methods,:);
        ylabel_txt = 'log_{10} (Runtime)';
    else
        plot_data = iter_counts{d}(:,plot_methods,:);
        ylabel_txt = 'log_{10} (Iterations)';
    end
    fail_counts = sum(plot_data < 0, 3);
    plot_data(plot_data < 0) = NaN;
    
    % Create base boxplot
    h = boxplot2( log10(plot_data), 'Whisker', Inf );
        
    % Color the different approximations
    for k = 1:length(plot_methods)
        plot_clr = method_clrs( plot_methods(k),: );
        structfun(@(x) set( x(k,:), 'color', plot_clr, 'markeredgecolor', plot_clr, 'linewidth', 2.5), h);
    end
        
    % Plot clean up
    xticks(L_list);
    xlabel('$d_{FR}$','Interpreter','LaTeX');
    % Add label to leftmost axis only
    if d == 1
        ylabel(ylabel_txt);
    end
    set(gca,'FontSize',24,'LineWidth',2);
        
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
           
    % Add legend on the rightmost figure
    if d == N_Ds
        
        % Add spacing after legend names for spacing
        for k = 1:N_methods-1
            legend_names{k} = [legend_names{k},'{   }'];
        end
        
        leg_obj = legend(leg_dummies,legend_names(plot_methods),'Location','SouthOutside','Orientation','Horizontal');
    end
        
    % Set figure and axis sizes after all else done
    set(gcf,'Position',figPos);
    set(gca,'Position',this_axPos);
    
    % Add title showing the dimension
    title(['$n = ',num2str(D_list(d)),'$'],'Fontsize',24, 'Interpreter', 'LaTeX');
        
    % Set the legend to be at the centre if only plotting one subplot
    if N_Ds == 1
        leg_obj.Position(1) = (1 - leg_obj.Position(3))/2;
    end
    
    % Output the failure data to the user for separate reporting
    for l = 1:N_Ls
        fprintf('Failure counts for the methods with d = %g and l = %g were:',D_list(d),L_list(l));
        if all(fail_counts(l,:) == 0)
            fprintf('  None!\n\n');
        else
            fprintf('\n');
            for m = 1:N_methods
                if init_flag(m)
                    fprintf('%g failures out of %g for method %s with initialisation\n',fail_counts(l,m),N_geo,method_list{m})
                else
                    fprintf('%g failures out of %g for method %s with no initialisation\n',fail_counts(l,m),N_geo,method_list{m})
                end
            end
            fprintf('\n');
        end
    end
                
end