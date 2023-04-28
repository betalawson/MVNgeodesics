function FIGURE7ab_shootingVectorChoiceFull(regenerate)
%
% This function serves as a self-contained script to produce the figure
% shown in the paper comparing the performance of each different choice of
% "residual" vector for a whole run of the shooting algorithm.

% Set the bounds for mean space
mu_bound = 10;
% Set the bounds for covariance eigenvalues (multiplicative factor)
SIGMA_bound = 10;

% Specify the number of samples to use in boxplot generation (in total)
N_samples = 5000;

% Define the colours
plot_colors = [
    0.00, 0.00, 0.00;        % Black - Euclidean
    0.40, 0.40, 0.40;        % Grey - Taylor
    0.85, 0.00, 0.00;        % Red - Eigen
    0.00, 0.00, 0.85;        % Blue - Projection
    1.00, 0.55, 0.55;        % Light Red - Eigen (Initialised)
    0.55, 0.55, 1.00;        % Light Blue - Projection (Initialised)
    ];

% Specify the maximum norm of the residual vector to use
max_vnorm = 0.5;

% Define the dimension
D = 2;

% List method names
methods = {'euclid','taylor','eigen','projection','eigen','projection'};

% Specify which methods are initialised using that approximation
init_flag = [false, false, false, false, true, true];

% Corresponding method names for display on legend
legend_names = {'Euclidean', 'Small-x', 'Component-based', 'Projection', 'Component-based (Initialised)', 'Projection (Initialised)'};

% List methods to actually show results for
plot_methods = [ [1,2,3,4];
                 [3,5,4,6] ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume, if no flag provided, not to regenerate data
if nargin < 1
    regenerate = false;
end

% Calculate numbers of options provided to abbreviate code
N_methods = length(methods);

% Regenerate data if requested/needed
if ~exist('data_shootingVectorChoiceFullrun.mat','file') || regenerate
    
    % Move back into the main folder
    cd('..');
    
    % Initialise length storage
    problemL = zeros(N_samples,1);
    % Initialise the iteration storage
    iters_data = zeros(N_samples,N_methods);
    
    % Define the origin in the requested dimension
    O = struct('mu',zeros(D,1),'SIGMA',eye(D));
    
    % Loop over the requested number of sample shooting problems
    parfor m = 1:N_samples
        
        % Create a random point within the given bounds
        p1 = randomPoint(mu_bound, SIGMA_bound, D);
        
        % Initialise storage for shooting method results
        Gs_here = cell(1,N_methods);
        symKL_here = NaN(1,N_methods);
        iters_here = NaN(1,N_methods);
        converged_here = false(1,N_methods);
        
        % Loop over each choice of residual, run shooting for each
        for M = 1:length(methods)
            
            % Call the shooting routine using current choice of vector
            if init_flag(M)
                [Gs_here{M}, diags] = onePointShooting(O,p1,struct('max_vnorm',max_vnorm,'approx_method',methods{M},'Ginit',approximateGeodesic(O,p1,methods{M})));
            else
                [Gs_here{M}, diags] = onePointShooting(O,p1,struct('max_vnorm',max_vnorm,'approx_method',methods{M}));
            end
            
            % If method converged, store the performance results
            if diags.converged
                converged_here(M) = true;
                symKL_here(M) = diags.symKL_history(end);
                iters_here(M) = diags.iterations;
                
            % If method didn't converge, store a -1 to mark this failure
            else
                iters_here(M) = -1;
            end
        end
        
        % Use the geodesic that got closest to record length to the point
        if any(converged_here)
            [~,loc] = min(symKL_here);
            problemL(m) = Gs_here{loc}.L;
        end
        
        % Store the iteration results calculated here
        iters_data(m,:) = iters_here;
                
    end
    
    % Move back to original folder
    cd('Figures');
    
    % Save the data
    save('data_shootingVectorChoiceFullrun.mat','iters_data','problemL');
    
else
    
    % Load in the data
    load('data_shootingVectorChoiceFullrun.mat','iters_data','problemL');
    
end


%%% CREATE BINNED DATA

% Define the bin spacing to use (affects plot ticks also so should be a
% nice number)
dL = 1;

% Find the edge of the bin that will contain all the data
Lmax = ceil(max(problemL) / dL) * dL;
% Split up the L values into bins with the given spacing and count how many
L_edges = 0:dL:Lmax;
N_bins = length(L_edges) - 1;
% Also store the middle of each bin for plotting purposes
L_bins = 0.5 * ( L_edges(1:end-1) + L_edges(2:end) );

% First loop over each bin and count how many L values fall in each
N_bin = zeros(1,N_bins);
for k = 1:N_bins
    N_bin(k) = sum( problemL > L_edges(k) & problemL <= L_edges(k+1) );
end

% Find the maximum number in a bin
N_max = max(N_bin);

% Initialise the data storage with room for the data in each bin (NaNs used
% for empty data)
full_plot_data = NaN(N_bins,size(iters_data,2),N_max);

% Loop over each bin and store the data in the third dimension
for k = 1:N_bins
    
    if N_bin(k) > 0
        % Grab out which values fall within this bin
        inbin = problemL > L_edges(k) & problemL <= L_edges(k+1);
        % Grab out the data corresponding to the L values falling in this bin
        store_data = iters_data(inbin,:);
        % Store this in the plot data 3-D array
        full_plot_data(k,:,1:N_bin(k)) = store_data';
    end
    
end


%%% LOOP OVER THE DIFFERENT PLOTS REQUESTED
for p = 1:size(plot_methods,1)

    %%% CUT DOWN THE PLOT DATA TO THE METHODS TO BE PLOTTED
    plot_data = full_plot_data(:,plot_methods(p,:),:);

    %%% FIND MAXIMUM ITERATION COUNT ACROSS RUNS (FOR AXIS LIMITS PURPOSES)
    max_iters = max(plot_data(:));

    %%% CREATE ACCOMPANYING FAILURE DATA IF NEEDED

    if any(plot_data(:) == -1)
        failures_present = true;
    else
        failures_present = false;
    end

    if failures_present
    
        % Set the fail marker to a value that's 1.2x the maximum
        fail_val = 1.2 * max_iters;
    
        % Create a new array in which all values are 100 x fail_val, initially
        fail_data = 100*fail_val*ones(size(plot_data));
        % Loop over each bit of the fail data so that if the method failed in this
        % bin, one value of fail_val (instead of 100xfail_val) is recorded
        for k = 1:size(plot_data,1)
            for j = 1:size(plot_methods,2)
                if any( squeeze(plot_data(k,j,:)) == -1)
                    fail_data(k,j,1) = fail_val;
                end
            end
        end
    
        % Now remove all the -1's that indicate a failure to converge from the
        % original plot data
        plot_data(plot_data == -1) = NaN;
    
    end

    %%% PLOT CODE

    %%% Initialise figure
    figobj = figure('units', 'normalized', 'OuterPosition', [0 0 1 1]);


    %%% Create base boxplot
    h = boxplot2( plot_data, L_bins );

    % Loop over each method's data and set its color
    for k = 1:size(plot_data,2)
        structfun(@(x) set( x(k,:), 'color', plot_colors(plot_methods(p,k),:), 'markeredgecolor', plot_colors(plot_methods(p,k),:), 'linewidth', 2.5), h);
    end

    % Set outlier markers
    set(h.out, 'marker', '.', 'markerSize', 20);

    % Set line style
    set([h.lwhis h.uwhis], 'linestyle', '-');


    %%% Create failure 'boxplot' (just shows x's at the failure value)
    if failures_present
        h2 = boxplot2( fail_data, L_bins );
    
        % Loop over each method's data and set its color
        for k = 1:size(fail_data,2)
            structfun(@(x) set( x(k,:), 'color', plot_colors(plot_methods(p,k),:), 'markeredgecolor', plot_colors(plot_methods(p,k),:), 'linewidth', 2.5), h2);
        end
    
        % Set outlier markers
        set(h2.out, 'marker', 'x', 'markerSize', 30);
    end


    %%% Plot lines separating each boxplot

    % Line locations are the bin edges
    line_locs = L_edges;

    % Plot these over a huge vertical range
    hold on;
    for k = 1:length(line_locs)
        plot([line_locs(k), line_locs(k)],[-50, 10*max_iters],'k--','LineWidth',1.25);
    end


    %%% Extra plot tinkering

    % Set up labels and legend and fontsize
    xlabel('Geodesic Length (Binned)');
    ylabel('Number of Iterations');
    legend(legend_names(plot_methods(p,:)), 'Fontsize', 24);
    set(gca,'Fontsize', 24);

    % Set up xlimits to cut off at the right spots
    xlim([0-0.025*Lmax, Lmax*1.025]);

    % Set up ylimits from 0 up to a bit above the maximum value (and the
    % failure value if being used)
    ylim([0 1.25*max_iters]);


    %%% Add x ticks between the actual ticks that we use to give the number of
    %%% samples in that bin

    Nedges = length(L_edges);

    xtick_vals = [];

    % Loop over all but the last edge, adding a tick for both the edge and the
    % centre value
    for k = 1:Nedges-1
    
        xtick_vals(end+1) = L_edges(k);
        xtick_vals(end+1) = 0.5 * (L_edges(k) + L_edges(k+1));
    
    end

    % Store also the last edge value
    xtick_vals(end+1) = L_edges(end);

    % Put these ticks on the plot
    xticks(xtick_vals);

    % Grab out the x tick labels from the plot
    xtick_labels = xticklabels();

    % Change each bin's tick to the number of samples used to form that bin
    for k = 1:Nedges-1
        xtick_labels{2*k} = ['\fontsize{18}\color{blue} (',num2str(N_bin(k)),')'];
    end

    % Put these modified ticks back on the figure
    xticklabels(xtick_labels);

    % Make sure the x labels don't get rotated
    xtickangle(0);

    
    
    %%% If needed, set up y-ticks so the failure value is displayed in text
    
    if failures_present
        % Get out the current ticks used
        yvals = get(gca,'ytick');
        % Keep only these up to one after the maximum value plotted
        yticks_use = yvals(yvals <= max_iters);
        yticks_use(end+1) = yvals(length(yticks_use)+1);
        % Append a final tick for the failure value
        yticks_use(end+1) = fail_val;
        % Put these on the plot
        set(gca,'ytick',yticks_use);
        % Now grab these labels out
        ytick_labels = get(gca,'yticklabels');
        % Set the last one to text
        ytick_labels{end} = 'FAIL';
        % Put these labels on the plot
        set(gca,'yticklabels',ytick_labels);
    end
    
end