function FIGURE5b_approximationPerformance(regenerate)
%
% This function picks out random geodesics of different length, then
% formulates the approximate solutions using different methods and
% evaluates the average performance of each

% Specify the dimension to operate in (D-variate MVNs)
D = 10;

% Set the bounds for mean space
mu_bound = 5;
% Set the bounds for covariance eigenvalues (multiplicative factor)
SIGMA_bound = 5;

% Specify number of samples to use in calculation
N_samples = 5000;

% Define the colours
plot_colors = [ 
                0.40, 0.40, 0.40;        % Grey - Euclidean
                1.00, 0.00, 0.00;        % Red - Taylor
                0.00, 0.00, 1.00;        % Blue - Eigen
              ];  
          
% List method names
methods = {'euclid','taylor','eigen'};

% Corresponding method names for display on legend
legend_names = {'Euclidean', 'Taylor', 'Eigen'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume if no flag provided, not to regenerate data
if nargin < 1
    regenerate = false;
end

% Calculate numbers of options provided to abbreviate code
N_methods = length(methods);

% Regenerate data if requested/needed
if ~exist('data_approximationPerformance.mat','file') || regenerate

    % Move back into the main folder
    cd('..');

    % Initialise length storage
    problemL = zeros(N_samples,1);
    % Initialise the approximation performance storage
    KL_data = zeros(N_samples,N_methods);
    
    % Define the origin in the requested dimension
    O = struct('mu',zeros(D,1),'SIGMA',eye(D));

    % Loop over the requested number of sample shooting problems
    parfor m = 1:N_samples
        
        % Create a random point within the given bounds
        p = randomPoint(mu_bound, SIGMA_bound, D);
        
        % Find the true geodesic to this point (provides length)
        G = onePointShooting(O,p,struct('approx_method','eigen','max_vnorm',0.5));
        problemL(m) = sqrt(innerProduct(G.v,G.v,eye(D)));
        
        % Loop over each choice of approximation
        symKL_here = zeros(1,N_methods)
        for M = 1:N_methods
            
            % Form the approximate geodesic to this target point
            Ga = approximateGeodesic(O, p, methods{M});
            
            % Fire this geodesic and check its distance to the target
            pa = fireGeodesic(Ga,1);
            symKL_here(M) = symKL(p,pa);
            
        end
        
        % Store the data for this sample
        KL_data(m,:) = symKL_here;
        
    end
    
    % Move back to original folder
    cd('Figures');
        
    % Save the data
    save('data_approximationPerformance.mat','KL_data','problemL');
     
else
    
    % Load in the data
    load('data_approximationPerformance.mat','KL_data','problemL');
    
end


% Use the plotting subfunction below to generate a boxplot of update
% direction accuracy versus (binned) distance from current state to target
plotFunction(problemL, log10(KL_data), plot_colors);

% Add label too
ylabel('Approximation Error (log_{10} D_{sKL})','FontSize',24);

% Append legend (done here so extra lines have no legend entry)
legend(legend_names, 'Fontsize', 24);





function figobj = plotFunction(L, data, colors)
%
% This function creates a boxplot by binning L values and splitting up the
% data associated with different L values

% Define the bin spacing to use (affects plot ticks also so should be a
% nice number)
dL = 1;

% Find the edge of the bin that will contain all the data
Lmax = ceil(max(L) / dL) * dL;
% Split up the L values into bins with the given spacing and count how many
L_edges = 0:dL:Lmax;
N_bins = length(L_edges) - 1;
% Also store the middle of each bin for plotting purposes
L_bins = 0.5 * ( L_edges(1:end-1) + L_edges(2:end) );

% First loop over each bin and count how many L values fall in each
N_bin = zeros(1,N_bins);
for k = 1:N_bins
    N_bin(k) = sum( L > L_edges(k) & L <= L_edges(k+1) );
end

% Find the maximum number in a bin
N_max = max(N_bin);

% Initialise the data storage with room for the data in each bin (NaNs used
% for empty data)
plot_data = NaN(N_bins,size(data,2),N_max);

% Loop over each bin and store the data in the third dimension
for k = 1:N_bins
    
    if N_bin(k) > 0
        % Grab out which values fall within this bin
        inbin = L > L_edges(k) & L <= L_edges(k+1);
        % Grab out the data corresponding to the L values falling in this bin
        store_data = data(inbin, :);
        % Store this in the plot data 3-D array
        plot_data(k,:,1:N_bin(k)) = store_data';
    end
    
end

% Initialise figure
figobj = figure('units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Create base boxplot
h = boxplot2( plot_data, L_bins );

% Loop over each method's data and set its color
for k = 1:size(plot_data,2)
    structfun(@(x) set( x(k,:), 'color', colors(k,:), 'markeredgecolor', colors(k,:), 'linewidth', 2.5), h);
end

% Set outlier markers
set(h.out, 'marker', '.', 'markerSize', 20);

% Set line style
set([h.lwhis h.uwhis], 'linestyle', '-');

%%% Plot lines separating each boxplot

% Line locations are the bin edges
line_locs = L_edges;

% Plot these over a huge vertical range
hold on;
for k = 1:length(line_locs)
    plot([line_locs(k), line_locs(k)],[-50, 50],'k--','LineWidth',1.25);
end

%%% Add x ticks between the actual ticks that we use to give the number of
%%% samples in that bin

Nedges = length(L_edges);

xtick_vals = [];

xlim([0-0.025*Lmax, Lmax*1.025]);

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
    

% Set up labels and fontsize
xlabel('Manifold Distance (Binned)');
set(gca,'Fontsize', 24);

% Set up xlimits to cut off at the right spots
xlim([0-0.025*Lmax, Lmax*1.025]);