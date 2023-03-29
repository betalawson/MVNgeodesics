function FIGURE7_shootingVectorChoice(regenerate)
%
% This function serves as a self-contained script to produce the figure
% shown in the paper comparing the performance of each different choice of
% "residual" vector for a single step of shooting across many different
% randomly generated scenarios. The user may specify whether to regenerate
% data, or plot using existing data (if it exists)

% Set the bounds for mean space
mu_bound = 10;
% Set the bounds for covariance eigenvalues (multiplicative factor)
SIGMA_bound = 10;

% Define the dimension
D = 2;

% Specify the number of samples to use in boxplot generation (in total)
N_samples = 10;

% Define the colours
plot_colors = [
    0.40, 0.40, 0.40;        % Grey - Euclidean
    1.00, 0.00, 0.00;        % Red - Taylor
    0.00, 0.00, 1.00;        % Blue - Eigen
    0.90, 0.00, 0.90;        % Purple - Transport
    ];

% Specify the maximum norm of the residual vector before velocity updates
% are scaled down - while finding true geodesics
max_vnorm_safe = 0.05;

% Specify the same maximum norm values to use in simulated shooting updates
max_vnorm = [Inf, 2, 0.5];

% List the different choices of vector to define between points
methods = {'euclid','taylor','eigen','transport'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume if no flag provided, not to regenerate data
if nargin < 1
    regenerate = false;
end

% Calculate numbers of options provided to abbreviate code
N_methods = length(methods);
N_vnorms = length(max_vnorm);

% Regenerate data if requested/needed
if ~exist('data_shootingVectorChoice.mat','file') || regenerate
    
    % Move back into the main folder
    cd('..');
    
    % Initialise length storage
    problemL = zeros(N_samples,1);
    % Initialise the angle and performance storage
    theta_data = zeros(N_samples,N_methods);
    KLimprovement_data = zeros(N_samples,N_methods,N_vnorms);
    
    % Define the origin in this dimension (starting point for all
    % geodesics fired)
    O = struct('mu',zeros(D,1),'SIGMA',eye(D));

    parfor m = 1:N_samples
        
        % Create a random current point for simulated shooting iteration
        p1 = randomPoint(mu_bound, SIGMA_bound, D);
        
        % Find the geodesic to this point from the origin
        [G1, diags1] = onePointShooting(O,p1,struct('amx_vnorm',max_vnorm_safe));
               
        % Create a random target point for the simulated shooting problem
        p2 = randomPoint(mu_bound, SIGMA_bound, D);

        % Find the geodesic connecting the first and second points
        % (provides the manifold distance between the two)
        [G2, diags2] = onePointShooting(p1,p2,struct('max_vnorm',max_vnorm_safe));
              
        % Formulate the Euclidean connection between the two points also
        R = struct('mu',p2.mu - p1.mu,'SIGMA',p2.SIGMA - p1.SIGMA);
        
        % Find the velocity update that would get to the target
        [Gtrue, diags3] = onePointShooting(O,p2,struct('max_vnorm',max_vnorm_safe));
        dvtrue = struct('mu',Gtrue.v.mu - G1.v.mu,'SIGMA',Gtrue.v.SIGMA - G1.v.SIGMA);
        
        % If all geodesics were found successfully, continue to process
        % this point
        if all([diags1.converged, diags2.converged, diags3.converged])
            
            % Store the length from current point to target point
            problemL(m) = G2.L;
            
            % Initialise storage for this simulated shooting problem
            this_theta = zeros(N_methods,1);
            this_KLimprovement = zeros(N_methods,N_vnorms);
            
            % Loop over each method (choice of connection vector)
            for M = 1:N_methods
                
                % Find the connection vector generated using this method
                Ga = approximateGeodesic(p1, p2, methods{M});
                va = Ga.v;
                
                % Find the norm of the connection vector
                norm_va = innerProduct(va, va, eye(D));

                % Vector va is for the space where origin is at p1, so
                % transform back to space where origin is at O
                va = affineTransform(va,Ga.P,zeros(D,1));
                
                % Parallel transport this vector backwards along G1 to get
                % a velocity correction (hence termed dva) at the origin
                dva = parallelTransport(va,G1,[1 0]);
                
                % Find the angle between the correction and the correction
                % that would produce the correct geodesic
                this_theta(M) = angleBetween(dva,dvtrue,eye(D));
                
                % Find the Jacobi field for this update applied to G1
                J = jacobiField(G1,dva,1);
            
                % Calculate step length via projection onto the Jacobi field
                invSIGMA = p1.SIGMA \ eye(D);
                s_base = innerProduct(R,J,invSIGMA) / innerProduct(J,J,invSIGMA);
            
                % Loop over the different maximum norm choices, taking the
                % appropriate step and evaluating performance
                for S = 1:length(max_vnorm)
                               
                    % Shrink stepsize if the connection vector norm exceeds maximum
                    s = min(1, max_vnorm(S) / norm_va) * s_base;
                
                    % Find the updated velocity with this stepsize
                    v_new = struct('mu',G1.v.mu + s*dva.mu,'SIGMA',G1.v.SIGMA + s*dva.SIGMA);
                    
                    % Create and fire the geodesic with this velocity
                    Gnew = struct('v',v_new,'P',eye(D),'r',zeros(D,1));
                    pnew = fireGeodesic(Gnew,1);
                    
                    % Calculate the improvement in symKL attained
                    this_KLimprovement(M,S) = log10( symKL(p1,p2) / symKL(pnew,p2) );
                    
                    % If the calculated symKL improvement is invalid,
                    % assign it dummy data so it doesn't show on boxplot
                    if ~isreal(this_KLimprovement(M,S))
                        this_KLimprovement(M,S) = NaN;
                    end
                    
                end
                                
            end
            
        else
                
            % Failed scenarios have dummy data assigned
            this_theta = nan(1,N_methods);
            this_KLimprovement = nan(N_methods,N_vnorms);
            problemL(m) = NaN;
        
        end
    
        % Store the data for this shooting scenario
        theta_data(m,:) = this_theta;
        KLimprovement_data(m,:,:) = this_KLimprovement;
        
    end
    
    % Save the data
    save('data_shootingVectorChoice.mat','problemL','theta_data','KLimprovement_data','-v7.3');
    
else
    
    % Load in the data
    load('data_shootingVectorChoice.mat','problemL','theta_data','KLimprovement_data');
    
end
    

%%% THETA PLOT

% Use the plotting subfunction below to generate a boxplot of update
% direction accuracy versus (binned) distance from current state to target
plotFunction(problemL, theta_data, plot_colors);

% y-axis is the angle, add ticks
ylim([0, pi]);
yticks([0 pi/8 pi/4 3*pi/8 pi/2, 5*pi/8, 3*pi/4, 7*pi/8, pi]);
yticklabels({'0','\pi/8','\pi/4','3\pi/8','\pi/2', '5*\pi/8', '3*\pi/4', '7*\pi/8', '\pi'});

% Add label too
ylabel('Angle between \Deltav and \Deltav_{true}','FontSize',24);

% Append legend (done here so extra lines have no legend entry)
legend({'Euclidean', 'Taylor', 'Eigen', 'Transport'}, 'Fontsize', 24);


%%% KL IMPROVEMENT PLOT

% Loop over the plots for different choices of maximum vector norm
for S = 1:size(KLimprovement_data,3)
    
    % Create the boxplot using the plotting subfunction
    plotFunction(problemL, squeeze(KLimprovement_data(:,:,S)), plot_colors);
    
    % Add the y label
    ylabel('log_{10} d_{SKL} Improvement','FontSize',24);
    
    % Note: y-limits could also be set here, but not knowing the values to
    % expect, this is not done and left to the user to re-specify manually
    % if required
    
    % Plot the zero improvement line as a reference
    xvals = get(gca,'xlim');
    xdiff = xvals(2) - xvals(1);
    xvals(1) = xvals(1) - xdiff * 0.05;
    xvals(2) = xvals(2) + xdiff * 0.05;
    plot(xvals, [0, 0], 'k', 'Linewidth', 2);

    % Append legend (done here so extra lines have no legend entry)
    legend({'Euclidean', 'Taylor', 'Eigen', 'Transport'}, 'Fontsize', 24);
    
end




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
xlabel('Geodesic Length (Binned)');
set(gca,'Fontsize', 24);

% Set up xlimits to cut off at the right spots
xlim([0-0.025*Lmax, Lmax*1.025]);