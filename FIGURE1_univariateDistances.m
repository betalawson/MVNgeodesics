function FIGURE1_univariateDistances()
%
% This function evaluates the distances between univariate normals under
% two illustrative scenarios, that demonstrate how Fisher-Rao distances
% behave differently to the other statistical distance metrics. Both
% scenarios involve the distance between two normals with the same
% variance, but differing means.
%     Scenario One: Moving the means different distances apart
%     Scenario Two: Changing the variance of both normals together
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PREPARATION

% Prepare list of distance functions taking two normals P1 and P2 as inputs
%   Note: Normals are specified as per other code, as structs with fields
%         'mu' and 'SIGMA', but functions here are univariate-specific
D_functions = { @(P1,P2) D_TVD(P1,P2), ...
                @(P1,P2) D_kolmogorovSmirnov(P1,P2), ...
                @(P1,P2) D_hellinger(P1,P2), ...
                @(P1,P2) D_sqrtKL(P1,P2), ...                                  % Not a formal metric
                @(P1,P2) D_sqrtKL(P2,P1), ...                                  % Not a formal metric
                @(P1,P2) 0.5*D_sqrtKL(P1,P2) + 0.5*D_sqrtKL(P2,P1), ...        % Not a formal metric
                @(P1,P2) D_sqrtJS(P1,P2), ...
                @(P1,P2) D_MMD(P1,P2,1), ...                                   % Kernel bandwidth h = 1
                @(P1,P2) D_energy(P1,P2), ...
                @(P1,P2) D_wasserstein(P1,P2), ...
                @(P1,P2) D_fisherRao(P1,P2)
                };

% Accompanying list of names for figure caption
D_names = {'d_{T\!V}', 'd_{K\!S}', 'd_{H}', '\sqrt{D_{K\!L}}', '(Reverse KL)', '(Sym. KL)', '\sqrt{D_{J\!S}}', 'd_{M\!M\!D}', 'd_{E}', 'd_{W}', 'd_{F\!R}'};

% List the indices of the distances listed to actually plot, in order
D_inds = [1,3,8,9,10,11,4];

% Specify the colours of the plotted distances, in order
plot_clrs = [ [0.0, 0.0, 0.6];            % TVD
              [0.4, 0.4, 1.0];            % Hellinger
              [0.8, 0.0, 0.0];            % MMD
              [1.0, 0.4, 0.4];            % E
              [0.9, 0.6, 0.1];            % W
              [1.0, 0.4, 1.0];            % F-R
              [0.4, 0.4, 0.4]  ];         % sqrtKL

% Specify whether to plot the scenarios on separate figures, or as subfigs
use_subfigs = false;
          


%%% SCENARIO ONE: SHIFTING THE MEAN

% Scenario specifications
max_mu = 3;
N_mu = 101;
mu_vec = linspace(0,max_mu,N_mu);
set_var = 0.25;

% Create the base normal, with mu = 0
P0 = struct('mu',0,'SIGMA',set_var);

% Pre-allocate storage for distance data
D = zeros(length(D_inds), N_mu);

% Loop over the listed discrepancies to be calculated
for d = 1:length(D_inds)
    
    % Loop over the different mu values and calculate discrepancy for each
    for m = 1:length(mu_vec)
       
        % Place second normal at current mu value (P0's variance inherited)
        P1 = setfield(P0, 'mu', mu_vec(m));
        
        % Calculate distance between P0 and P1 for current distance measure
        D(d,m) = D_functions{D_inds(d)}(P0,P1);
        
    end
end
          
% Figure/subfigure preparation
if use_subfigs
    figure('units','Normalized','OuterPosition',[0 0 1 1]);
    tiles_obj = tiledlayout(1,2);
    nexttile; 
    hold on;
else
    figure('units','Normalized','OuterPosition',[0.2 0.2 0.4 0.7]);
    hold on;
end

% Plot curves of d vs mu for the different distance measures
for d = 1:length(D_inds)
    if ismember(D_inds(d),[4,5,6])
        plot(mu_vec,D(d,:),'--','LineWidth',2.5,'Color',plot_clrs(d,:));
    else
        plot(mu_vec,D(d,:),'LineWidth',2.5,'Color',plot_clrs(d,:));
    end
    leg_txt{d} = ['$',D_names{(D_inds(d))},'$'];
    legend(leg_txt,'Box','Off','FontSize',24,'Interpreter','LaTeX','Location','NorthWest');
end

% Clean up plot
set(gca,'FontSize',24, 'LineWidth', 2.5);
xlabel('$\Delta\mu$', 'FontSize', 32, 'Interpreter', 'LaTeX');
ylabel('Distance', 'FontSize', 32);
ylim([0 5]);
axis square;


%%% SCENARIO TWO: ADJUSTING VARIANCE, MEAN SEPARATION FIXED

% Scenario specifications
mean_sep = 2.5;
N_var = 201;
logvar_min = -5;
logvar_max = 5;
var_vec = exp( linspace(logvar_min, logvar_max, N_var) );

% Pre-allocate storage for distance data
D = zeros(length(D_inds), N_var);

% Loop over the listed discrepancies to be calculated
for d = 1:length(D_inds)
    
    % Loop over the different variances and calculate discrepancy for each
    for m = 1:length(var_vec)
       
        % Create the base normal with current variance
        P0 = struct('mu',0,'SIGMA',var_vec(m));
        
        % Place second normal at given mean separation (inherit variance)
        P1 = setfield(P0, 'mu', mean_sep);
        
        % Calculate distance between P0 and P1 for current distance measure
        D(d,m) = D_functions{D_inds(d)}(P0,P1);
        
    end
end
          
% Figure/subfigure preparation
if use_subfigs
    nexttile;
    hold on;
else
    figure('units','Normalized','OuterPosition',[0.2 0.2 0.4 0.7]);
    hold on;
end

% Plot curves of d vs variance for the different distance measures
for d = 1:length(D_inds)
    if ismember(D_inds(d),[4,5,6])
        plot(-log(var_vec),D(d,:),'--','LineWidth',2.5,'Color',plot_clrs(d,:));
    else
        plot(-log(var_vec),D(d,:),'LineWidth',2.5,'Color',plot_clrs(d,:));
    end
end

% Clean up plot
set(gca,'FontSize',24, 'LineWidth', 2.5);
xlabel('$\log(\sigma^{-2})$', 'FontSize', 32, 'Interpreter', 'LaTeX');
ylabel('Distance', 'FontSize', 32);
ylim([0 5]);
xlim([-logvar_max,-logvar_min]);
axis square;


%%% SCHEMATIC DIAGRAM - FIRST SCENARIO

% Prepare the figure
figure('units','Normalized','OuterPosition',[0.4 0.3 0.6 0.7]);
hold on;

% Set the x limits to be centred on the midpoint between the Gaussians
x_limits = mean_sep/2 + mean_sep * [-2, 2];

% Plot the shifted second Gaussians
for i = 4:-1:1
    mu = mean_sep + 0.25*i;
    fplot(@(x) normpdf(x,mu,0.5), 'LineWidth', 2, 'Color', [1.0 0.0 0.0] + i * [0 0.2 0.2]);
end
% Plot the second Gaussian base
fplot(@(x) normpdf(x,mean_sep,0.5), 'LineWidth',4.5, 'Color', [1 0 0]);
% Plot the first Gaussian
fplot(@(x) normpdf(x,0,0.5), 'LineWidth',4.5, 'Color', [0 0 0]);

% Clean plot
xlim(x_limits);
ylim([0 1]);
axis off;


%%% SCHEMATIC DIAGRAM - SECOND SCENARIO

% Prepare the figure
figure('units','Normalized','OuterPosition',[0.4 0.3 0.6 0.7]);
hold on;

% Set the x limits to be centred on the midpoint between the Gaussians
x_limits = mean_sep/2 + mean_sep * [-2, 2];

% Plot the Gaussians of different variance levels
for i = 4:-1:1
    sigma = exp(log(0.5) + 0.25*i);
    fplot(@(x) normpdf(x,0,sigma), x_limits, 'LineWidth', 2, 'Color', [0.3 0.3 0.3] + i * [0.15 0.15 0.15]);
    fplot(@(x) normpdf(x,mean_sep,sigma), x_limits, 'LineWidth', 2, 'Color', [1.0 0.0 0.0] + i * [0.0 0.2 0.2]);
end
% Plot the second Gaussian base
fplot(@(x) normpdf(x,mean_sep,0.5), 'LineWidth',4.5, 'Color', [1 0 0]);
% Plot the first Gaussian base
fplot(@(x) normpdf(x,0,0.5), 'LineWidth',4.5, 'Color', [0 0 0]);

% Clean plot
xlim(x_limits);
ylim([0 1]);
axis off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the functions measuring statistical distance between two normals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% TOTAL VARIATION (ESTIMATED)
function D = D_TVD(P1,P2)

% Prepare integration bounds (most distant points, mu plus/minus 10 sigma)
sigma_points = [P1.mu - 10*sqrt(P1.SIGMA), P1.mu + 10*sqrt(P1.SIGMA), P2.mu - 10*sqrt(P2.SIGMA), P2.mu + 10*sqrt(P2.SIGMA)];
x_start = min(sigma_points);
x_end = max(sigma_points);

% Use a tiny fraction of the smaller variance as the discretisation width
dx = 0.0001 * min([P1.SIGMA, P2.SIGMA]);
% Split the domain into enough points to roughly have this separation
x_pts = linspace(x_start, x_end, ceil((x_end - x_start)/dx)+1);

% Numerically integrate the 1-norm between the normal pdfs
D = 0.5 * trapz(x_pts, abs( normpdf(x_pts,P1.mu,sqrt(P1.SIGMA)) - normpdf(x_pts,P2.mu,sqrt(P2.SIGMA))));

end

%%% KOLMOGOROV-SMIRNOV
function D = D_kolmogorovSmirnov(P1,P2)

% Use the negative squared difference between CDFs as a way to smoothly
% find the point where they maximally differ (initialised at zero)
x_pt = fminsearch( @(x) -(normcdf(x,P1.mu,sqrt(P1.SIGMA)) - normcdf(x,P2.mu,sqrt(P2.SIGMA)))^2, 0);

% Evaluate the actual absolute difference in CDFs at this point
D = abs( normcdf(x_pt,P1.mu,sqrt(P1.SIGMA)) - normcdf(x_pt,P2.mu,sqrt(P2.SIGMA)));

end

%%% MAXIMUM MEAN DISCREPANCY (Gaussian kernel)
function D = D_MMD(P1,P2,h)

% Closed formula for Gaussian-kernel MMD between two univariate normals
%   (Provided by generative AI but numerically validated)
D = sqrt(h ./ sqrt(h^2 + 2*P1.SIGMA) + h ./ sqrt(h^2 + 2*P2.SIGMA) - 2 * h / sqrt(h^2 + P1.SIGMA + P2.SIGMA) * exp( -(P1.mu - P2.mu)^2 / ( 2 * (h^2 + P1.SIGMA + P2.SIGMA) ) ));

end

%%% ENERGY DISTANCE (Euclidean metric in sample space)
function D = D_energy(P1,P2)

% Closed formula for energy distance between two univariate normals
%   (Provided by generative AI but numerically validated)
D = sqrt( 2 * (sqrt( P1.SIGMA + P2.SIGMA) * sqrt(2/pi) * exp( -(P1.mu - P2.mu)^2 / (2 * (P1.SIGMA + P2.SIGMA))) + abs(P1.mu-P2.mu) * erf( abs(P1.mu-P2.mu)/sqrt(2*(P1.SIGMA+P2.SIGMA)))) - 2*sqrt(P1.SIGMA/pi) - 2*sqrt(P2.SIGMA/pi));

end

%%% SQUARE ROOT KULLBACK-LEIBLER DIVERGENCE
function D = D_sqrtKL(P1,P2)

% Closed formula for the KL divergence between univariate normals
D = sqrt( 0.5 * ( P1.SIGMA / P2.SIGMA + (P2.mu - P1.mu)^2 / P2.SIGMA - 1 + log( P2.SIGMA / P1.SIGMA ) ) );

end

%%% SQUARE ROOT JENSEN-SHANNON DIVERGENCE (ESTIMATED)
function D = D_sqrtJS(P1,P2)

% Number of Monte Carlo samples
N = 1e7;
% Generate samples from the two distributions
x1 = P1.mu + sqrt(P1.SIGMA) * randn(N,1);
x2 = P2.mu + sqrt(P2.SIGMA) * randn(N,1);
% Monte Carlo estimation of the Jensen-Shannon divergence formula:
% D_JS = 1/2 E_p [log(p) - log(p/2+q/2)] + 1/2 E_q [log(q) - log(p/2+q/2)]
D = sqrt( 1/2 * 1/N * sum( log(normpdf(x1,P1.mu,sqrt(P1.SIGMA))) - log(0.5*normpdf(x1,P1.mu,sqrt(P1.SIGMA))+0.5*normpdf(x1,P2.mu,sqrt(P2.SIGMA)))) + 1/2 * 1/N * sum( log(normpdf(x2,P2.mu,sqrt(P2.SIGMA))) - log(0.5*normpdf(x1,P1.mu,sqrt(P1.SIGMA))+0.5*normpdf(x1,P2.mu,sqrt(P2.SIGMA))) )  );

end

%%% HELLINGER DISTANCE
function D = D_hellinger(P1,P2)

% Closed formula for the Hellinger distance between univariate normals
D = sqrt(1 - sqrt( 2 * sqrt(P1.SIGMA) * sqrt(P2.SIGMA) / (P1.SIGMA + P2.SIGMA) ) * exp( -0.25 * (P2.mu - P1.mu)^2 / ( P1.SIGMA + P2.SIGMA ) ));

end

%%% 2 - WASSERSTEIN DISTANCE (Euclidean metric in sample space)
function D = D_wasserstein(P1,P2)

% 2-Wasserstein is Euclidean distance in (mu, sigma) space
D = sqrt((P2.mu - P1.mu)^2 + (sqrt(P2.SIGMA) - sqrt(P1.SIGMA))^2);

end

%%% FISHER-RAO DISTANCE
function D = D_fisherRao(P1,P2)

% Closed formula for the Fisher-Rao distance between univariate normals
D = 2 * sqrt(2) * atanh(sqrt( ( 0.5*(P2.mu - P1.mu)^2 + (sqrt(P2.SIGMA) - sqrt(P1.SIGMA))^2 ) / ( 0.5*(P2.mu - P1.mu)^2 + (sqrt(P1.SIGMA)+sqrt(P2.SIGMA))^2) ) );

end
