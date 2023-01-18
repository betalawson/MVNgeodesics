function FIGURE_multiPointPerformance(regenerate)
% This function behaves essentially as a script that generates figures
% showing the performance of different versions of the multi-point shooting
% algorithm

% Define the test problem
p1 = struct('mu', [1; 2] ,'SIGMA', [1, 0.1; 0.1, 10]);
p2 = struct('mu', [70; 35] ,'SIGMA', [10, -0.8; -0.8, 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume data is not to be regenerated unless specifically requested
if nargin < 1
    regenerate = false;
end

if ~exist('data_multiPointPerformance.mat','file') || regenerate

    % Move into the base folder
    cd('..');

    % Collect performance versus iterations for the multi-point algorithm run
    % with different settings
    [G_ES,diags_ES] = multiPointShooting(p1, p2, struct('approx_method','eigen','initial_path','euclid','verbose',true));
    [G_aES,diags_aES] = multiPointShootingAdaptive(p1, p2, struct('approx_method', 'eigen','initial_path','euclid','verbose',true));
    [G_EH,diags_EH] = multiPointShooting(p1, p2, struct('approx_method','eigen','initial_path','hybrid','verbose',true));
    [G_aEH,diags_aEH] = multiPointShootingAdaptive(p1, p2, struct('approx_method', 'eigen','initial_path','hybrid','verbose',true));

    % Use a convergent method to find the length
    success = [diags_ES.converged, diags_aES.converged, diags_EH.converged, diags_aEH.converged];
    Ls_found = [G_ES.L, G_aES.L, G_EH.L, G_aEH.L];
    use_method = find(success,1);
    if ~isempty(use_method)
        L_true = Ls_found(use_method);
    else
        warning('No methods converged. Will not be able to plot a comparison of path length to true geodesic length.');
    end
        
    % Move back to original folder
    cd('Figures');
    
    % Save the generated data
    save('data_multiPointPerformance.mat','p1','p2','L_true','diags_ES','diags_aES','diags_EH','diags_aEH','G_ES','G_aES','G_EH','G_aEH','-v7.3');
    
else
    
    % Load the saved data
    load('data_multiPointPerformance.mat','L_true','diags_ES','diags_aES','diags_EH','diags_aEH','G_ES','G_aES','G_EH','G_aEH');
    
end
    
%%% SYMMETERISED KL DIVERGENCE PLOT

% Initialise figure
figure; hold on;

% Plot performance of each method
plot_symKL(diags_ES,'--',[1.0 0.3 0.3]);
plot_symKL(diags_aES,'--',[0.3 0.3 1.0]);
plot_symKL(diags_EH,'-',[1.0 0.3 0.3]);
plot_symKL(diags_aEH,'-',[0.3 0.3 1.0]);

% Add labels and set fontsize
xlabel('Shooting Iterations');
ylabel('log10 d_s_K_L');
set(gca,'FontSize',24);

% Add legend
legend({'S-S','adapt-S-S','S-H','adapt-S-H','E-S','adapt-E-S','E-H','adapt-E-H'},'FontSize',24);
%legend({'Euclidean', 'Euclidean (Adaptive)', 'Hybrid', 'Hyrid (Adaptive)'},'FontSize',24);


%%% LENGTH PLOT

% Plot the log10 of the length
figure; hold on;

L_true

plot_deltaL(diags_ES,L_true,'--',[1.0 0.3 0.3]);
plot_deltaL(diags_aES,L_true,'--',[0.3 0.3 1.0], G_aES);
plot_deltaL(diags_EH,L_true,'-',[1.0 0.3 0.3]);
plot_deltaL(diags_aEH,L_true,'-',[0.3 0.3 1.0], G_aEH);

% Add labels and set fontsize
xlabel('Shooting Iterations');
ylabel('log10 Wasted Length');
set(gca,'FontSize',24);

% Add legend
legend({'Euclidean', 'Euclidean (Adaptive)', 'Hybrid', 'Hyrid (Adaptive)'},'FontSize',24);


function plot_symKL(data,lineSpec,color_vec)
% Subfunction that plots iterations versus symmeterised KL

% Read out the data to plot here
x = data.iters_history;
y = data.symKL_history;

% Remove any data where the symKL could not be numerically calculated
valid = ~isinf(y) & ~isnan(y) & (y >= 0);
x = x(valid);
y = y(valid);

% Plot the log10 of symKL versus the iterations
plot(x, log10(y), lineSpec, 'LineWidth',2,'Color',color_vec);


function plot_deltaL(data,L_true,lineSpec,color_vec,G)
% Subfunction that plots iterations versus "wasted length" (difference
% between sum of individual geodesic segments, and the true geodesic
% distance)

y = data.L_history;
x = data.iters_history(1:length(y));     % Cut off any iterations after leaving multi-point mode

% If a geodesic was provided, this corresponds to an adaptive method that
% will have used single-point shooting to refine its length estimate, so
% add on another data point corresponding to the final result
if nargin > 4
    y(end+1) = G.L;
    x(end+1) = data.iters_history(end);
end

% Plot the log10 of symKL versus the iterations
plot(x, log10(y-L_true), lineSpec, 'LineWidth',2,'Color',color_vec);