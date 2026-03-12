function HPC_Figure10Data
% This function performs mass testing of geodesic shooting, generating data
% that is saved as a .mat file for use with the FIGURE10 visualisation code
% (this function is just a way to run on HPC)


%%% TESTING SPECIFICATIONS

% Specify the filename where the run data is stored (will be loaded from if
% it exists, and the regenerate flag was not manually set)
filename = 'DATA_shootingResidualsHPCtemp';

% Specify the random seed for reproducibility
rng(7);

% Specify list of dimensions to perform the test over
D_list = [2, 10];

% Number of geodesics per length to include in the test
N_geo = 5000;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Append .mat to the filename if it wasn't given already
if ~strcmp(filename(end-3:end),'.mat')
    filename = [filename, '.mat'];
end


%%% DATA LOADING/GENERATION

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


