function [G, diagnostics] = multiPointShooting(p1, p2, supplied_options)
%
%     [G, diags] = multiPointShooting(p1, p2)
%     [G, diags] = multiPointShooting(p1, p2, options)
%
% This function uses the multi-point shooting method to find the geodesic
% between two provided points p1 and p2, which should be specified as
% structs with a 'mu' and 'SIGMA' component.
%
% The method functions by initialising a set of points along a curve
% between the two input points p1 and p2. Geodesics connecting the
% even-numbered and odd-numbered pairs of points are then used to
% successively refine all point locations until the points lie on the
% overall geodesic path between p1 and p2.
%
% The user may additionally supply options in the form of a struct, with
% fields matching the names given in OPTIONSdefaults.m and values as
% desired. Any options not provided will take the values in OPTIONSdefaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LOAD IN ALGORITHM OPTIONS
% Defaults will be used where custom options were not specified

% Initialise with the default options
options = OPTIONSdefaults();
% Overwrite any of the default options with options provided by user
if nargin > 2 && isstruct(supplied_options)
    options_list = fieldnames(supplied_options);
    for k = 1:length(options_list)
        options.(options_list{k}) = supplied_options.(options_list{k});
    end
end

% Prepare options for the one-point shooting used inside this algorithm
one_point_options = options;
% Turn off visualisation and text output for the single-point connections
one_point_options.visualise = false;
one_point_options.verbose = false;



%%% CODE INITIALISATION

% Validate the inputs and read out the dimension
[pts, D] = validatePoints({p1,p2});
[p1, p2] = deal(pts{:});

% Prepare figure if requested
if options.visualise
    figure('units','normalized','OuterPosition',[0.1 0.1 0.8 0.8]);
    
    if options.make_animation
        vid_obj = VideoWriter('output_anim.avi');
        vid_obj.FrameRate = 5;
        open(vid_obj);
    end
    
end

% Loop and diagnostics initialisation
iters = 0;
multi_iters = 0;
multi_time = tic;
symKL_vec = [];
L_vec = [];
symKL_cur = Inf;

%%% INITIALISE PATH

% Read out the number of points as it is referred to often
Npts = options.N_points;
% Convert this to an odd number if it is not (avoids misleading plots and
% extra shooting to calculate lengths)
if mod(Npts,2) == 0
    Npts = Npts + 1;
end

% Place the points along the requested path (specified in options)
path = closedFormPath(p1,p2,Npts,options.initial_path);

% Geodesics between points are not yet known, so create an empty cell to
% hold them
Gs = cell(1,Npts-1);

%%% MAIN LOOP

looping = true;
while looping
    
    %%% UPDATE ITERATION COUNT FOR MULTI-SHOOT LOOP
    multi_iters = multi_iters + 1;
    loop_iter_count = 0;
    
    %%% PREPARE LOOP FOR UPDATING ODD-INDEX TARGET POINTS
    update_indices = 3:2:Npts-1;
    
    for k = 1:length(update_indices)
        start_points{k} = path{ update_indices(k) - 1 };
        end_points{k} = path{ update_indices(k) + 1 };
        Gs_here{k} = Gs{ update_indices(k) - 1 };
    end
    
    %%% RUN PARALLEL LOOP FOR UPDATING ODD-INDEX TARGET POINTS
    fail_flags = zeros(1,length(update_indices));
    iter_counts = zeros(1,length(update_indices));
    parfor k = 1:length(update_indices)
        
        % Initialise shooting with the current geodesic
        options_here = one_point_options;
        options_here.Ginit = Gs_here{k};
        
        % Perform shooting, then new even point location is halfway along
        [G_found, diags] = onePointShooting( start_points{k}, end_points{k}, options_here );
        point_found = fireGeodesic(G_found,0.5);
        
        % Store geodesic and new point
        Gs_found{k} = G_found;
        points_found{k} = point_found;
        
        % Check that this point is valid
        if any([isnan(point_found.mu);isinf(point_found.mu)]) || any([isnan(point_found.SIGMA(:));isinf(point_found.SIGMA(:))])
            fail_flags(k) = true;
        end
        
        % Update iteration count
        iter_counts(k) = diags.iterations;
        
    end
    
    %%% CLEAN UP PARALLEL STRUCTURE BACK TO OVERALL STRUCTURE
    
    % Store the geodesics and new points found in the main arrays
    for k = 1:length(update_indices)
        Gs{ update_indices(k) - 1 } = Gs_found{k};
        path{ update_indices(k) } = points_found{k};
    end
    
    % Trip the loop to stop if any failure case was hit
    if any(fail_flags)
        looping = false;
    end
    
    % Update iteration count
    loop_iter_count = loop_iter_count + sum(iter_counts);
    
    
    %%% ONLY LOOP OVER EVEN-INDEX POINTS IF ODD-INDEX UPDATES SUCCESSFUL
    if looping
        
        %%% PREPARE LOOP FOR UPDATING EVEN-INDEX TARGET POINTS
        update_indices = 2:2:Npts-1;
        
        for k = 1:length(update_indices)
            start_points{k} = path{ update_indices(k) - 1 };
            end_points{k} = path{ update_indices(k) + 1 };
            Gs_here{k} = Gs{ update_indices(k) - 1 };
        end
        
        %%% RUN PARALLEL LOOP FOR UPDATING EVEN-INDEX TARGET POINTS
        fail_flags = zeros(1,length(update_indices));
        iter_counts = zeros(1,length(update_indices));
        parfor k = 1:length(update_indices)
            
            % Initialise shooting with the current geodesic
            options_here = one_point_options;
            options_here.Ginit = Gs_here{k};
            
            % Perform shooting, then new even point location is halfway along
            [G_found, diags] = onePointShooting( start_points{k}, end_points{k}, options_here );
            point_found = fireGeodesic(G_found,0.5);
            
            % Store geodesic and new point
            Gs_found{k} = G_found;
            points_found{k} = point_found;
            
            % Check that this point is valid
            if any([isnan(point_found.mu);isinf(point_found.mu)]) || any([isnan(point_found.SIGMA(:));isinf(point_found.SIGMA(:))])
                fail_flags(k) = true;
            end
            
            % Update iteration count
            iter_counts(k) = diags.iterations;
            
        end
        
        %%% CLEAN UP PARALLEL STRUCTURE BACK TO OVERALL STRUCTURE
        
        % Store the geodesics and new points found in the main arrays
        for k = 1:length(update_indices)
            Gs{ update_indices(k) - 1 } = Gs_found{k};
            path{ update_indices(k) } = points_found{k};
        end
        
        % Trip the loop to stop if any failure case was hit
        if any(fail_flags)
            looping = false;
        end
        
        % Update iteration count
        loop_iter_count = loop_iter_count + sum(iter_counts);
        
    end
    
    % Update total iteration count
    iters = iters + loop_iter_count;
    
    
    %%% PLOT IF REQUESTED
    % Iterations used for plotting are not counted
    
    if options.visualise
        
        % Figure preparation
        clf; hold on;
        
        try
            
            % Trace out the geodesics between all the new point locations
            for k = 1:2:Npts-1
                plotGeodesic( Gs{k}, [0 1], 'Npts', 10 );
            end
            
            % Also plot the full geodesic obtained from the first geodesic
            
            %plotGeodesic(Gs{1},[0 0.5*(Npts-1)],'Npts',25, 'pathColor', [1 0 0]);
            
            % Plot the positions of all points on top, including the target
            for k = 1:Npts
                plotMVN(path{k});
            end
            plotMVN(p2);
            drawnow;
            
        catch
            warning('Plotting failed. Current geodesic may numerically reach a negative eigenvalue');
        end
        
    end
    
    
    %%% CHECK TERMINATION CRITERIA
    
    % Fire out the geodesic approximated by the method
    pcheck = fireGeodesic( Gs{1}, 0.5*(Npts-1) );
    
    % Check the symmeterised KL between the end point and the target
    symKL_cur = symKL(pcheck, p2);
    
    % Add this to the vector
    symKL_vec(multi_iters) = symKL_cur;
    
    % Also store the number of iterations in total currently used
    iters_vec(multi_iters) = iters;
    
    % Also store the current length of the total path
    L_tot = 0;
    for k = 1:2:length(Gs)
        L_tot = L_tot + Gs{k}.L;
    end
    L_vec(multi_iters) = L_tot;
    
    % Output performance if on verbose mode
    if options.verbose
        fprintf('Current length is L = %g, used %g iterations this loop\n', L_tot, loop_iter_count);
    end
    
    %%% Terminate loop if...
    %
    % ...if tolerance is reached
    if symKL_vec(multi_iters) <= options.symKL_tol && symKL_vec(multi_iters) >= 0
        looping = false;
        converged = true;
        
    % ...if iterations limit exceeded
    elseif multi_iters >= options.max_multi_iters
        looping = false;
        converged = false;
    end
    
end

% Grab out the multi-point method's idea of the geodesic between start/end
G = struct('v',struct('mu',Gs{1}.v.mu * 0.5*(Npts-1),'SIGMA',Gs{1}.v.SIGMA * 0.5*(Npts-1)),'P',Gs{1}.P,'r',Gs{1}.r);    
G.L = sqrt(innerProduct(G.v,G.v,eye(length(G.v.mu)));
    
% Store diagnostics
diagnostics.runtime = toc(multi_time);
diagnostics.iterations = iters;
diagnostics.multi_iterations = multi_iters;
diagnostics.symKL_history = symKL_vec;
diagnostics.iters_history = iters_vec;
diagnostics.L_history = L_vec;
diagnostics.converged = converged;    