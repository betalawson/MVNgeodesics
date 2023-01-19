function [G, diagnostics] = onePointShooting(p1, p2, supplied_options)
%
%     [G, diags] = onePointShooting( p1, p2 )
%     [G, diags] = onePointShooting( p1, p2, options )
%
% This function uses the shooting method for MVNs to find a geodesic G that
% connects the two provided points p1 and p2, which should be specified as
% structs with a 'mu' and 'SIGMA' component.
%
% The user may additionally supply options in the form of a struct, with
% fields matching the names given in OPTIONSdefaults.m and values as
% desired

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

% Diagnostics initialisation
tic;
s_vec = [];


%%% TRANSFORM START POINT TO ORIGIN

[O,pt,P,r] = affineToOrigin(p1,p2);


%%% INITIALISATION OF GEODESIC

% If a geodesic provided, use it as the initial geodesic - with initial
% velocity transformed appropriately
if isstruct(options.Ginit)
    
    % Grab out the input geodesic's velocity
    v = options.Ginit.v;
    
    % Transform this into the space we are now in
    v = affineTransform(v, inv(P) * options.Ginit.P, zeros(D,1));
    
    % Create a geodesic starting at the origin with this new velocity
    G = struct('v',v,'r',zeros(D,1),'P',eye(D));
    
% If a geodesic not provided, intialise at zero    
else
    
    % Zero velocity
    v = struct('mu',zeros(D,1),'SIGMA',zeros(D));
    % Geodesic in "origin space", so the transform's P and r are not used
    G = struct('v',v,'r',zeros(D,1),'P',eye(D));
    
end

%%% SHOOTING PROCEDURE

% Continually refine the geodesic until the symmeterised KL divergence
% between the end point and the target is sufficiently small, or maximum
% iteration count reached
looping = true;
iters = 0;
while looping
    
    % Fire the current geodesic out to its endpoint
    p = fireGeodesic(G,1);
    % Also find the precision it reached
    invSIGMA = p.SIGMA \ eye(D);
    
    % Plot the geodesic if requested
    if options.visualise
        
        % Figure preparation
        clf; hold(gca, 'on');
        
        % Plot the geodesic in the original space
        Gplot = struct( 'v', G.v, 'P', P, 'r', r );
        plotGeodesic(Gplot, 'scale', 0.15, 'Npts', 200);
        
        % Plot the start and end points
        plotMVN( p1 );
        plotMVN( p2 );
        
        if options.make_animation
            frame = getframe(gcf);
            writeVideo(vid_obj,frame);
        else
            % Pause so screen updates for user
            pause(0.1);
        end
    end
    
    % Calculate the current proxy distance to the target
    err_cur = sqrt( 2 * symKL( p, pt ) );
    err_vec(iters + 1) = err_cur;
    
    % Shooting iteration proceeds if current error is not sufficiently
    % small, and number of iterations not exceeded
    if iters < options.max_iters && err_cur > options.err_tol
        
        % Find the approximate geodesic connecting these points
        Gconnect = approximateGeodesic( p, pt, options.approx_method );
        
        % Find the norm of the velocity (easier calculation before transform)
        vconnect_norm = sqrt(innerProduct(Gconnect.v,Gconnect.v,eye(D)));
        
        % The velocity of Gconnect is in a space where p is the origin, so
        % transform this velocity back into our "origin space"
        vconnect = affineTransform( Gconnect.v, Gconnect.P, zeros(D,1) );
        
        % Find the euclidean connection from p to pt too
        R.mu = pt.mu - p.mu;
        R.SIGMA = pt.SIGMA - p.SIGMA;
        
        % Use backwards parallel transport to shift vconnect back to the origin
        % along the current G
        dv = parallelTransport(vconnect,G,[1 0]);
        
        % Find the Jacobi field associated with this shift to the geodesic G at
        % its endpoint
        J = jacobiField(G,dv,1);
        
        % Project the residual onto the Jacobi field to select a step length
        s = innerProduct(R,J,invSIGMA) / innerProduct(J,J,invSIGMA);
        
        % Reduce the step length when the velocity used for the update has too
        % large a magnitude (implying we should not over-index on this update)
        s = min( 1, options.max_vnorm / vconnect_norm ) * s;
        
        % Apply this velocity update
        G.v = struct('mu',G.v.mu + s*dv.mu,'SIGMA',G.v.SIGMA + s*dv.SIGMA);
        
        % Update the iteration count and store stepsize used
        iters = iters + 1;
        s_vec(iters) = s;
        
        % If verbose flag set, output the error for this iteration
        if options.verbose
            fprintf('Shooting... current proxy manifold distance is %g\n', err_cur);
        end
        
    % Otherwise, terminate loop and output a warning if terminating due to
    % iteration count
    else
        
        if iters >= options.max_iters 
            fprintf('\n -- WARNING! -- \n Geodesic shooting routine reached maximum number of iterations before achieving tolerance. \n Returned geodesic may not be correct!\n');            
            converged = false;
        elseif ~(err_vec(end) < options.err_tol)
            fprintf('\n -- WARNING! -- \n Geodesic shooting routine returned a geodesic that does not meet the tolerance. \n Blow-up has likely occurred (e.g. velocity went to NaN)\n');
            converged = false;    
        else
            converged = true;
        end
        looping = false;
        
    end
end

%%% OUTPUTS AND CLEANUP

% Close video 
if options.visualise && options.make_animation
    close(vid_obj);
end

% Geodesic was defined in origin space, so supply P and r so that it
% represents a geodesic in original space
G = struct('v',G.v,'P',P,'r',r,'L',sqrt(innerProduct(G.v,G.v,eye(D))));

% Output diagnostics
diagnostics.iterations = iters;
diagnostics.runtime = toc;
diagnostics.final_err = err_vec(end);
diagnostics.err_history = err_vec;
diagnostics.s_history = s_vec;
diagnostics.converged = converged;