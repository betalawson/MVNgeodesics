function G = approximateGeodesic(p1, p2, method)
%
%     G = approximateGeodesic(p1, p2, method)
%
% This function returns a geodesic object that begins at the first input
% p1, a point in the space of multivariate normals specified as a struct
% defining the points mu and SIGMA co-ordiantes. The geodesic ends at a
% point that is (hopefully) close to the second input multivariate normal,
% p2 (specified similarly). Performance will typically be better when the 
% two input points p1 and p2 are close together.
%
% Geodesics in MVN space are available in closed form in terms of the
% initial velocity (tangent vector), but this is not available in closed
% form. Instead, an approximation to this initial velocity can be selected
% via a number of methods, specified by an input string:
%
%     'euclidean' : The velocity is that of the Euclidean path between the
%                   two points, such that it reaches p2 at t=1.
%
%        'taylor' : The velocity is selected according to the assumption
%                   that the velocity in mean space is small.
%
%         'eigen' : The velocity is selected by decomposing the required
%                   movement in mu-space into separate movements along each
%                   direction of a target covariance eigenvector.
%
%     'transport' : The velocity is calculated by forming two perfectly
%                   geodesic connections, from p1 to p* and then p* to p2.
%                   The velocity from p* to p2 is parallel transported
%                   backwards along the geodesic between p1 and p* to
%                   become a velocity correction from p1 to approximately
%                   hit p2.    ** COMPARABLY EXPENSIVE **
%
% This function is used inside the function onePointshooting() for numeric
% calculation of the geodesic between two points. If called directly, with
% no 'method' input specified, the 'eigen' method will be chosen.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default to the eigenvector-based method if none provided
if nargin < 3
    method = 'eigen';
end

% Check the two input points are valid, and output the dimension if so
[pts, D] = validatePoints({p1, p2});
% Re-assign the validated points to p1 and p2
[p1, p2] = deal(pts{:});

% As geodesic objects are specified in terms of their velocity at the
% origin, first affine transform the two points so the first is at the
% origin
[O, pt, P, r] = affineToOrigin(p1,p2);

% Select the initial velocity v = (x,B) for the geodesic connecting p1 and
% a point near p2 using the requested method
switch lower(method)
    
    % Euclidean approach
    case {'euclid','euclidean','straight'}
        
        % The velocity is the Euclidean difference between the transformed
        % points
        v.mu = pt.mu - O.mu;
        v.SIGMA = pt.SIGMA - O.SIGMA;

        
    % Taylor-series based approximation
    case {'taylor','expansion'}
                
        % Find the velocities given by the small-x approximation
        logSIGMA = logm(pt.SIGMA);
        v.mu = ( eye(D) - ( eye(D) - pt.SIGMA ) \ eye(D) ) * logSIGMA * ( pt.SIGMA \ pt.mu );
        v.SIGMA = logSIGMA;
        
        
    % Eigenvector-based approximation
    case {'eigen','decomposition'}
        
        %%% Further rotate so the target mean (canonical) aligns with a
        %%% standard basis vector
        
        % Find the canonical mean at the target
        DELTA_t = pt.SIGMA \ eye(D);
        delta_t = DELTA_t * pt.mu;
        
        % Find a rotation that will align this with a standard basis vector
        S = randomOrthogonalMatrix( delta_t );
        
        % Apply this rotation
        delta_t = S' * delta_t;
        DELTA_t = S' * DELTA_t * S;
        % Also update the tracking of all rotations used:
        %     S' * inv(P)   in origin space
        %     P * S         to get back to original
        P = P * S;
        
        % Find the eigendecomposition of the new target DELTA
        [V, LAM] = eig( DELTA_t );
        
        % Initialise the velocity at zero
        v.mu = zeros(D,1);
        v.SIGMA = zeros(D);
        
        % Add on the contribution from the different eigenvalue "pieces"
        for k = 1:D
            
            % Find amount of target mean in this direction
            dk = V(:,k)' * delta_t;
            nd2 = dk^2;
            % Read out current eigenvalue
            lambda = LAM(k,k);
            
            % Calculate the auxilliary quantities
            alpha = nd2 - 2*lambda^2 + 2*lambda;
            gamma = sqrt( alpha^2 + 8*nd2*lambda^2 );
            g = acosh( 1 + gamma^2 / (8*lambda^3) );
        
            % Only add a correction if its non-zero (avoids divide by 0)
            if abs(g) > 1e-13
                v.mu = v.mu + sign(dk) * g * (1 / sqrt(2)) * sqrt( 1 - (alpha/gamma)^2 ) * V(:,k);
                v.SIGMA = v.SIGMA + ( alpha/gamma * g ) * V(:,k) * V(:,k)';
            end
            
        end
        
        % Apply numerical protection to the covariance velocity
        v.SIGMA = numericalProtection(v.SIGMA);
        
        
    % Parallel transport based velocity approximation
    case {'transport','parallel'}
    
        %%% FIRST GEODESIC -  O = (0,I)  to  p* = (mu*,SIGMA_t)        
        
        % Find the eigendecomposition of the target covariance
        [V, LAM] = eig( pt.SIGMA );
        
        % Apply the rotation that diagonalises the covariance
        pt.mu = V' * pt.mu;
        pt.SIGMA = LAM;
        % Also note that this transformation was used - update P:
        %      V' * inv(P)   -   transformation to get to origin space
        %      P * V         -   transform back to original space
        P = P * V;
        
        % Get variance-weighted components of the mean in each direction
        components = sqrt(1./diag(LAM)) .* pt.mu;
        
        % Find the maximal direction
        [~, maxloc] = max( abs(components) );
                
        % Use this to get the intermediary point, p* = (mu*, SIGMAt)
        p_star.mu = zeros(D,1);
        p_star.mu(maxloc) =  pt.mu(maxloc);
        p_star.SIGMA = pt.SIGMA;
        
        % Find the canonical target (scalars as this is a 1-D solution)
        DELTA1 = 1 / LAM(maxloc,maxloc);
        delta1 = p_star.mu(maxloc) * DELTA1;
        
        % Get the 1-D solution - calculate auxillary quantities
        alpha = delta1^2 - 2 * DELTA1^2 + 2 * DELTA1;
        gamma = sqrt( alpha^2 + 8 * delta1^2 * DELTA1^2 );
        g = acosh( 1 + gamma^2 / (8 * DELTA1^3) );
        
        % Initialise the velocity - other directions just use logm(SIGMA)
        v1.mu = zeros(D,1);
        v1.SIGMA = diag( log( diag(LAM) ) );        
        
        % Add on the 1-D solution - if it's a nonzero update
        if norm(g) > 1e-13
            v1.mu(maxloc) = v1.mu(maxloc) + sign( delta1 ) * g * (1 / sqrt(2)) * sqrt( 1 - (alpha/gamma)^2 );
            v1.SIGMA(maxloc,maxloc) = v1.SIGMA(maxloc,maxloc) +  ( alpha/gamma * g + log(DELTA1) );
        end
        
        % Create this geodesic as a formal object for later use. This sits
        % in the origin space, so don't need the P and r that map back
        G1 = struct('v',v1,'P',eye(D),'r',zeros(D,1));
        
        %%% SECOND GEODESIC -  p* = (mu*,SIGMA_t)  to  pt = (mu_t,SIGMA_t)
        
        % Create a second transform of the space that shifts p* to origin
        [~,pt_tilde,P_tilde,~] = affineToOrigin( p_star, pt );

        % Find the known geodesic connection velocity in this new space
        nmu = norm(pt_tilde.mu);
        acosh_val = acosh( 1 + nmu^2 + 1/8 * nmu^4 ) / ( nmu * sqrt( nmu^2 + 8 ) );
        v2.mu = 2 * acosh_val * pt_tilde.mu;
        v2.SIGMA = acosh_val * ( pt_tilde.mu * pt_tilde.mu' );
        
        % Convert this back into a velocity in the original space
        v2.mu = P_tilde * v2.mu;
        v2.SIGMA = P_tilde * v2.SIGMA * P_tilde';
        
        %%% PARALLEL TRANSPORT TO UPDATE VELOCITY OF G1
        
        % Parallel transport the second geodesic's velocity back along G1
        dv = parallelTransport(v2, G1, [1 0]);
        
        % Find the Jacobi field at geodesic endpoint (t=1)
        J = jacobiField( G1, dv, 1 );
        
        % Calculate the step length by projecting the Euclidean residual
        % onto the Jacobi field (Jacobi is also a "Euclidean" approximation
        % to the derivative
        R.mu = pt.mu - p_star.mu;
        R.SIGMA = zeros(D);
        s = innerProduct(R, J, diag( 1 ./ diag(LAM) )) / innerProduct( J, J, diag( 1./ diag(LAM) ));
        
        % Final velocity approximation is G1 velocity plus the update
        v.mu = v1.mu + s * dv.mu;
        v.SIGMA = v1.SIGMA + s * dv.SIGMA;       
                    
    otherwise
        error('Incorrect specification of approximation method');

end

% Create a geodesic obejct using the calculated velocity, and also noting
% the transform to move back from "origin space" to true space
G = struct('v',v,'P',P,'r',r,'L',sqrt(innerProduct(v,v,eye(D))));

