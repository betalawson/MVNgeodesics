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
%    'projection' : The velocity is calculated using the embedding into a
%                   higher-dimensional manifold, for which geodesic paths
%                   are known but leave the MVN manifold (Nielsen 2023)
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
[pts, d] = validatePoints({p1, p2});
% Re-assign the validated points to p1 and p2
[p1, p2] = deal(pts{:});

% As geodesic objects are specified in terms of their velocity at the
% origin, first affine transform the two points so the first is at the
% origin
[O, pt, P, r] = affineToOrigin(p1,p2, true);

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
        v.mu = ( eye(d) - ( eye(d) - pt.SIGMA ) \ eye(d) ) * logSIGMA * ( pt.SIGMA \ pt.mu );
        v.SIGMA = logSIGMA;
                        
    % Eigenvector-based approximation
    case {'eigen','decomposition'}
        
        %%% Further rotate so the target mean (canonical) aligns with a
        %%% standard basis vector
        
        % Find the canonical mean at the target
        LAMBDA_t = pt.SIGMA \ eye(d);
        eta_t = LAMBDA_t * pt.mu;
                
        % Find the eigendecomposition of the new target LAMBDA
        [V, LAM] = eig( numericalProtection(LAMBDA_t) );        
        
        % Initialise the velocity at zero
        v.mu = zeros(d,1);
        v.SIGMA = zeros(d);
        v.mu = zeros(d,1);
        v.SIGMA = zeros(d);
        
        % Extract the unique eigenvalues, to some small level of tolerance,
        % so that we can process only one update for each projection of the
        % target mean onto the eigenbasis
        LAMBDAeigs = diag(LAM);
        [Ueigs, ~, Ugroup] = unique( round(LAMBDAeigs,10), 'stable' );        
        
        % Loop over each eigenbasis
        for k = 1:length(Ueigs)
            
            % Extract the actual eigenvalue for this basis
            lambda = Ueigs(k);
            
            % Extract out the matrix of eigenvectors associated with it
            Vk = V(:,Ugroup==k);
            
            % Calculate how much of the target mean lies in this direction
            pk = Vk' * eta_t;
            nd2 = norm(pk)^2;
            
            % Calculate the auxilliary quantities
            alpha = nd2 - 2*lambda^2 + 2*lambda;
            beta = sqrt( alpha^2 + 8*nd2*lambda^2 );
            gamma = acosh( 1 + beta^2 / (8*lambda^3) ) / beta;
        
            % Add on the "axis-aligned solution" for this eigenspace
            v.mu = v.mu + 2 * gamma * lambda * Vk*pk;
            v.SIGMA = v.SIGMA + Vk * ( ( alpha * gamma + log(lambda) ) * (pk*pk') / nd2  - log(lambda)*eye(length(pk)) ) * Vk';
                        
        end
        
        % Apply numerical protection to the covariance velocity
        v.SIGMA = numericalProtection(v.SIGMA);
        
    case {'path'} %%% NOT A REAL GEODESIC APPROXIMATOR, BUT HERE TO MAKE CODE EASIER - THIS ONE GENERATES AN APPROXIMATE PATH LENGTH ONLY
        
        % Generate the closed-form, embedding-based path to the target
        path = closedFormPath(O,pt,10001,'projection');
        % Calculate its length
        L = pathLength(path);
        % Store this as a geodesic velocity
        v.mu = L*[1;zeros(length(pt.mu)-1,1)];
        v.SIGMA = zeros(size(pt.SIGMA));
        
    case {'projection','con','tracemetric'}
        
        % Find the rate of change in the projected space
        P_t = [ pt.SIGMA + pt.mu * pt.mu', pt.mu; pt.mu', 1 ];
        dPdt = logm( P_t );
        % Convert this back into a velocity in the MVN space
        v.mu = dPdt(1:end-1,end);
        v.SIGMA = dPdt(1:end-1,1:end-1);
                
    case {'ansatz'}
        
        % Calculate the commutator that measures eigenvector misalignment
        C = 1/(12.5 + 4*log(1+1/6*(pt.mu'*(pt.SIGMA\pt.mu)))) * ( logm(pt.SIGMA) * pt.mu * pt.mu' - pt.mu * pt.mu' * logm(pt.SIGMA));
        % Add the symmetric component
        C = C - 0.5 * pt.mu * pt.mu';
        % Put this into the solution generator matrix
        Dt = [pt.SIGMA \ eye(d), zeros(d,1), zeros(d); zeros(1,d), 1, zeros(1,d); zeros(d), zeros(d,1), pt.SIGMA];
        Lt = [eye(d), zeros(d,1), zeros(d); pt.mu', 1, zeros(1,d); C, -pt.mu, eye(d)];
        % Take the logarithm of the full matrix, and read out the blocks
        A = logm( Lt * Dt * Lt' );
        v.mu = A(1:d,d+1);
        v.SIGMA = -A(1:d,1:d);
        
        
    otherwise
        error('Incorrect specification of approximation method');

end

% Create a geodesic obejct using the calculated velocity, and also noting
% the transform to move back from "origin space" to true space
G = struct('v',v,'P',P,'r',r,'L',sqrt(innerProduct(v,v,eye(d))));

