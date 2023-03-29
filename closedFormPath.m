function path = closedFormPath(p1, p2, Npts, path_type)
%
%     path = closedFormPath(p1, p2, path_type)
%
% This function outputs a non-geodesic path (cell array of Npts points)
% connecting the two input multivariate normals (MVNs) p1 and p2, each
% given as a struct defining mu and SIGMA for the distribution
%
% Several path types are available:
%
%     'euclidean' : A linear interpolation between the initial and final
%                   mu and SIGMA values
%
%     'geometric' : A linear interpolation between the log densities of the
%                   two input multivariate normals P1 and P2. Also given by
%                   the linear interpolation of their canonical (natural)
%                   parameters. Also known as an 'annealing' path.
%
%        'moment' : A path that linearly interpolates between the first and
%                   second moments of the two multivariate normals
%                   (see Grosse et al. 2013, cited in paper)
%
%        'hybrid' : A pragmatic path that simply averages the co-ordinates
%                   of the 'geometric' and 'moment' paths in attempt to
%                   balance their pros and cons
%
%     'transport' : Geodesics as defined by the 2-Wasserstein metric
%                   (optimal transport)
%
% If the path_type is not specified, the 'hybrid' path will be chosen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameterise the path using t = [0, 1]
tvec = linspace(0,1,Npts);

% Initialise the path object
path = cell(1,Npts);

% Read out the dimension (used for some methods)
d = size(p1.SIGMA);

switch lower(path_type)
    
    case {'euclid','euclidean'}
        
        for k = 1:Npts
            % Read out current time value (for notational ease)
            t = tvec(k);
            % Find linear interpolation of start and endpoints
            path{k}.mu = (1 - t) * p1.mu + t * p2.mu;
            path{k}.SIGMA = (1 - t) * p1.SIGMA + t * p2.SIGMA;
        end
        
    case {'geometric','geo','anneal','annealing','log','logdensity'}
        
        % Find natural parameters of start and end point
        DELTA1 = p1.SIGMA \ eye(d);
        DELTA2 = p2.SIGMA \ eye(d);
        delta1 = DELTA1 * p1.mu;
        delta2 = DELTA2 * p2.mu;
        
        for k = 1:Npts
            % Read out current time value (for notational ease)
            t = tvec(k);
            % Find linear interpolation of start and endpoints' natural
            % parameters
            deltat = (1 - t) * delta1 + t * delta2;
            DELTAt = (1 - t) * DELTA1 + t * DELTA2;
            % Convert back to regular parameters
            SIGMAt = DELTAt \ eye(d);
            path{k}.SIGMA = SIGMAt;
            path{k}.mu = SIGMAt * deltat;
            
        end
        
    case {'moment','moment_averaged','moments'}
        
        for k = 1:Npts
            % Read out current time value (for notational ease)
            t = tvec(k);
            % Find linear interpolation of mean
            path{k}.mu = (1 - t) * p1.mu + t * p2.mu;
            % Variance is a linear interpolation coupled with a stretch
            % along the direction between the means
            path{k}.SIGMA = (1 - t) * p1.SIGMA + t * p2.SIGMA + t * (1-t) * ( (p2.mu - p1.mu) * (p2.mu - p1.mu)' );
        end
        
    case {'hybrid','mixed','balanced'}
        
        % Most easily defined using recursive calls to this function
        path_moment = closedFormPath(p1,p2,Npts,'moment');
        path_anneal = closedFormPath(p1,p2,Npts,'anneal');
        
        for k = 1:length(path_moment)
            
            path{k}.mu = 0.5 * ( path_moment{k}.mu + path_anneal{k}.mu );
            path{k}.SIGMA = 0.5 * ( path_moment{k}.SIGMA + path_anneal{k}.SIGMA );
            
        end
        
        
    case {'transport','wasserstein','w2'}
        
        % Pre-calculate the matrix used in "affine interpolation"
        p2half = p2.SIGMA^(1/2);
        C = p2half * ( p2half * p1.SIGMA * p2half )^(-1/2) * p2half;
        
        % Has a closed form, given for example in Delon & Desolneux (2020)
        for k = 1:Npts
            % Read out current time value (for notational ease)
            t = tvec(k);
            % Optimal transport is a linear interpolation of mean
            path{k}.mu = (1 - t) * p1.mu + t * p2.mu;
            % Covariance matrix transforms by a set formula
            path{k}.SIGMA = ( (1-t) * eye(d) + t * C ) * p1.SIGMA * ( (1 - t) * eye(d) + t * C );
        end
        
end
