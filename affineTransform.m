function object = affineTransform( object, P, r )
%
%     object = affineTransform(object, P, r)
%
% This function applies the affine transform
%      mu    -> P mu + r
%      SIGMA -> P SIGMA P^T
% to the input object. The metric defining the statistical manifold of
% multivariate normal distributions is invariant to this affine
% transformation.
%
% Acceptable objects include:
%   - A single vector or point, a struct defining v = (v_mu, v_SIGMA)
%   - A collection of such vectors or points, provided as a cell array
%     (e.g. a path defined by an ordered collection of points)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if a singular object was provided
if ~iscell(object)
   
    object.mu = P * object.mu + r;
    object.SIGMA = P * object.SIGMA * P';
    
% Otherwise if a cell array was provided, try to transform all elements
else 
    
    for k = 1:length(object)
        object{k}.mu = P * object{k}.mu + r;
        object{k}.SIGMA = P * object{k}.SIGMA * P';
    end
    
end