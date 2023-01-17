function path = traceGeodesic(G,range,Npts)
%
%     p = traceGeodesic(G, range, Npts)
%     p = traceGeodesic(G, range)
%     p = traceGeodesic(G)
%
% This function traces the input geodesic G over the range of values given
% in the input argument 'range', split up into the number of points
% requested. The path is returned as a cell array of points, each a struct
% with a 'mu' and 'SIGMA' component.
%
% If the number of points is not specified, fifty will be used. If no range
% is specified, the geodesic will be traced out between t=0 and t=1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default values for any inputs not provided
if nargin < 2
    range = [0, 1];
end
if nargin < 3
    Npts = 50;
end

% Grab out the list of t values at which the geodesic should be evaluated
t_vals = linspace( range(1), range(2), Npts );
  
% Initialise path array
path = cell(1,Npts);
  
% Loop over each time value and evaluate the geodesic there
for k = 1:length(t_vals)
    
    path{k} = fireGeodesic( G, t_vals(k) );
  
end
