% This MATLAB script tests the MEX interface of the gsTHBSpline class.
% Author: O. Chanon

% One should add the geopdes library to the path TODO
addpath( genpath('/Users/ondine/Documents/MATLAB/GeoPDEs-full/geopdes/') )

% Construct a truncated hierarchical basis by reading the specified file
filename = join([filedata,'surfaces/thbs_face_3levels.xml']);
fprintf('Reading THB spline from file: %s.\n',filename)
hbs = gsTHBSpline(filename);

% Print support
para = hbs.support();
fprintf('The parameter range is: [%f %f %f %f]\n', para);

% Print points
pts = uniformPointGrid(para(1:2),para(3:4),1000);

% Print evaluations
ev  = hbs.eval(pts);
fprintf('ev\n')
disp(size(ev))
figure;
scatter3(ev(1,:), ev(2,:), ev(3,:), '+')

% Print jacobian
jac = hbs.deriv(pts);
fprintf('jac')
disp(size(jac))

% Build GeoPDEs geometry structure
geometry = geo_load(hbs);

% Plot GeoPDEs geometry coming from G+smo
[X,Y] = ndgrid(0:0.01:1, 0:0.01:1);
reshape(X,1,[]);
X = reshape(X,1,[]);
Y = reshape(Y,1,[]);
pts = [X;Y];

ev2 = geometry.map(pts);
scatter3(ev2(1,:),ev2(2,:),ev2(3,:))

