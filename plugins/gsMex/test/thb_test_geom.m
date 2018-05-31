% This MATLAB script tests the MEX interface of the gsTHBSpline class.
% Author: O. Chanon
close all

% One should add the geopdes library to the path TODO
addpath( genpath('/Users/ondine/Documents/MATLAB/GeoPDEs-full/geopdes/') )

% Construct a truncated hierarchical basis by reading the specified file
filename = join([filedata, 'surfaces/thbs_face_3levels.xml']); % 'domain2d/squareTHB.xml']);
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
%plot(ev(1,:),ev(2,:),'+')

% Print jacobian
jac = hbs.deriv(pts);
fprintf('jac')
disp(size(jac))

% Print hessian on direction 1 % TODO once it's implemented in G+smo
% hess = hbs.hess(pts,1);
% fprintf('hess')
% disp(size(hess))

% Build GeoPDEs geometry structure
geometry = geo_load(hbs);

% Get the knot vector corresponding to the first level, 2nd direction
kts12 = geometry.knots{1}{2};
fprintf('Knots level 1, direction 2\n')
disp(kts12);

% Plot GeoPDEs geometry coming from G+smo
[X,Y] = ndgrid(0:0.01:1, 0:0.01:1);
reshape(X,1,[]);
X = reshape(X,1,[]);
Y = reshape(Y,1,[]);
pts = [X;Y];

ev2 = geometry.map(pts);
scatter3(ev2(1,:),ev2(2,:),ev2(3,:))

% Get the gsTHBSplineBasis from which hbs is built
basis = hbs.basis();
% Get the degree of the underlying B-splines in the first direction
deg1 = basis.degree(1);
fprintf('Degree in the first direction: %d\n', deg1);

