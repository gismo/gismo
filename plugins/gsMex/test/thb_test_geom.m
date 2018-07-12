% This MATLAB script tests the MEX interface of the gsTHBSpline class.
% Author: O. Chanon
close all

% One should add the geopdes library to the path, and run set_path

%% TEST CONSTRUCTORS
% Construct a truncated hierarchical geometry by reading the specified file
filename = join([filedata, 'surfaces/thbs_face_3levels.xml']); % 'domain2d/squareTHB.xml']); %
fprintf('Reading THB spline from file: %s.\n',filename)
hbs = gsTHBSpline(filename);

% Get the gsTHBSplineBasis from which hbs is built
basis = hbs.basis;
% Get the control points from which hbs is built
coefs = hbs.coefs;

% Construct another truncated hierarchical geometry from the basis and
% the control points of the previous one. 
fprintf('Loading THB spline from basis and control points.\n')
hbs2 = gsTHBSpline(basis, coefs);

%% TEST ACCESSORS
fprintf('Dimension of the parametric space 1: %d\n',hbs.parDim);      % rdim in geopdes
fprintf('Dimension of the parametric space 2: %d\n',hbs2.parDim);
fprintf('Dimension of the physical space 1: %d\n',hbs.geoDim);        % ndim in geopdes
fprintf('Dimension of the physical space 2: %d\n',hbs2.geoDim);
fprintf('Size of the basis 1: %d\n',hbs.size);                        % ncomp in geopdes? TODO
fprintf('Size of the basis 2: %d\n',hbs2.size);

% Support
para = hbs.support;
fprintf('The parameter range 1 is: [%f %f %f %f]\n', para);
para2 = hbs2.support;
fprintf('The parameter range 2 is: [%f %f %f %f]\n', para2);

% Basis
deg1 = basis.degree(basis.dim());
fprintf('Degree 1 in the last direction: %d\n', deg1);
basis2 = hbs2.basis;
deg2 = basis2.degree(basis2.dim());
fprintf('Degree 2 in the last direction: %d\n', deg2);

% Control points
coefs_size = size(coefs);
fprintf('Number of basis functions 1: %d\n', prod(coefs_size(2:end)));
coefs2_size = size(hbs2.coefs());
fprintf('Number of basis functions 2: %d\n', prod(coefs2_size(2:end)));

%% TEST OTHER METHODS
% Print evaluations at pts
pts = uniformPointGrid(para(1:2),para(3:4),1000);
ev  = hbs.eval(pts);
fprintf('Evaluation 1: cf figure.\n')
figure;
subplot(1,2,1)
if (hbs.geoDim() == 3)
    scatter3(ev(1,:), ev(2,:), ev(3,:), '+')
elseif (hbs.geoDim() == 2)
    plot(ev(1,:),ev(2,:),'+')
end

ev2  = hbs2.eval(pts);
fprintf('Evaluation 12: cf figure.\n')
subplot(1,2,2);
if (hbs2.geoDim() == 3)
    scatter3(ev2(1,:), ev2(2,:), ev2(3,:), '+')
elseif (hbs2.geoDim() == 2)
    plot(ev2(1,:),ev2(2,:),'+')
end

%% TO FINISH
% Print jacobian
jac = hbs.deriv(pts);
fprintf('jac')
disp(size(jac))

% Print hessian on direction 1 
hess = hbs.hess(pts,1);
fprintf('hess')
disp(size(hess))

% Build GeoPDEs geometry structure
geometry = geo_load(hbs);
geometry2 = geo_load(hbs2);

% Get the knot vector corresponding to the first level, 2nd direction
kts12 = geometry.knots{1}{2};
fprintf('Knots level 1, direction 2\n')
disp(kts12);

kts12 = geometry2.knots{1}{2};
fprintf('Knots level 1, direction 2\n')
disp(kts12);

% Plot GeoPDEs geometry coming from G+smo
[X,Y] = ndgrid(0:0.01:1, 0:0.01:1);
reshape(X,1,[]);
X = reshape(X,1,[]);
Y = reshape(Y,1,[]);
pts = [X;Y];

figure;
ev21 = geometry.map(pts);
if (hbs.geoDim() == 3)
    scatter3(ev21(1,:),ev21(2,:),ev21(3,:))
elseif (hbs.geoDim() == 2)
    plot(ev21(1,:),ev21(2,:),'+')
end

figure;
ev22 = geometry2.map(pts);
if (hbs2.geoDim() == 3)
    scatter3(ev22(1,:),ev22(2,:),ev22(3,:))
elseif (hbs2.geoDim() == 2)
    plot(ev22(1,:),ev22(2,:),'+')
end
