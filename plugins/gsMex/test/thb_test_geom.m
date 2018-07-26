% This MATLAB script tests the MEX interface of the gsTHBSpline class.
% Author: O. Chanon
close all

% One should add the geopdes and nurbs library to the path, and run set_path

%% TEST CONSTRUCTORS
% Construct a truncated hierarchical geometry by reading the specified file
filename = join([filedata, 'domain2d/rectangleTHB.xml']); % 'surfaces/thbs_face_3levels.xml']); % 'domain2d/squareTHB.xml']); %
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

fprintf('Slicing along the first direction, fixing it to 0:\n')
sl = hbs.sliceCoefs(1,0.);

%% TEST ACCESSORS
fprintf('Dimension of the parametric space 1: %d\n',hbs.parDim);      % rdim in geopdes
fprintf('Dimension of the parametric space 2: %d\n',hbs2.parDim);
fprintf('Dimension of the physical space 1: %d\n',hbs.geoDim);        % ndim in geopdes
fprintf('Dimension of the physical space 2: %d\n',hbs2.geoDim);
fprintf('Size of the function set 1: %d\n',hbs.size);                 % ncomp in geopdes? TODO
fprintf('Size of the function set 2: %d\n',hbs2.size);

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
fprintf('Number of basis functions 1: %d\n', coefs_size(1));
coefs2_size = size(hbs2.coefs);
fprintf('Number of basis functions 2: %d\n\n', coefs2_size(1));
assert(isequal(coefs,hbs2.coefs));

% Uniforimly refine basis and change coefficients
new_coefs = basis.uniformRefine_withCoefs(coefs,1,1);
new_coefs2 = basis2.uniformRefine_withCoefs(hbs2.coefs,1,1);
fprintf('Refinement adds a single knot with multiplicity 1 on each knot span.\n');
fprintf('Number of basis functions 1 after refinement: %d\n', size(new_coefs,1));
fprintf('Number of basis functions 2 after refinement: %d\n', size(new_coefs2,1));
fprintf('Number of knots 1 at level 1 direction 1 after refinement: %d\n', length(basis.knots(1,1)));
fprintf('Number of knots 2 at level 1 direction 1 after refinement: %d\n', length(basis2.knots(1,1)));

%% TEST OTHER METHODS
% Print evaluations at pts
pts = uniformPointGrid(para(1:2),para(3:4),1000);
ev  = hbs.eval(pts);
fprintf('Evaluation 1 on %d pts: cf figure.\n',length(pts))
figure;
subplot(1,2,1)
if (hbs.geoDim() == 3)
    scatter3(ev(1,:), ev(2,:), ev(3,:), '+')
elseif (hbs.geoDim() == 2)
    plot(ev(1,:),ev(2,:),'+')
end

ev2  = hbs2.eval(pts);
fprintf('Evaluation 2: cf figure.\n')
subplot(1,2,2);
if (hbs2.geoDim() == 3)
    scatter3(ev2(1,:), ev2(2,:), ev2(3,:), '+')
elseif (hbs2.geoDim() == 2)
    plot(ev2(1,:),ev2(2,:),'+')
end
savefig('thb_test_geom_eval')

% Print jacobian at pts
jac = hbs.jacobian(pts);
fprintf('Jacobian (:,1:5) 1 of total size %d x %d: \n', size(jac,1),size(jac,2))
disp(jac(:,1:5))

jac2 = hbs2.jacobian(pts);
fprintf('Jacobian (:,1:5) 2 of total size %d x %d: \n', size(jac2,1),size(jac2,2))
disp(jac2(:,1:5))

% Print hessian on direction 1 %% TODO !!!! NOT WELL IMPLEMENTED IN GISMO!!
hess = hbs.hess(pts,1);
fprintf('Hessian (:, 1:5) 1 in direction 1 of total size %d x %d: \n', size(hess,1),size(hess,2))
disp(hess(:,1:5))

hess2 = hbs2.hess(pts,1);
fprintf('Hessian (:, 1:5) 2 in direction 1 of total size %d x %d: \n', size(hess2,1),size(hess2,2))
disp(hess2(:,1:5))

hess2last = hbs2.hess(pts,hbs2.parDim);
fprintf('Hessian 2 in direction parDim has total size %d x %d. \n', size(hess2last,1),size(hess2last,2))

% Print active functions on pts %% TODO!! does what we want? what does it mean?
act = hbs.active(pts(:,131:133));
fprintf('Active functions 1 on three pts:\n')
disp(act)

act2 = hbs2.active(pts(:,131:133));
fprintf('Active functions 2 on three pts:\n')
disp(act2)
 
%% TEST GEOPDES LOADING OF A GISMO GEOMETRY
% Build GeoPDEs geometry structures
geometry = geo_load(hbs);
geometry2 = geo_load('/Users/ondine/Documents/geopdes/geopdes/inst/examples/geometry_files/geo_rectangle.txt');

% Get the knot vector corresponding to the first level, 2nd direction
kts12 = geometry.knots{1}{2};
fprintf('Knots level 1, direction 2 from G+smo to GeoPdes:\n')
disp(kts12);

kts12 = geometry2.nurbs.knots{2};
fprintf('Knots direction 2 in GeoPdes 2:\n')
disp(kts12);

% Plot GeoPDEs geometry coming from G+smo
[X,Y] = ndgrid(0:0.01:1, 0:0.01:1);
reshape(X,1,[]);
X = reshape(X,1,[]);
Y = reshape(Y,1,[]);
pts = [X;Y];

fprintf('Plot GeoPDEs geometries coming from G+smo: see figure.\n')

figure;
subplot(1,2,1)
ev21 = geometry.map(pts);
if (hbs.geoDim() == 3)
    scatter3(ev21(1,:),ev21(2,:),ev21(3,:))
elseif (hbs.geoDim() == 2)
    plot(ev21(1,:),ev21(2,:),'+')
end

subplot(1,2,2)
ev22 = geometry2.map(pts);
if (hbs2.geoDim() == 3)
    scatter3(ev22(1,:),ev22(2,:),ev22(3,:))
elseif (hbs2.geoDim() == 2)
    plot(ev22(1,:),ev22(2,:),'+')
end
savefig('thb_test_geom_map')
fprintf('Max difference: %f\n', max(max(abs(ev21-ev22))))

der21 = geometry.map_der(pts);
der22 = geometry2.map_der(pts);
fprintf('GeoPDEs 1st derivative geometries coming from G+smo: max difference: %f\n', max(max(max(abs(der21-der22)))))

hess21 = geometry.map_der2(pts);
hess22 = geometry2.map_der2(pts);
fprintf('GeoPDEs hessian coming from G+smo: max difference: %f\n', max(max(max(max(abs(hess21-hess22))))))
