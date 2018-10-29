% This MATLAB script tests the MEX interface of the gsTensorBSpline class.
% Author: O. Chanon
close all

% One should add the geopdes and nurbs library to the path, and run set_path

%% TEST CONSTRUCTORS
% Construct a tensor B-spline geometry by reading the specified file
filename = join([filedata, 'domain2d/lake.xml']); %'domain2d/rectangleTHB.xml']); % 
fprintf('Reading gsTensorBSpline from file: %s.\n',filename)
bsp = gsTensorBSpline(filename);

% Get the gsTensorBSplineBasis from which hbs is built
basis = bsp.basis;
% Get the control points from which hbs is built
coefs = bsp.coefs;

% Construct another tensor B-spline geometry from the basis and
% the control points of the previous one.
fprintf('Loading gsTensorBSpline from basis and control points.\n')
bsp2 = gsTensorBSpline(basis, coefs);

fprintf('Save to xml file.')
bsp.save('bspgeom_test');

%% TEST ACCESSORS
fprintf('Dimension of the parametric space 1: %d\n',bsp.parDim);      % rdim in geopdes
fprintf('Dimension of the parametric space 2: %d\n',bsp2.parDim);
fprintf('Dimension of the physical space 1: %d\n',bsp.geoDim);        % ndim in geopdes
fprintf('Dimension of the physical space 2: %d\n',bsp2.geoDim);
fprintf('Size of the function set 1: %d\n',bsp.size);                 % ncomp in geopdes? TODO
fprintf('Size of the function set 2: %d\n',bsp2.size);

% Support
para = bsp.support;
fprintf('The parameter range 1 is: [%f %f %f %f]\n', para);
para2 = bsp2.support;
fprintf('The parameter range 2 is: [%f %f %f %f]\n', para2);

% Basis
deg1 = basis.degree(basis.domainDim());
fprintf('Degree 1 in the last direction: %d\n', deg1);
basis2 = bsp2.basis;
deg2 = basis2.degree(basis2.domainDim());
fprintf('Degree 2 in the last direction: %d\n', deg2);

% Control points
coefs_size = size(coefs);
fprintf('Number of basis functions 1: %d\n', coefs_size(1));
coefs2_size = size(bsp2.coefs);
fprintf('Number of basis functions 2: %d\n\n', coefs2_size(1));
assert(isequal(coefs,bsp2.coefs));

% Uniforimly refine basis and change coefficients
new_coefs = basis.uniformRefine_withCoefs(coefs,1,1);
new_coefs2 = basis2.uniformRefine_withCoefs(bsp2.coefs,1,1);
fprintf('Refinement adds a single knot with multiplicity 1 on each knot span.\n');
fprintf('Number of basis functions 1 after refinement: %d\n', size(new_coefs,1));
fprintf('Number of basis functions 2 after refinement: %d\n', size(new_coefs2,1));
fprintf('Number of knots 1 in direction 1 after refinement: %d\n', length(basis.knots(1)));
fprintf('Number of knots 2 in direction 1 after refinement: %d\n', length(basis2.knots(1)));
bsp = gsTensorBSpline(basis, new_coefs);
bsp2 = gsTensorBSpline(basis2, new_coefs2);

%% TEST OTHER METHODS
% Print evaluations at pts
pts = uniformPointGrid(para(1:2),para(3:4),1000);
ev  = bsp.eval(pts);
fprintf('Evaluation 1 on %d pts: cf figure.\n',length(pts))
figure;
subplot(1,2,1)
if (bsp.geoDim() == 3)
    scatter3(ev(1,:), ev(2,:), ev(3,:), '+')
elseif (bsp.geoDim() == 2)
    plot(ev(1,:),ev(2,:),'+')
end

ev2  = bsp2.eval(pts);
fprintf('Evaluation 2: cf figure.\n')
subplot(1,2,2);
if (bsp2.geoDim() == 3)
    scatter3(ev2(1,:), ev2(2,:), ev2(3,:), '+')
elseif (bsp2.geoDim() == 2)
    plot(ev2(1,:),ev2(2,:),'+')
end
savefig('thb_test_geom_eval')

% Print jacobian at pts
jac = bsp.jacobian(pts);
fprintf('Jacobian (:,1:5) 1 of total size %d x %d: \n', size(jac,1),size(jac,2))
disp(jac(:,1:5))

jac2 = bsp2.jacobian(pts);
fprintf('Jacobian (:,1:5) 2 of total size %d x %d: \n', size(jac2,1),size(jac2,2))
disp(jac2(:,1:5))

% Print hessian on direction 1 %% TODO !!!! NOT WELL IMPLEMENTED IN GISMO!!
hess = bsp.hess(pts,1);
fprintf('Hessian (:, 1:5) 1 in direction 1 of total size %d x %d: \n', size(hess,1),size(hess,2))
disp(hess(:,1:5))

hess2 = bsp2.hess(pts,1);
fprintf('Hessian (:, 1:5) 2 in direction 1 of total size %d x %d: \n', size(hess2,1),size(hess2,2))
disp(hess2(:,1:5))

hess2last = bsp2.hess(pts,bsp2.parDim);
fprintf('Hessian 2 in direction parDim has total size %d x %d. \n', size(hess2last,1),size(hess2last,2))

% Print active functions on pts %% TODO!! does what we want? what does it mean?
act = bsp.active(pts(:,131:133));
fprintf('Active functions 1 on three pts:\n')
disp(act)

act2 = bsp2.active(pts(:,131:133));
fprintf('Active functions 2 on three pts:\n')
disp(act2)
 
%% TEST GEOPDES LOADING OF A GISMO GEOMETRY
% Build GeoPDEs geometry structures
%bsp = gsTensorBSpline(join([filedata, 'domain2d/lake.xml.xml']));
geometry = geo_load(bsp);

% Get the knot vector corresponding to the 2nd direction
kts2 = geometry.knots{2};
fprintf('Knots in direction 2 from G+smo to GeoPdes:\n')
disp(kts2);

% Plot GeoPDEs geometry coming from G+smo
[X,Y] = ndgrid(0:0.01:1, 0:0.01:1);
reshape(X,1,[]);
X = reshape(X,1,[]);
Y = reshape(Y,1,[]);
pts = [X;Y];

fprintf('Plot GeoPDEs geometries coming from G+smo: see figure.\n')

figure;
ev2 = geometry.map(pts);
if (bsp.geoDim() == 3)
    scatter3(ev2(1,:),ev2(2,:),ev2(3,:))
elseif (bsp.geoDim() == 2)
    plot(ev2(1,:),ev2(2,:),'+')
end
savefig('thb_test_geom_map')

der21 = geometry.map_der(pts);
hess21 = geometry.map_der2(pts);
