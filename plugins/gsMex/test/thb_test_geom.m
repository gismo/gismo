% This MATLAB script tests the MEX interface of the gsTHBSpline class.
% Author: O. Chanon

%% TEST CONSTRUCTORS
% Construct a truncated hierarchical geometry by reading the specified file
filename = join([filedata, 'domain1d/thbcurve.xml']);
hbs1d = gsTHBSpline(filename, 1);
filename = join([filedata, 'domain2d/rectangleTHB.xml']);
hbs2d = gsTHBSpline(filename, 2);

% Copy constructor of a THB spline geometry.
hbs1d_copy = gsTHBSpline(hbs1d, 1);
assert(isequal(hbs1d_copy.coefs, hbs1d.coefs))
hbs2d_copy = gsTHBSpline(hbs2d, 2);
assert(isequal(hbs2d_copy.coefs, hbs2d.coefs))

% Get the gsTHBSplineBasis from which hbs is built
basis1d = hbs1d.basis;
basis2d = hbs2d.basis;
% Get the control points from which hbs is built
coefs1d = hbs1d.coefs;
coefs2d = hbs2d.coefs;

% Construct another truncated hierarchical geometry from the basis and
% the control points of the previous one.
hbs1d_copy = gsTHBSpline(basis1d, coefs1d, 1);
assert(isequal(hbs1d.support, hbs1d_copy.support))
hbs2d_copy = gsTHBSpline(basis2d, coefs2d, 2);
assert(isequal(hbs2d.support, hbs2d_copy.support))

fprintf('Test on constructors: passed.\n')

%% TEST ACCESSORS
assert(hbs1d.parDim==1); 
assert(hbs1d.geoDim==3);
assert(hbs1d.size==1);
assert(hbs2d.parDim==2); 
assert(hbs2d.geoDim==2); 
assert(hbs2d.size==1); 

% Support
para1d = hbs1d.support;
assert(isequal(para1d,[0 1]))
para2d = hbs2d.support;
assert(isequal(para2d,[0 1; 0 1]))

% Basis
deg1 = basis1d.degree(basis1d.dim());
assert(deg1==1);
deg1 = basis2d.degree(basis2d.dim());
assert(deg1==2);

% Control points
c2 = [0    0    0;
     0.6  0.1  0.1;
     1    1    0.5;
     2    0.4  0.45;
     2.2  0.3  0.4 ];
assert(prod(ismembertol(c2,coefs1d,'ByRows',true))==1)
clear c
[c(:,:,1),c(:,:,2)] = ndgrid(linspace(0,2,4), linspace(0,1,3));
c = reshape(c, [12 2]);
assert(prod(ismembertol(c,coefs2d,'ByRows',true))==1)

fprintf('Test on accessors: passed.\n\n')

%% TEST OTHER METHODS
% Print evaluations at pts
ev1 = hbs1d.eval(linspace(0,1,50));
pts = uniformPointGrid(para2d(1:2),para2d(3:4),1000);
ev2 = hbs2d.eval(pts);
figure;
subplot(1,3,1)
scatter3(ev1(1,:), ev1(2,:), ev1(3,:), '+')
subplot(1,3,2)
plot(ev2(1,:),ev2(2,:),'+')

% Compute jacobian
jac = hbs1d.jacobian(0.4);
assert(prod(ismembertol(jac,[1.6; 3.6; 1.6],'ByRows',true))==1)
jac = hbs2d.jacobian([0.5;0.2]);
assert(prod(ismembertol(jac,[4/3 0; 0 1],'ByRows',true))==1)

% Print hessian on direction 1
hess = hbs1d.hess(0.4,1);
assert(prod(ismembertol(hess,[0; 0; 0],'ByRows',true))==1)
hess = hbs2d.hess([0.5;0.2],1);
assert(prod(ismembertol(hess,[8/3; 0; 0; 0],'ByRows',true))==1)

hess = hbs2d.hess([0.5;0.2],hbs2d.parDim);
assert(prod(ismembertol(hess,zeros(4,1),'ByRows',true))==1)

% Slicing along the some direction (not possible for THBSplines of
% parametric dimension 1)
sl = hbs2d.sliceCoefs(1,0);
assert(prod(ismembertol(sl,[0 0;0 0.5;0 1],'ByRows',true))==1)

% Save geometry to xml file.
hbs1d.save('hbsgeom1d');
hbs2d.save('hbsgeom2d');
fprintf('Geometries saved to hbsgeom1d, hbsgeom2d and hbsgeom3d xml files.\n')

% Uniformly refine basis and change coefficients: obtain the same geometry
new_coefs = basis1d.uniformRefine_withCoefs(coefs1d,1,1);
hbs1d = gsTHBSpline(basis1d, new_coefs, 1);
new_coefs = basis2d.uniformRefine_withCoefs(coefs2d,1,1);
hbs2d = gsTHBSpline(basis2d, new_coefs, 2);

% Refine the basis by defining boxes and change coefficients: obtain the
% same geometry
boxes = [2,1,2];
basis1d = hbs1d_copy.basis; coefs1d = hbs1d_copy.coefs;
new_coefs1d = basis1d.refineElements_withCoefs(coefs1d, boxes);
hbs1d2 = gsTHBSpline(basis1d, new_coefs1d, 1);
boxes = [2,1,1,2,2];
basis2d = hbs2d_copy.basis; coefs2d = hbs2d_copy.coefs;
new_coefs2 = basis2d.refineElements_withCoefs(coefs2d,boxes);
hbs2d2 = gsTHBSpline(basis2d, new_coefs2, 2);

fprintf('All tests: passed.\n')
