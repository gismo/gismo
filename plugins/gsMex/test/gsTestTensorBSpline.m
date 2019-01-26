% This MATLAB script tests the MEX interface of the gsTensorBSpline class.
% Author: O. Chanon

%% TEST CONSTRUCTORS
% Construct a tensor B-spline geometry by reading the specified file
filename = join([filedata, 'domain2d/lake.xml']);
bsp2d = gsTensorBSpline(filename, 2);
filename = join([filedata, 'volumes/GshapedVolume.xml']);
bsp3d = gsTensorBSpline(filename, 3);

% Copy constructor of a THB spline geometry.
bsp2d_copy = gsTensorBSpline(bsp2d, 2);
assert(isequal(bsp2d_copy.coefs, bsp2d.coefs))
bsp3d_copy = gsTensorBSpline(bsp3d, 3);
assert(isequal(bsp3d_copy.coefs, bsp3d.coefs))

% Get the gsTHBSplineBasis from which hbs is built
basis2d = bsp2d.basis;
basis3d = bsp3d.basis;
% Get the control points from which hbs is built
coefs2d = bsp2d.coefs;
coefs3d = bsp3d.coefs;

% Construct another truncated hierarchical geometry from the basis and
% the control points of the previous one.
bsp2d_copy = gsTensorBSpline(basis2d, coefs2d, 2);
assert(isequal(bsp2d.support, bsp2d_copy.support))
bsp3d_copy = gsTensorBSpline(basis3d, coefs3d, 3);
assert(isequal(bsp3d.support, bsp3d_copy.support))

fprintf('Test on constructors: passed.\n')

%% TEST ACCESSORS
assert(bsp2d.parDim==2); 
assert(bsp2d.geoDim==2); 
assert(bsp2d.size==1); 
assert(bsp3d.parDim==3); 
assert(bsp3d.geoDim==3); 
assert(bsp3d.size==1); 

% Support
para2d = bsp2d.support;
assert(isequal(para2d,[0 1; 0 1]))
para3d = bsp3d.support;
assert(isequal(para3d,[0 1; 0 1; 0 1]))

% Basis
deg1 = basis2d.degree(basis2d.domainDim);
assert(deg1==5);
deg1 = basis3d.degree(basis3d.domainDim);
assert(deg1==2);

% Control points
clear c
[c(:,:,1),c(:,:,2)] = ndgrid(linspace(0,2,4), linspace(0,1,3));
c = reshape(c, [12 2]);
assert(ismembertol(coefs2d(end-10,:),[0.2 0.8],'ByRows',true))
assert(ismembertol(coefs3d(1,:),[0.737 0.631 0.0],'ByRows',true))

fprintf('Test on accessors: passed.\n\n')

%% TEST OTHER METHODS
% Print evaluations at pts
pts = uniformPointGrid(para2d(1:2),para2d(3:4),1000);
ev2 = bsp2d.eval(pts);
[c3(:,:,:,1),c3(:,:,:,2),c3(:,:,:,3)] = ndgrid(linspace(0,1,40),...
    linspace(0,1,40),linspace(0,1,40));
ev3 = bsp3d.eval(reshape(c3,40^3, 3)');
figure;
subplot(1,2,1)
plot(ev2(1,:),ev2(2,:),'+')
subplot(1,2,2)
scatter3(ev3(1,:), ev3(2,:), ev3(3,:), '+')

% Compute jacobian
jac = bsp2d.jacobian([0.5;0.2]);
assert(prod(ismembertol(jac,[1.570660462256091 0.119909691384325;...
                            -0.360644358371707 1.077517973680908],...
                            'ByRows',true))==1)
jac = bsp3d.jacobian([0.7;0.7;0.7]);
assert(prod(ismembertol(jac,[0.436723000005131   0.124442999999963  0
                             1.725562999999736   0.008482999999526  0
                             0 0 1],'ByRows',true))==1)

% Print hessian on direction 1
hess = bsp2d.hess([0.5;0.2],1);
assert(prod(ismembertol(hess,[3.007106769205564; -0.942702710917056;...
    -0.942702710917056; -0.840231601625828],'ByRows',true))==1)
hess = bsp3d.hess([0.7;0.7;0.7],1);
assert(ismembertol(hess(1),-19.006609999967452))

% Save geometry to xml file.
bsp2d.save('tgeom2d');
bsp3d.save('tgeom3d');
fprintf('Geometries saved to tgeom2d and tgeom3d xml files.\n')

% Uniformly refine basis and change coefficients: obtain the same geometry
new_coefs = basis2d.uniformRefine_withCoefs(coefs2d,1,1);
bsp2d = gsTensorBSpline(basis2d, new_coefs, 2);
new_coefs = basis3d.uniformRefine_withCoefs(coefs3d,1,1);
bsp3d = gsTensorBSpline(basis3d, new_coefs, 3);

fprintf('All tests: passed.\n')

