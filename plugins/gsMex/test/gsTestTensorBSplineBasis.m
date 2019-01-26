% This MATLAB script tests the MEX interface of the gsTensorBSplineBasis class.
% Author: O. Chanon

%% TEST CONSTRUCTORS
% Construct a tensor B-spline basis by reading the specified file
filename = join([filedata,'bspbasis/tpBSpline2_02.xml']);
bspb2d = gsTensorBSplineBasis(filename, 2);
filename = join([filedata,'bspbasis/tpBSpline3.xml']);
bspb3d = gsTensorBSplineBasis(filename, 3);

% Copy constructor test
bspb2_copy = gsTensorBSplineBasis(bspb2d, 2);
assert(isequal(bspb2_copy.support, bspb2d.support))
bspb3_copy = gsTensorBSplineBasis(bspb3d, 3);
assert(isequal(bspb3_copy.support, bspb3d.support))

% Construct a tensor B-spline basis from its knot vector.
knots_2build = {[0,0,0,1,2,2,2],[4,4,4,5,6,6,6]};
bspb2 = gsTensorBSplineBasis(knots_2build, 2);
knots_2build = {[0,0,0,1,2,2,2],[4,4,4,5,6,6,6], [0,0,1,1]};
bspb3 = gsTensorBSplineBasis(knots_2build, 3);

fprintf('Test on constructors: passed.\n')

%% TEST ACCESSORS
assert(bspb2d.domainDim==2)
assert(bspb3d.domainDim==3)
assert(bspb2d.numElements==16)
assert(bspb3d.numElements==64)
assert(bspb2d.size==42)
assert(bspb3d.size==252)
para = bspb2d.support;
assert(isequal(para,[0 1; 0 1]))
para = bspb3d.support;
assert(isequal(para,[0 1; 0 1; 0 1]))
fprintf('Test on accesssors: passed.\n')

%% TEST OTHER METHODS
% Check degree
assert(bspb2d.degree(1)==2)
assert(bspb3d.degree(bspb3d.domainDim)==2)

% Check active functions and evaluations at pts
pt = [0.5;0.5];
assert(isequal(bspb2d.active(pt),[15;16;17;21;22;23;27;28;29;33;34;35]))
assert(isequal(bspb2d.eval(pt),[1/12;1/12;0;1/3;1/3;0;1/12;1/12;0;0;0;0]))
pt2 = [0.5;0.5;0.5];
act3 = bspb3d.active(pt2);
assert(isequal(act3(1:4),[99;100;101;105]))
eval3 = bspb3d.eval(pt2);
assert(ismembertol(eval3(1:4)',[1/24 1/24 0 1/6],'ByRows',1))

% Check evaluations at a single active basis function
assert(bspb2d.evalSingle(15,pt)==1/12);
assert(bspb3d.evalSingle(105,pt2)==1/6);

% Save to output file
bspb2d.save('tbasis2d');
bspb3d.save('tbasis3d');
fprintf('Check the xml output files tbasis2d and tbasis3d.\n')
% tbasis2d.xml
%
% <xml>
%  <Basis type="TensorBSplineBasis2" id="0">
%   <Basis type="BSplineBasis" index="0">
%    <KnotVector degree="2">0 0 0 0.25 0.5 0.75 1 1 1 </KnotVector>
%   </Basis>
%   <Basis type="BSplineBasis" index="1">
%    <KnotVector degree="3">0 0 0 0 0.25 0.5 0.75 1 1 1 1 </KnotVector>
%   </Basis>
%  </Basis>
% </xml>
%
%
% tbasis3d.xml
%
% <xml>
%  <Basis type="TensorBSplineBasis3" id="0">
%   <Basis type="BSplineBasis" index="0">
%    <KnotVector degree="2">0 0 0 0.25 0.5 0.75 1 1 1 </KnotVector>
%   </Basis>
%   <Basis type="BSplineBasis" index="1">
%    <KnotVector degree="3">0 0 0 0 0.25 0.5 0.75 1 1 1 1 </KnotVector>
%   </Basis>
%   <Basis type="BSplineBasis" index="2">
%    <KnotVector degree="2">0 0 0 0.25 0.5 0.75 1 1 1 </KnotVector>
%   </Basis>
%  </Basis>
% </xml>
%

% Check knot vectors
assert(isequal(bspb2d.knots(1),[0,0,linspace(0,1,5),1,1]));
assert(isequal(bspb3.knots(bspb3.domainDim),knots_2build{3}));

% Uniformily refine the basis
bspb2d.uniformRefine(1,1);
assert(isequal(bspb2d.knots(1),[0,0,linspace(0,1,9),1,1]));
bspb3d.uniformRefine(1,1);
assert(isequal(bspb3d.knots(2),[0,0,0,linspace(0,1,9),1,1,1]));

% Uniformly refine the basis with update of coefficients
clear coefs
[coefs(:,:,1),coefs(:,:,2)] = ndgrid(1:6,1:7);
coefs = reshape(coefs,42,2);
coefs2 = bspb2_copy.uniformRefine_withCoefs(coefs,1,1);
clear coefs3
assert(isequal(coefs2(50,:),[6 3.5]))

% Refine the basis by defining boxes
bspb2_copy = gsTensorBSplineBasis(bspb2, 2);
boxes = [2,1,1,3,3];
bspb2.refineElements(boxes);
assert(bspb2.numElements==16);
boxes = [3,1,1,1,3,3,3];
bspb3d.refineElements(boxes);
assert(bspb3d.numElements==810);

fprintf('All tests: passed.\n')

