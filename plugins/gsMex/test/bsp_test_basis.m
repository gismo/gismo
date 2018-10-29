% This MATLAB script tests the MEX interface of the gsTensorBSplineBasis class.
% Author: O. Chanon

%% TEST CONSTRUCTORS
% Construct a tensor B-spline basis by reading the specified file
filename = join([filedata,'bspbasis/tpBSpline2_02.xml']);
fprintf('Reading gsTensorBSplineBasis from file: %s.\n',filename)
bspb = gsTensorBSplineBasis(filename);

% Construct a tensor B-spline basis from its knot vector.
knots_2build = {[0,0,0,1,2,2,2],[4,4,4,5,6,6,6]};
fprintf('Reading gsTensorBSpline basis from knot vector.\n')
bspb2 = gsTensorBSplineBasis(knots_2build);

%% TEST ACCESSORS
assert(bspb.domainDim==2)
assert(bspb2.domainDim==2)
assert(bspb.numElements==16)
assert(bspb2.numElements==4)
assert(bspb.size==42)
assert(bspb2.size==16)
para = bspb.support;
assert(isequal(para,[0 1; 0 1]))
para2 = bspb2.support;
assert(isequal(para2,[0 2; 4 6]))
fprintf('Test on accesssors: passed.\n')

%% TEST OTHER METHODS
% Check degree
assert(bspb.degree(1)==2) % order = degree-1 in geopdes
assert(bspb2.degree(bspb2.domainDim)==2)

% Check active functions and evaluations at pts
pt = [0.5;0.5];
assert(isequal(bspb.active(pt),[15;16;17;21;22;23;27;28;29;33;34;35]))
assert(isequal(bspb.eval(pt),[1/12;1/12;0;1/3;1/3;0;1/12;1/12;0;0;0;0]))

pt2 = [1;5];
ev2  = bspb2.eval(pt2);
assert(isequal(bspb2.active(pt2),[6;7;8;10;11;12;14;15;16]))
assert(isequal(bspb2.eval(pt2),[0.25;0.25;0;0.25;0.25;0;0;0;0]))

% Check evaluations at a single active basis function
assert(bspb.evalSingle(15,pt)==1/12);
assert(bspb2.evalSingle(6,pt2)==0.25);

% Save to output file
bspb.save('output1');
bspb2.save('output2');
fprintf('Check the output files output1.xml and output2.xml\n\n')
% output1.xml
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
% output2.xml
%
% <xml>
%  <Basis type="TensorBSplineBasis2" id="0">
%   <Basis type="BSplineBasis" index="0">
%    <KnotVector degree="2">0 0 0 1 2 2 2 </KnotVector>
%   </Basis>
%   <Basis type="BSplineBasis" index="1">
%    <KnotVector degree="2">4 4 4 5 6 6 6 </KnotVector>
%   </Basis>
%  </Basis>
% </xml>

% Check knot vectors
assert(isequal(bspb.knots(1),[0,0,linspace(0,1,5),1,1]));
assert(isequal(bspb2.knots(bspb2.domainDim),knots_2build{2}));

% Uniformily refine the basis
bspb.uniformRefine(1,1);
assert(isequal(bspb.knots(1),[0,0,linspace(0,1,9),1,1]));

fprintf('All tests: passed.\n')

