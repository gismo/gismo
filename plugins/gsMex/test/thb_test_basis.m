% This MATLAB script tests the MEX interface of the gsTHBSplineBasis class.
% Author: O. Chanon, P. Noertoft

%% TEST CONSTRUCTORS
% Construct a truncated hierarchical basis by reading the specified file
filename = join([filedata,'thbbasis/basis1d.xml']);
hbs1d = gsTHBSplineBasis(filename, 1);
filename = join([filedata,'thbbasis/simple.xml']);
hbs2d = gsTHBSplineBasis(filename, 2);
filename = join([filedata,'thbbasis/basis3d.xml']);
hbs3d = gsTHBSplineBasis(filename, 3);

% Copy constructor test
hbs1_copy = gsTHBSplineBasis(hbs1d, 1);
assert(isequal(hbs1_copy.support, hbs1d.support))
hbs2_copy = gsTHBSplineBasis(hbs2d, 2);
assert(isequal(hbs2_copy.support, hbs2d.support))
hbs3_copy = gsTHBSplineBasis(hbs3d, 3);
assert(isequal(hbs3_copy.support, hbs3d.support))

% Construct a truncated hierarchical basis from its knot vector in a 
% carthesian product way.
knots_2build = {[0,0,0,1,2,2,2]};
hbs1 = gsTHBSplineBasis(knots_2build, 1);
knots_2build = {[0,0,0,1,2,2,2],[4,4,4,5,6,6,6]};
hbs2 = gsTHBSplineBasis(knots_2build, 2);
knots_2build = {[0,0,0,1,2,2,2],[4,4,4,5,6,6,6],[4,4,4,5,6,6,6]};
hbs3 = gsTHBSplineBasis(knots_2build, 3);

fprintf('Test on constructors: passed.\n')

%% TEST ACCESSORS
assert(hbs1d.dim==1)
assert(hbs1d.numElements==4)
assert(hbs1d.size==5)
assert(hbs1d.treeSize==1) 
assert(hbs1d.treeLeafSize==1) 
para1 = hbs1d.support;
assert(isequal(para1,[0 1]))
assert(hbs1d.maxLevel==1)
assert(hbs2d.dim==2)
assert(hbs2d.numElements==55)
assert(hbs2d.size==55)
assert(hbs2d.treeSize==13) 
assert(hbs2d.treeLeafSize==7) 
para2 = hbs2d.support;
assert(isequal(para2,[0 1; 0 1]))
assert(hbs2d.maxLevel==3)
assert(hbs3d.dim==3)
assert(hbs3d.numElements==15)
assert(hbs3d.size==34)
assert(hbs3d.treeSize==7) 
assert(hbs3d.treeLeafSize==4) 
para3 = hbs3d.support;
assert(isequal(para3,[0 1; 0 1; 0 1]))
assert(hbs3d.maxLevel==2)
fprintf('Test on accesssors: passed.\n\n')

%% TEST OTHER METHODS
% Print trees
hbs1.treePrintLeaves;
fprintf('\nCheck that the output is:\nLeaf node (0), (16384). level=0  \n\n')
hbs2.treePrintLeaves;
fprintf('\nCheck that the output is:\nLeaf node (0 0), (16384 16384). level=0 \n\n')
hbs3.treePrintLeaves;
fprintf('\nCheck that the output is:\nLeaf node (0 0 0), (16384 16384 16384). level=0  \n\n')

% Check degree (dir)
assert(hbs1d.degree(1)==1)
assert(hbs2d.degree(1)==2)
assert(hbs3d.degree(1)==1)

% Check active functions and evaluations at pts
pt1 = 0.4;
assert(isequal(hbs1d.active(pt1),[2;3]))
assert(ismembertol(hbs1d.eval(pt1)',[0.4 0.6],'ByRows',true))
pt = [0.5;0.5];
assert(isequal(hbs2d.active(pt),[15;16;17;21;22;23;27;28;29;37;38;39;50;51;54;55]))
assert(isequal(hbs2d.eval(pt),[0;0;0;0;0;0;0;0;0;0;0;0;0.25;0.25;0.25;0.25]))
pt3 = [0.7;0.7;0.7];
assert(isequal(hbs3d.active(pt3),[13;14;16;17;22;23;25;26]))
assert(ismembertol(hbs3d.eval(pt3)',[0.216 0.144 0.144 0.096 0.144 0.096 0.096 0.064],...
                   'ByRows',true))

% Check evaluations at a single active basis function
assert(ismembertol(hbs1d.evalSingle(2,pt1),0.4));
assert(ismembertol(hbs2d.evalSingle(55,pt),0.25));
assert(ismembertol(hbs3d.evalSingle(13,pt3),0.216));

% Save to output file
hbs1d.save('basis1d');
hbs2d.save('basis2d');
hbs3d.save('basis3d');
fprintf('Check the xml output files basis1d, basis2d and basis3d.\n')
% basis1d.xml
%
% <xml>
%  <Basis type="THBSplineBasis" id="0">
%   <Basis type="BSplineBasis">
%    <KnotVector degree="1">0 0 0.25 0.5 0.75 1 1 </KnotVector>
%   </Basis>
%  </Basis>
% </xml>
%
% basis2d.xml
%
% <xml>
%  <Basis type="THBSplineBasis2" id="0">
%   <Basis type="TensorBSplineBasis2">
%    <Basis type="BSplineBasis" index="0">
%     <KnotVector degree="2">0 0 0 0.25 0.5 0.75 1 1 1 </KnotVector>
%    </Basis>
%    <Basis type="BSplineBasis" index="1">
%     <KnotVector degree="2">0 0 0 0.25 0.5 0.75 1 1 1 </KnotVector>
%    </Basis>
%   </Basis>
%   <box level="1">5 2 6 6 
% </box>
%   <box level="1">2 5 5 6 
% </box>
%   <box level="2">4 4 10 10 
% </box>
%  </Basis>
% </xml>
%
%
% basis3d.xml
%
% <xml>
%  <Basis type="THBSplineBasis3" id="0">
%   <Basis type="TensorBSplineBasis3">
%    <Basis type="BSplineBasis" index="0">
%     <KnotVector degree="1">0 0 0.5 1 1 </KnotVector>
%    </Basis>
%    <Basis type="BSplineBasis" index="1">
%     <KnotVector degree="1">0 0 0.5 1 1 </KnotVector>
%    </Basis>
%    <Basis type="BSplineBasis" index="2">
%     <KnotVector degree="1">0 0 0.5 1 1 </KnotVector>
%    </Basis>
%   </Basis>
%   <box level="1">0 0 0 2 2 2 
% </box>
%  </Basis>
% </xml>
% 

% Check knot vectors
assert(isequal(hbs1d.knots(1,1),[0,linspace(0,1,5),1]));
assert(isequal(hbs2d.knots(1,1),[0,0,linspace(0,1,5),1,1]));
assert(isequal(hbs3d.knots(2,3),[0,linspace(0,1,5),1]));
% Check number of breaks
assert(hbs1d.numBreaks(1,1)==5)
assert(hbs2d.numBreaks(1,1)==5)
assert(hbs3d.numBreaks(2,3)==5)

% Check getBoxes
[b1,~,~] = hbs1d.getBoxes;
assert(b1==0);
[~,b2,~] = hbs2d.getBoxes;
assert(isequal(b2(3,:),[12 12]));
[~,~,lev] = hbs3d.getBoxes;
assert(isequal(lev,[1;1;1;2]));

% Uniformly refine the basis
hbs1d.uniformRefine(1,1);
assert(isequal(hbs1d.knots(1,1),[0,linspace(0,1,9),1]));
hbs2d.uniformRefine(1,1);
assert(isequal(hbs2d.knots(1,1),[0,0,linspace(0,1,9),1,1]));
hbs3d.uniformRefine(1,1);
assert(isequal(hbs3d.knots(2,3),[0,linspace(0,1,9),1]));

% Uniformly refine the basis with update of coefficients
hbs2_copy = gsTHBSplineBasis(hbs2,2);
clear coefs
[coefs(:,:,1),coefs(:,:,2)] = ndgrid(1:4,1:4);
coefs = reshape(coefs,16,2);
coefs2 = hbs2_copy.uniformRefine_withCoefs(coefs,1,1);
clear coefs3
[coefs3(:,:,1),coefs3(:,:,2)] = ndgrid([1 1.5 2.25 2.75 3.5 4],[1 1.5 2.25 ...
    2.75 3.5 4]);
coefs3 = reshape(coefs3,36,2);
assert(isequal(coefs2,coefs3))

% Refine the basis by defining boxes
boxes = [3,1,3];
hbs1d.refineElements(boxes);
assert(hbs1d.numElements==10);
hbs2_copy = gsTHBSplineBasis(hbs2,2);
boxes = [2,1,3,3,5];
hbs2.refineElements(boxes);
assert(hbs2.numElements==7);
boxes = [3,1,1,1,3,3,3];
hbs3d.refineElements(boxes);
assert(hbs3d.numElements==127);

% Refine the basis by defining boxes with update of coefficients
boxes = [2 1 1 3 3];
coefs4 = hbs2_copy.refineElements_withCoefs(coefs, boxes);
coefs =[ 2 1; 3 1; 4 1; 1 2; 2 2; 3 2; 4 2; 1 3; 2 3; 3 3; 4 3; 1 4; 2 4; ...
    3 4; 4 4; 1 1; 1.5 1; 1 1.5; 1.5 1.5];
assert(isequal(coefs,coefs4));

% Slice the basis
slice = hbs2_copy.basisSlice(1, 0.5);
assert(slice.numElements == 3)
slice = hbs3_copy.basisSlice(1, 0.5);
assert(slice.numElements == 7)

fprintf('All tests: passed.\n')
