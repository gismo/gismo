% This MATLAB script tests the MEX interface of the gsTHBSplineBasis class.
% Author: O. Chanon, P. Noertoft

%% TEST CONSTRUCTORS
% Construct a truncated hierarchical basis by reading the specified file
filename = join([filedata,'thbbasis/simple.xml']);
fprintf('Reading THB spline basis from file: %s.\n',filename)
hbs = gsTHBSplineBasis(filename);

% Construct a truncated hierarchical basis from its knot vector in a carthesian product way.
knots_2build = {[0,0,0,1,2,2,2],[4,4,4,5,6,6,6]};
fprintf('Reading THB spline basis from knot vector, cartesian product way.\n')
hbs2 = gsTHBSplineBasis(knots_2build);

%% TEST ACCESSORS
assert(hbs.dim==2)
assert(hbs2.dim==2)
assert(hbs.numElements==55)
assert(hbs2.numElements==4)
assert(hbs.size==55)
assert(hbs2.size==16)
assert(hbs.treeSize==13) 
assert(hbs2.treeSize==1)
assert(hbs.treeLeafSize==7) 
assert(hbs2.treeLeafSize==1)
para = hbs.support;
assert(isequal(para,[0 1; 0 1]))
para2 = hbs2.support;
assert(isequal(para2,[0 2; 4 6]))
assert(hbs.maxLevel==3) % nlevels in GeoPDEs
assert(hbs2.maxLevel==1)
fprintf('Test on accesssors: passed.\n')

%% TEST OTHER METHODS
% Print trees
fprintf('\n Tree 1:\n')
hbs.treePrintLeaves;
fprintf('\nCheck that the output is:\nLeaf node (0 0), ( 8192 32768). level=0 \nLeaf node (8192    0), (24576  8192). level=0 \nLeaf node (8192 8192), (20480 20480). level=2 \nLeaf node ( 8192 20480), (20480 24576). level=1 \nLeaf node (20480  8192), (24576 24576). level=1 \nLeaf node ( 8192 24576), (24576 32768). level=0 \nLeaf node (24576     0), (32768 32768). level=0\n')
fprintf('\n Tree 2:\n')
hbs2.treePrintLeaves;
fprintf('\nCheck that the output is:\nLeaf node (0 0), (16384 16384). level=0 \n\n')

% Check degree
assert(hbs.degree(1)==2) % order = degree-1 in geopdes
assert(hbs2.degree(hbs2.dim)==2)

% Check active functions and evaluations at pts
pt = [0.5;0.5];
assert(isequal(hbs.active(pt),[15;16;17;21;22;23;27;28;29;37;38;39;50;51;54;55]))
assert(isequal(hbs.eval(pt),[0;0;0;0;0;0;0;0;0;0;0;0;0.25;0.25;0.25;0.25]))

pt2 = [1;5];
ev2  = hbs2.eval(pt2);
assert(isequal(hbs2.active(pt2),[6;7;8;10;11;12;14;15;16]))
assert(isequal(hbs2.eval(pt2),[0.25;0.25;0;0.25;0.25;0;0;0;0]))

% Check evaluations at a single active basis function
assert(hbs.evalSingle(55,pt)==0.25);
assert(hbs2.evalSingle(6,pt2)==0.25);

%%% assert(hbs.elementIndex([0;0])==1) TODO

% Save to output file
hbs.save('output1');
hbs2.save('output2');
fprintf('Check the output files output1.xml and output2.xml\n\n')
% output1.xml
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
% output2.xml
%
% <xml>
%  <Basis type="THBSplineBasis2" id="0">
%   <Basis type="TensorBSplineBasis2">
%    <Basis type="BSplineBasis" index="0">
%     <KnotVector degree="2">0 0 0 1 2 2 2 </KnotVector>
%    </Basis>
%    <Basis type="BSplineBasis" index="1">
%     <KnotVector degree="2">4 4 4 5 6 6 6 </KnotVector>
%    </Basis>
%   </Basis>
%  </Basis>
% </xml>

% Check knot vectors
assert(isequal(hbs.knots(1,1),[0,0,linspace(0,1,5),1,1]));
assert(isequal(hbs2.knots(hbs2.maxLevel,hbs2.dim),knots_2build{2}));
% Check number of breaks
assert(hbs.numBreaks(1,1)==5)
assert(hbs2.numBreaks(hbs2.maxLevel,hbs2.dim)==length(unique(knots_2build{2})))

% Check getBoxes
[b1,b2,lev] = hbs.getBoxes;
assert(isequal(b1,[12 0; 4 12; 10 4; 4 10; 4 4; 4 0; 0 0]));
assert(isequal(b2,[16 16; 12 16; 12 12; 10 12; 10 10; 12 4; 4 16]));
assert(isequal(lev,[1;1;2;2;3;1;1]));

% Uniformily refine the basis
hbs.uniformRefine(1,1);
assert(isequal(hbs.knots(1,1),[0,0,linspace(0,1,9),1,1]));

% Refine the basis by defining boxes
boxes = [2,1,3,3,5];
hbs2.refineElements(boxes);
assert(hbs2.numElements==7);

fprintf('All tests: passed.\n')

