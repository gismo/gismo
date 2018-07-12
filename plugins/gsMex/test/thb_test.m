% This MATLAB script tests the MEX interface of the gsTHBSplineBasis class.
% Author: O. Chanon, P. Noertoft

%% TEST CONSTRUCTORS
% Construct a truncated hierarchical basis by reading the specified file
filename = join([filedata,'thbbasis/simple.xml']);
fprintf('Reading THB spline basis from file: %s.\n',filename)
hbs = gsTHBSplineBasis(filename);

% Construct a truncated hierarchical basis from its knot vectors.
knots_2build = {[0,0,0,1,2,2,2],[4,4,5,5.5,6,6]};
fprintf('Reading THB spline basis from knot vector.\n')
hbs2 = gsTHBSplineBasis(knots_2build);

%% TEST ACCESSORS
fprintf('Dimension of the parameter space 1: %d\n',hbs.dim);        % rdim in geopdes
fprintf('Dimension of the parameter space 2: %d\n',hbs2.dim);
fprintf('Number of elements 1: %d\n',hbs.numElements);              % nel in geopdes
fprintf('Number of elements 2: %d\n',hbs2.numElements); 
fprintf('Nb of degrees of freedom of the basis 1: %d\n',hbs.size);  % ndof in geopdes
fprintf('Nb of degrees of freedom of the basis 2: %d\n',hbs2.size);
fprintf('The tree 1 has %d nodes.\n',hbs.treeSize);
fprintf('The tree 2 has %d nodes.\n',hbs2.treeSize);
fprintf('Num. of leaves 1: %d.\n',hbs.treeLeafSize);
fprintf('Num. of leaves 2: %d.\n',hbs.treeLeafSize);
para = hbs.support;
fprintf('The parameter range 1 is: [%f %f %f %f]\n', para);
para2 = hbs2.support;
fprintf('The parameter range 2 is: [%f %f %f %f]\n', para2);
fprintf('Number of levels 1: %d\n',hbs.maxLevel);                    % nlevels in geopdes
fprintf('Number of levels 2: %d\n',hbs2.maxLevel);

%% TEST OTHER METHODS
% Print trees
fprintf('Tree 1:\n')
hbs.treePrintLeaves;
fprintf('Tree 2:\n')
hbs2.treePrintLeaves;

% Print degree
fprintf('Degree in the first direction 1: %d\n',hbs.degree(1));     % order = degree-1 in geopdes
fprintf('Degree in the last direction 2: %d\n',hbs2.degree(hbs2.dim));

% Print evaluations at pts
pts = uniformPointGrid(para(1:2),para(3:4),5);
ev  = hbs.eval(pts');
fprintf('Evaluation 1\n')
disp(ev)

pts2 = uniformPointGrid(para2(1:2),para2(3:4),5);
ev2  = hbs2.eval(pts2);
fprintf('Evaluation 2\n')
disp(ev2)

% Print evaluations at a single active basis function
pt = [0.7;0.7];
evs = hbs.evalSingle(16,pt);
fprintf('Evaluation 1 at basis function 16\n')
disp(evs)

pt2 = [0.7;5.7];
evs2 = hbs2.evalSingle(16,pt2);
fprintf('Evaluation 2 at basis function 16\n')
disp(evs2)

% Save to output file
hbs.save('output1');
hbs2.save('output2');
fprintf('Saved.\n')

% Print knot vectors
knots = hbs.knots(1,1);
fprintf('Knots 1, level 1, direction 1\n')
disp(knots)

knots2 = hbs2.knots(hbs2.maxLevel,hbs2.dim);
fprintf('Knots 2 level max, direction max\n')
disp(knots)

% Print active functions
act = hbs.active(pt);
fprintf('Active functions 1 at point:\n')
disp(act)

act2 = hbs2.active(pt2);
fprintf('Active functions 2 at point:\n')
disp(act2)



