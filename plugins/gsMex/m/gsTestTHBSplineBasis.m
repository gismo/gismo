% This MATLAB script tests the MEX interface of the gsTHBSplineBasis class.

% Author: Peter Noertoft

% Initiliazations
clear all classes
close all

% Update path
addpath('stable/plugins/gsMex/');

% Construct a truncated hierarchical basis by reading the specified file
filename = '../stable/filedata/basis_thbs.xml';
fprintf('Reading THB spline basis from file: %s.\n',filename)
hbs = gsTHBSplineBasis(filename);

% Print sizes
fprintf('Size of the basis %d\n',hbs.size());
fprintf('The tree has %d nodes.\n',hbs.treeSize());

% Print support
para = hbs.support();
fprintf('The parameter range is: [%f %f %f %f]\n', para);

% Print tree
fprintf('Num. of leaves: %d.\n',hbs.treeLeafSize());
hbs.treePrintLeaves();

% Print points
pts = gsUniformPointGrid(para(1:2),para(3:4),11);
fprintf('pts\n')
disp(pts)

% Print evaluations
ev  = hbs.eval(pts);
fprintf('ev\n')
disp(ev)

% Print sum of evaluations
fprintf('sum(ev)\n')
disp(sum(ev,1))

% Print active functions
act = hbs.active(pts);
fprintf('act\n')
disp(act)

% Print single evaluation
evs = hbs.evalSingle(0,pts);
fprintf('evs\n')
disp(evs)
