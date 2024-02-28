 % This MATLAB script tests the MEX interface of the gsTensorBSplineBasis class.
% Author: H.M.Verhelst, A.Mantzaflaris

% Construct a truncated hierarchical basis by reading the specified file
% kv = [0.0,0.0,1.0,1.0]
kv = gsKnotVector([0.0,0.0,1.0,1.0])
tbasis2 = gsTensorBSplineBasis(kv,kv)

% Print sizes
fprintf('Size of the basis %d\n',tbasis2.size());
fprintf('Dimension of the parameter space: %d\n',tbasis2.dim());

% Print support
para = tbasis2.support();
fprintf('The parameter range is: [%f %f %f %f]\n', para);

% Print points
pts = uniformPointGrid(para(1:2),para(3:4),5);
fprintf('pts\n')
disp(pts)

% Print evaluations
ev  = tbasis2.eval(pts);
fprintf('ev\n')
disp(ev)

% Print sum of evaluations
fprintf('sum(ev)\n')
disp(sum(ev,1))

% Print active functions
act = tbasis2.active(pts);
fprintf('act\n')
disp(act)

% Print single evaluation
evs = tbasis2.evalSingle(1,pts);
fprintf('evs\n')
disp(evs)

% % Save to output file
% tbasis2.save('output');

% % Print knot vector of level 1, direction 1
% knots = tbasis2.knots(1);
% fprintf('knots level 1, direction 1\n')
% disp(knots)

tbasis3 = gsTensorBSplineBasis(kv,kv,kv)

% Print sizes
fprintf('Size of the basis %d\n',tbasis3.size());
fprintf('Dimension of the parameter space: %d\n',tbasis3.dim());

% Print support
para = tbasis3.support();
fprintf('The parameter range is: [%f %f %f %f %f %f]\n', para);

% Print points
pts = uniformPointGrid(para(1:3),para(4:6),5);
fprintf('pts\n')
disp(pts)

% Print evaluations
ev  = tbasis3.eval(pts);
fprintf('ev\n')
disp(ev)

% Print sum of evaluations
fprintf('sum(ev)\n')
disp(sum(ev,1))

% Print active functions
act = tbasis3.active(pts);
fprintf('act\n')
disp(act)

% Print single evaluation
evs = tbasis3.evalSingle(1,pts);
fprintf('evs\n')
disp(evs)

% % Save to output file
% tbasis3.save('output');

% Print knot vector of direction 1
knots = tbasis3.knots(1);
fprintf('knots direction 1\n')
disp(knots)