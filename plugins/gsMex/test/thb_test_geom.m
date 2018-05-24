% This MATLAB script tests the MEX interface of the gsTHBSplineBasis class.
% Author: Peter Noertoft

% Construct a truncated hierarchical basis by reading the specified file
filename = join([filedata,'surfaces/thbs_face_3levels.xml']);
fprintf('Reading THB spline from file: %s.\n',filename)
hbs = gsTHBSpline(filename);

% Print support
para = hbs.support();
fprintf('The parameter range is: [%f %f %f %f]\n', para);

% Print points
pts = uniformPointGrid(para(1:2),para(3:4),1000);

% Print evaluations
ev  = hbs.eval(pts);
fprintf('ev\n')
disp(size(ev))

scatter3(ev(1,:), ev(2,:), ev(3,:))