
clear all;
close all;
clc;

addpath(genpath('../src/'));
addpath(genpath('../../../optional/gsCInterface/'));

KV = gsKnotVector([0,0,0,1,1,1]);
TBB = gsTensorBSplineBasis(KV,KV);
TBB.uniformRefine();
THB = gsTHBSplineBasis(TBB);

% THB.refine([0,0.5;0.0,0.5]);
THB.refineElements([1,0,0,2,2]); % equivalent

disp(THB)

coefs = rand(TBB.size(),3);
TB = gsTensorBSpline(TBB,coefs);

points1D = linspace(0,1,100);
[xi,eta] = meshgrid(points1D,points1D);
N = length(points1D);
Nflat = length(points1D)*length(points1D);
points2D = [reshape(xi,[1,Nflat]); reshape(eta,[1,Nflat])];

ev = TB.eval(points2D);

X = reshape(ev(1,:),[N,N]);
Y = reshape(ev(2,:),[N,N]);
Z = reshape(ev(3,:),[N,N]);

figure;
hold on;
axis equal;
xlim([0,1]);
ylim([0,1]);
zlim([0,1]);
surf(X,Y,Z)
scatter3(coefs(:,1),coefs(:,2),coefs(:,3))

