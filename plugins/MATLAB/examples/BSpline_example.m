clear all;
close all;
clc;

addpath(genpath('../src/'));
addpath(genpath('../../../optional/gsCInterface/'));

KV = gsKnotVector([0,0,0,1,1,1]);
BB = gsBSplineBasis(KV);
BB.uniformRefine();

coefs = rand(BB.size(),2);
B = gsBSpline(BB,coefs);

points1D = linspace(0,1,100);
ev = transpose(B.eval(points1D));

[dist,xi] = B.closest(coefs(3,:)');
evxi = transpose(B.eval(xi));

normal = B.normal(xi);

figure;
hold on;
axis equal;
xlim([0,1]);
ylim([0,1]);
plot(ev(:,1),ev(:,2))
scatter(coefs(:,1),coefs(:,2))
scatter(evxi(1,1),evxi(1,2))
quiver(evxi(1,1),evxi(1,2),normal(1,1),normal(2,1))

