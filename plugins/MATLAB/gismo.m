clear all;
close all;
clc;


if (libisloaded('libgismo'))
unloadlibrary libgismo
end

path_to_main = '/home/hverhelst/Documents/code/gismo_mex';
path_to_gismo = [path_to_main,'/build/lib/libgismo'];
path_to_Cinterface = '/optional/gsCInterface/src/Cgismo.h';
% path_to_Minterface = './optional/gsCInterface/MATLAB/';
path_to_gsCore = '/home/hverhelst/Documents/code/gismo_mex/build';
path_to_gsCInterface = '/optional/gsCInterface/src/';

addpath(genpath(path_to_gsCInterface));

% addpath(path_to_Minterface)

% The following will create the library 'libgismo'
loadlibrary(path_to_gismo,...
            path_to_Cinterface,... Path where to find the C interface
            'addheader','gsCTypes.h',...
            'addheader','gsCMatrix.h',...
            'addheader','gsCMatrixInt.h',...
            'addheader','gsCVector.h',...
            'addheader','gsCVectorInt.h',...
            'addheader','gsCKnotVector.h',...
            'addheader','gsCFunctionSet.h',...
            'addheader','gsCBasis.h',...
            'addheader','gsCGeometry.h',...
            'addheader','gsCMultiPatch.h',...
            'addheader','gsCMultiPatch.h',...
            'addheader','gsCReadFile.h',...
            'includepath', ...
                path_to_gsCore... Path where to find gsCore
                );

libfunctions('libgismo')



% mat = calllib('libgismo','gsMatrix_create');
% calllib('libgismo','gsMatrix_delete',mat);

% A=zeros(10);
% Aptr=libpointer("doublePtr",A);

% size(A)

% mat = EigenMatrix(size(A,1),size(A,2),Aptr);

% % get(Aptr)

% Gismo.rows(mat)
% Gismo.cols(mat)
% ptr = Gismo.data(mat);

% filename = 'filedata/surfaces/simple.xml';
% g = Geometry(filename)

% Gismo.domainDim(g)
% Gismo.targetDim(g)

% pts = [0,1; 0,1];
% Gismo.eval(g,pts)

% % ptr.DataType
% % ptr.Value
