% G+Smo library definition for Matlab

% Setup platform specific variables
if isunix
    % Set up compiler - g++
    mex -setup:g++;

    srcpath = "/user/amantzaf/home/Workspace/gismo/stable"
    binpath = "/user/amantzaf/home/Workspace/gismo/stable/build_dbg"
    
    definedMacros = ["EIGEN_MATRIXBASE_PLUGIN=<gsMatrix/gsMatrixAddons.h>" "EIGEN_PLAINOBJECTBASE_PLUGIN=<gsMatrix/gsPlainObjectBaseAddons.h>" "EIGEN_DEFAULT_DENSE_INDEX_TYPE=index_t" "EIGEN_DEFAULT_SPARSE_INDEX_TYPE=index_t" "EIGEN_DEFAULT_TO_COL_MAJOR" "EIGEN_NO_STATIC_ASSERT" "GISMO_BUILD_LIB" "_XOPEN_SOURCE"];
    additionalCompilerFlags = ["-std=c++14"];
    undefinedMacros = ["_XOPEN_SOURCE" "__USE_XOPEN2K8"];
    additionalLinkerFlags = "LFLAGgismo";
    
    interfaceGenerationFiles = [...
    % fullfile(binpath,"gsCore/gsConfig.h"), ...
    % fullfile(binpath,"gsCore/gsConfigExt.h"), ...
    % fullfile(binpath,"gsCore/gsExport.h"), ...
    % fullfile(srcpath,"src/gsCore/gsLinearAlgebra.h"), ...
    % fullfile(srcpath,"src/gsCore/gsLinearAlgebra.h"), ...
        fullfile(srcpath,"src/gsCore/gsMath.h"), ...
        fullfile(srcpath,"src/gsUtils/gsCombinatorics.h") ...
       ];

    includePath = [binpath, fullfile(srcpath,"src"), fullfile(srcpath,"external"), fullfile(srcpath,"optional")];
    libraries = fullfile(binpath,"lib/libgismo.so");
    outputFolder = fullfile(srcpath,"plugins/gsMex/lib");
    overwriteExistingDefinitionFiles = true;
else
    error('Live task not configured for current platform')
end

% Check if the library definition exists
if ~overwriteExistingDefinitionFiles && isfile(fullfile(outputFolder,'definegismo.m')) && evalin('base','exist(''libraryDefinitionFromTask'',''var'') == 0')
    disp('Click this link to restore missing workspace variable for existing definition file definegismo.m.');
    disp('<a href="matlab: internal.matlab.desktop.commandwindow.executeCommandForUser(''addpath(outputFolder);libraryDefinitionFromTask = feval(''''definegismo'''');'')">Restore library definition workspace variable</a>');
end

% Generate definition file for C++ library
clibgen.generateLibraryDefinition(interfaceGenerationFiles, ...
    "IncludePath",includePath, ...
    "Libraries",libraries, ...
    "OutputFolder",outputFolder, ...
    "PackageName","gismo", ...
    "OverwriteExistingDefinitionFiles",overwriteExistingDefinitionFiles, ...
    "AdditionalCompilerFlags",additionalCompilerFlags, ...
    "Verbose",true);
%    "UndefinedMacros",undefinedMacros, ...

% Create the library definition object
addpath(outputFolder);
libraryDefinitionFromTask = feval("definegismo");

% Show available constructs
summary(libraryDefinitionFromTask);

% Clear temporary variables
clear('interfaceGenerationFiles','includePath','libraries','outputFolder','overwriteExistingDefinitionFiles');
