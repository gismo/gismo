%% About definegismo.m
% This file defines the MATLAB interface to the library |gismo|.
%
% Commented sections represent C++ functionality that MATLAB cannot automatically define. To include
% functionality, uncomment a section and provide values for <SHAPE>, <DIRECTION>, etc. For more
% information, see helpview(fullfile(docroot,'matlab','helptargets.map'),'cpp_define_interface') to "Define MATLAB Interface for C++ Library".



%% Setup
% Do not edit this setup section.
function libDef = definegismo()
libDef = clibgen.LibraryDefinition("gismoData.xml");

%% OutputFolder and Libraries 
libDef.OutputFolder = "/home/amantzaf/Workspace/gismo/stable/plugins/gsMex/lib";
libDef.Libraries = "/user/amantzaf/home/Workspace/gismo/stable/build_dbg/lib/libgismo.so";

%% C++ class |gismo::binomialT<0, 0>| with MATLAB name |clib.gismo.gismo.binomialT_0_0_| 
binomialT_0_0_Definition = addClass(libDef, "gismo::binomialT<0, 0>", "MATLABName", "clib.gismo.gismo.binomialT_0_0_", ...
    "Description", "clib.gismo.gismo.binomialT_0_0_    Representation of C++ class gismo::binomialT<0, 0>."); % Modify help description values as needed.

%% C++ class constructor for C++ class |gismo::binomialT<0, 0>| 
% C++ Signature: gismo::binomialT<0, 0>::binomialT(gismo::binomialT<0, 0> const & input1)

binomialT_0_0_Constructor1Definition = addConstructor(binomialT_0_0_Definition, ...
    "gismo::binomialT<0, 0>::binomialT(gismo::binomialT<0, 0> const & input1)", ...
    "Description", "clib.gismo.gismo.binomialT_0_0_ Constructor of C++ class gismo::binomialT<0, 0>."); % Modify help description values as needed.
defineArgument(binomialT_0_0_Constructor1Definition, "input1", "clib.gismo.gismo.binomialT_0_0_", "input");
validate(binomialT_0_0_Constructor1Definition);

%% C++ class constructor for C++ class |gismo::binomialT<0, 0>| 
% C++ Signature: gismo::binomialT<0, 0>::binomialT()

binomialT_0_0_Constructor2Definition = addConstructor(binomialT_0_0_Definition, ...
    "gismo::binomialT<0, 0>::binomialT()", ...
    "Description", "clib.gismo.gismo.binomialT_0_0_ Constructor of C++ class gismo::binomialT<0, 0>."); % Modify help description values as needed.
validate(binomialT_0_0_Constructor2Definition);

%% C++ function |gismo::factorial| with MATLAB name |clib.gismo.gismo.factorial|
% C++ Signature: unsigned int gismo::factorial(unsigned int n)

factorialDefinition = addFunction(libDef, ...
    "unsigned int gismo::factorial(unsigned int n)", ...
    "MATLABName", "clib.gismo.gismo.factorial", ...
    "Description", "clib.gismo.gismo.factorial Representation of C++ function gismo::factorial."); % Modify help description values as needed.
defineArgument(factorialDefinition, "n", "uint32");
defineOutput(factorialDefinition, "RetVal", "uint32");
validate(factorialDefinition);

%% C++ function |gismo::numCubeElements| with MATLAB name |clib.gismo.gismo.numCubeElements|
% C++ Signature: int gismo::numCubeElements(int const k,int const d)

numCubeElementsDefinition = addFunction(libDef, ...
    "int gismo::numCubeElements(int const k,int const d)", ...
    "MATLABName", "clib.gismo.gismo.numCubeElements", ...
    "Description", "clib.gismo.gismo.numCubeElements Representation of C++ function gismo::numCubeElements." + newline + ...
    "Returns the number of elements (faces) of dimension \a k " + newline + ...
    " of a \a d-cube" + newline + ...
    " \ingroup combinatorics"); % Modify help description values as needed.
defineArgument(numCubeElementsDefinition, "k", "int32");
defineArgument(numCubeElementsDefinition, "d", "int32");
defineOutput(numCubeElementsDefinition, "RetVal", "int32");
validate(numCubeElementsDefinition);

%% C++ function template instantiation |gismo::binomial<int>| with MATLAB name |clib.gismo.gismo.binomial|
% C++ Signature: int gismo::binomial<int>(int const n,int const r)

binomialDefinition = addFunction(libDef, ...
    "int gismo::binomial<int>(int const n,int const r)", ...
    "MATLABName", "clib.gismo.gismo.binomial", ...
    "TemplateUniqueName", "clib.gismo.gismo.binomial_int_", ...
    "Description", "clib.gismo.gismo.binomial Representation of C++ function gismo::binomial."); % Modify help description values as needed.
defineArgument(binomialDefinition, "n", "int32");
defineArgument(binomialDefinition, "r", "int32");
defineOutput(binomialDefinition, "RetVal", "int32");
validate(binomialDefinition);

%% C++ function |gismo::numCompositions| with MATLAB name |clib.gismo.gismo.numCompositions|
% C++ Signature: unsigned int gismo::numCompositions(int sum,int dim)

numCompositionsDefinition = addFunction(libDef, ...
    "unsigned int gismo::numCompositions(int sum,int dim)", ...
    "MATLABName", "clib.gismo.gismo.numCompositions", ...
    "Description", "clib.gismo.gismo.numCompositions Representation of C++ function gismo::numCompositions." + newline + ...
    "Number of compositions of \a sum into \a dim integers" + newline + ...
    " \ingroup combinatorics"); % Modify help description values as needed.
defineArgument(numCompositionsDefinition, "sum", "int32");
defineArgument(numCompositionsDefinition, "dim", "int32");
defineOutput(numCompositionsDefinition, "RetVal", "uint32");
validate(numCompositionsDefinition);

%% Validate the library definition
validate(libDef);

end
