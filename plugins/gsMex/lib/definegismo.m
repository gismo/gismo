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

%% C++ class |gismo::gsDomain<double>| with MATLAB name |clib.gismo.gismo.gsDomain_double_| 
gsDomain_double_Definition = addClass(libDef, "gismo::gsDomain<double>", "MATLABName", "clib.gismo.gismo.gsDomain_double_", ...
    "Description", "clib.gismo.gismo.gsDomain_double_    Representation of C++ class gismo::gsDomain<double>."); % Modify help description values as needed.

%% C++ class method |dim| for C++ class |gismo::gsDomain<double>| 
% C++ Signature: int gismo::gsDomain<double>::dim() const

dimDefinition = addMethod(gsDomain_double_Definition, ...
    "int gismo::gsDomain<double>::dim() const", ...
    "MATLABName", "dim", ...
    "Description", "dim Method of C++ class gismo::gsDomain<double>."); % Modify help description values as needed.
defineOutput(dimDefinition, "RetVal", "int32");
validate(dimDefinition);

%% C++ class method |minMeshSize| for C++ class |gismo::gsDomain<double>| 
% C++ Signature: double gismo::gsDomain<double>::minMeshSize()

minMeshSizeDefinition = addMethod(gsDomain_double_Definition, ...
    "double gismo::gsDomain<double>::minMeshSize()", ...
    "MATLABName", "minMeshSize", ...
    "Description", "minMeshSize Method of C++ class gismo::gsDomain<double>."); % Modify help description values as needed.
defineOutput(minMeshSizeDefinition, "RetVal", "double");
validate(minMeshSizeDefinition);

%% C++ class method |breaks| for C++ class |gismo::gsDomain<double>| 
% C++ Signature: std::vector<double, std::allocator<double>> gismo::gsDomain<double>::breaks() const

breaksDefinition = addMethod(gsDomain_double_Definition, ...
    "std::vector<double, std::allocator<double>> gismo::gsDomain<double>::breaks() const", ...
    "MATLABName", "breaks", ...
    "Description", "breaks Method of C++ class gismo::gsDomain<double>."); % Modify help description values as needed.
defineOutput(breaksDefinition, "RetVal", "clib.array.gismo.Double");
validate(breaksDefinition);

%% C++ class method |clone| for C++ class |gismo::gsDomain<double>| 
% C++ Signature: gismo::gsDomain<double> * gismo::gsDomain<double>::clone() const

%cloneDefinition = addMethod(gsDomain_double_Definition, ...
%    "gismo::gsDomain<double> * gismo::gsDomain<double>::clone() const", ...
%    "MATLABName", "clone", ...
%    "Description", "clone Method of C++ class gismo::gsDomain<double>."); % Modify help description values as needed.
%defineOutput(cloneDefinition, "RetVal", "clib.gismo.gismo.gsDomain_double_", <SHAPE>);
%validate(cloneDefinition);

%% C++ class method |unique| for C++ class |gismo::gsDomain<double>| 
% C++ Signature: std::vector<double, std::allocator<double>> gismo::gsDomain<double>::unique() const

uniqueDefinition = addMethod(gsDomain_double_Definition, ...
    "std::vector<double, std::allocator<double>> gismo::gsDomain<double>::unique() const", ...
    "MATLABName", "unique", ...
    "Description", "unique Method of C++ class gismo::gsDomain<double>."); % Modify help description values as needed.
defineOutput(uniqueDefinition, "RetVal", "clib.array.gismo.Double");
validate(uniqueDefinition);

%% C++ class method |merge| for C++ class |gismo::gsDomain<double>| 
% C++ Signature: void gismo::gsDomain<double>::merge(gismo::gsDomain<double> * other)

%mergeDefinition = addMethod(gsDomain_double_Definition, ...
%    "void gismo::gsDomain<double>::merge(gismo::gsDomain<double> * other)", ...
%    "MATLABName", "merge", ...
%    "Description", "merge Method of C++ class gismo::gsDomain<double>."); % Modify help description values as needed.
%defineArgument(mergeDefinition, "other", "clib.gismo.gismo.gsDomain_double_", "input", <SHAPE>); % <MLTYPE> can be "clib.gismo.gismo.gsDomain_double_", or "clib.array.gismo.gismo.gsDomain_double_"
%validate(mergeDefinition);

%% C++ class |gismo::gsKnotVector<double>| with MATLAB name |clib.gismo.gismo.gsKnotVector_double_| 
gsKnotVector_double_Definition = addClass(libDef, "gismo::gsKnotVector<double>", "MATLABName", "clib.gismo.gismo.gsKnotVector_double_", ...
    "Description", "clib.gismo.gismo.gsKnotVector_double_    Representation of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.

%% C++ class constructor for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::gsKnotVector()

gsKnotVector_double_Constructor1Definition = addConstructor(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::gsKnotVector()", ...
    "Description", "clib.gismo.gismo.gsKnotVector_double_ Constructor of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
validate(gsKnotVector_double_Constructor1Definition);

%% C++ class constructor for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::gsKnotVector(gismo::gsKnotVector<double>::knotContainer knots,int degree)

gsKnotVector_double_Constructor2Definition = addConstructor(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::gsKnotVector(gismo::gsKnotVector<double>::knotContainer knots,int degree)", ...
    "Description", "clib.gismo.gismo.gsKnotVector_double_ Constructor of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(gsKnotVector_double_Constructor2Definition, "knots", "clib.array.gismo.Double");
defineArgument(gsKnotVector_double_Constructor2Definition, "degree", "int32");
validate(gsKnotVector_double_Constructor2Definition);

%% C++ class method |swap| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::swap(gismo::gsKnotVector<double> & other)

swapDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::swap(gismo::gsKnotVector<double> & other)", ...
    "MATLABName", "swap", ...
    "Description", "swap Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(swapDefinition, "other", "clib.gismo.gismo.gsKnotVector_double_", "input");
validate(swapDefinition);

%% C++ class method |clone| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double> * gismo::gsKnotVector<double>::clone() const

%cloneDefinition = addMethod(gsKnotVector_double_Definition, ...
%    "gismo::gsKnotVector<double> * gismo::gsKnotVector<double>::clone() const", ...
%    "MATLABName", "clone", ...
%    "Description", "clone Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
%defineOutput(cloneDefinition, "RetVal", "clib.gismo.gismo.gsKnotVector_double_", <SHAPE>);
%validate(cloneDefinition);

%% C++ class method |insert| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::insert(double knot,gismo::gsKnotVector<double>::mult_t mult)

insertDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::insert(double knot,gismo::gsKnotVector<double>::mult_t mult)", ...
    "MATLABName", "insert", ...
    "Description", "insert Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(insertDefinition, "knot", "double");
defineArgument(insertDefinition, "mult", "int32");
validate(insertDefinition);

%% C++ class method |remove| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::remove(double const knot,gismo::gsKnotVector<double>::mult_t mult)

removeDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::remove(double const knot,gismo::gsKnotVector<double>::mult_t mult)", ...
    "MATLABName", "remove", ...
    "Description", "remove Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(removeDefinition, "knot", "double");
defineArgument(removeDefinition, "mult", "int32");
validate(removeDefinition);

%% C++ class method |multiplicity| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::multiplicity(double u) const

multiplicityDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::multiplicity(double u) const", ...
    "MATLABName", "multiplicity", ...
    "Description", "multiplicity Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(multiplicityDefinition, "u", "double");
defineOutput(multiplicityDefinition, "RetVal", "int32");
validate(multiplicityDefinition);

%% C++ class method |multFirst| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::multFirst() const

multFirstDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::multFirst() const", ...
    "MATLABName", "multFirst", ...
    "Description", "multFirst Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(multFirstDefinition, "RetVal", "int32");
validate(multFirstDefinition);

%% C++ class method |multLast| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::multLast() const

multLastDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::multLast() const", ...
    "MATLABName", "multLast", ...
    "Description", "multLast Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(multLastDefinition, "RetVal", "int32");
validate(multLastDefinition);

%% C++ class method |maxInteriorMultiplicity| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::maxInteriorMultiplicity() const

maxInteriorMultiplicityDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::maxInteriorMultiplicity() const", ...
    "MATLABName", "maxInteriorMultiplicity", ...
    "Description", "maxInteriorMultiplicity Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(maxInteriorMultiplicityDefinition, "RetVal", "int32");
validate(maxInteriorMultiplicityDefinition);

%% C++ class method |minInteriorMultiplicity| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::minInteriorMultiplicity() const

minInteriorMultiplicityDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::minInteriorMultiplicity() const", ...
    "MATLABName", "minInteriorMultiplicity", ...
    "Description", "minInteriorMultiplicity Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(minInteriorMultiplicityDefinition, "RetVal", "int32");
validate(minInteriorMultiplicityDefinition);

%% C++ class method |multiplicityIndex| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::multiplicityIndex(gismo::gsKnotVector<double>::mult_t i) const

multiplicityIndexDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::mult_t gismo::gsKnotVector<double>::multiplicityIndex(gismo::gsKnotVector<double>::mult_t i) const", ...
    "MATLABName", "multiplicityIndex", ...
    "Description", "multiplicityIndex Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(multiplicityIndexDefinition, "i", "int32");
defineOutput(multiplicityIndexDefinition, "RetVal", "int32");
validate(multiplicityIndexDefinition);

%% C++ class method |affineTransformTo| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::affineTransformTo(double newBeg,double newEnd)

affineTransformToDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::affineTransformTo(double newBeg,double newEnd)", ...
    "MATLABName", "affineTransformTo", ...
    "Description", "affineTransformTo Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(affineTransformToDefinition, "newBeg", "double");
defineArgument(affineTransformToDefinition, "newEnd", "double");
validate(affineTransformToDefinition);

%% C++ class method |size| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: size_t gismo::gsKnotVector<double>::size() const

xSizeDefinition = addMethod(gsKnotVector_double_Definition, ...
    "size_t gismo::gsKnotVector<double>::size() const", ...
    "MATLABName", "xSize", ...
    "Description", "xSize Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(xSizeDefinition, "RetVal", "uint64");
validate(xSizeDefinition);

%% C++ class method |uSize| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: size_t gismo::gsKnotVector<double>::uSize() const

uSizeDefinition = addMethod(gsKnotVector_double_Definition, ...
    "size_t gismo::gsKnotVector<double>::uSize() const", ...
    "MATLABName", "uSize", ...
    "Description", "uSize Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(uSizeDefinition, "RetVal", "uint64");
validate(uSizeDefinition);

%% C++ class method |uValue| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: double gismo::gsKnotVector<double>::uValue(size_t const & i) const

uValueDefinition = addMethod(gsKnotVector_double_Definition, ...
    "double gismo::gsKnotVector<double>::uValue(size_t const & i) const", ...
    "MATLABName", "uValue", ...
    "Description", "uValue Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(uValueDefinition, "i", "uint64", "input");
defineOutput(uValueDefinition, "RetVal", "double");
validate(uValueDefinition);

%% C++ class method |numElements| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: size_t gismo::gsKnotVector<double>::numElements() const

numElementsDefinition = addMethod(gsKnotVector_double_Definition, ...
    "size_t gismo::gsKnotVector<double>::numElements() const", ...
    "MATLABName", "numElements", ...
    "Description", "numElements Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(numElementsDefinition, "RetVal", "uint64");
validate(numElementsDefinition);

%% C++ class method |unique| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::knotContainer gismo::gsKnotVector<double>::unique() const

uniqueDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::knotContainer gismo::gsKnotVector<double>::unique() const", ...
    "MATLABName", "unique", ...
    "Description", "unique Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(uniqueDefinition, "RetVal", "clib.array.gismo.Double");
validate(uniqueDefinition);

%% C++ class method |multiplicities| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::multContainer gismo::gsKnotVector<double>::multiplicities() const

multiplicitiesDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::multContainer gismo::gsKnotVector<double>::multiplicities() const", ...
    "MATLABName", "multiplicities", ...
    "Description", "multiplicities Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(multiplicitiesDefinition, "RetVal", "clib.array.gismo.Int");
validate(multiplicitiesDefinition);

%% C++ class method |data| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: double const * gismo::gsKnotVector<double>::data() const

%dataDefinition = addMethod(gsKnotVector_double_Definition, ...
%    "double const * gismo::gsKnotVector<double>::data() const", ...
%    "MATLABName", "data", ...
%    "Description", "data Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
%defineOutput(dataDefinition, "RetVal", "clib.array.gismo.Double", <SHAPE>); % <MLTYPE> can be "clib.array.gismo.Double", or "double"
%validate(dataDefinition);

%% C++ class method |multSumData| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::mult_t const * gismo::gsKnotVector<double>::multSumData() const

%multSumDataDefinition = addMethod(gsKnotVector_double_Definition, ...
%    "gismo::gsKnotVector<double>::mult_t const * gismo::gsKnotVector<double>::multSumData() const", ...
%    "MATLABName", "multSumData", ...
%    "Description", "multSumData Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
%defineOutput(multSumDataDefinition, "RetVal", "clib.array.gismo.Int", <SHAPE>); % <MLTYPE> can be "clib.array.gismo.Int", or "int32"
%validate(multSumDataDefinition);

%% C++ class method |check| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: bool gismo::gsKnotVector<double>::check() const

checkDefinition = addMethod(gsKnotVector_double_Definition, ...
    "bool gismo::gsKnotVector<double>::check() const", ...
    "MATLABName", "check", ...
    "Description", "check Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(checkDefinition, "RetVal", "logical");
validate(checkDefinition);

%% C++ class method |inDomain| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: bool gismo::gsKnotVector<double>::inDomain(double const u) const

inDomainDefinition = addMethod(gsKnotVector_double_Definition, ...
    "bool gismo::gsKnotVector<double>::inDomain(double const u) const", ...
    "MATLABName", "inDomain", ...
    "Description", "inDomain Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(inDomainDefinition, "u", "double");
defineOutput(inDomainDefinition, "RetVal", "logical");
validate(inDomainDefinition);

%% C++ class method |erase| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::erase(gismo::gsKnotVector<double>::mult_t const first,gismo::gsKnotVector<double>::mult_t const last)

eraseDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::erase(gismo::gsKnotVector<double>::mult_t const first,gismo::gsKnotVector<double>::mult_t const last)", ...
    "MATLABName", "erase", ...
    "Description", "erase Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(eraseDefinition, "first", "int32");
defineArgument(eraseDefinition, "last", "int32");
validate(eraseDefinition);

%% C++ class method |trimLeft| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::trimLeft(gismo::gsKnotVector<double>::mult_t const numKnots)

trimLeftDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::trimLeft(gismo::gsKnotVector<double>::mult_t const numKnots)", ...
    "MATLABName", "trimLeft", ...
    "Description", "trimLeft Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(trimLeftDefinition, "numKnots", "int32");
validate(trimLeftDefinition);

%% C++ class method |trimRight| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::trimRight(gismo::gsKnotVector<double>::mult_t const numKnots)

trimRightDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::trimRight(gismo::gsKnotVector<double>::mult_t const numKnots)", ...
    "MATLABName", "trimRight", ...
    "Description", "trimRight Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(trimRightDefinition, "numKnots", "int32");
validate(trimRightDefinition);

%% C++ class method |numLeftGhosts| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: int gismo::gsKnotVector<double>::numLeftGhosts() const

numLeftGhostsDefinition = addMethod(gsKnotVector_double_Definition, ...
    "int gismo::gsKnotVector<double>::numLeftGhosts() const", ...
    "MATLABName", "numLeftGhosts", ...
    "Description", "numLeftGhosts Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(numLeftGhostsDefinition, "RetVal", "int32");
validate(numLeftGhostsDefinition);

%% C++ class method |numRightGhosts| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: int gismo::gsKnotVector<double>::numRightGhosts() const

numRightGhostsDefinition = addMethod(gsKnotVector_double_Definition, ...
    "int gismo::gsKnotVector<double>::numRightGhosts() const", ...
    "MATLABName", "numRightGhosts", ...
    "Description", "numRightGhosts Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(numRightGhostsDefinition, "RetVal", "int32");
validate(numRightGhostsDefinition);

%% C++ class method |isConsistent| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: static bool gismo::gsKnotVector<double>::isConsistent(gismo::gsKnotVector<double>::knotContainer const & repKnots,gismo::gsKnotVector<double>::multContainer const & multSums)

isConsistentDefinition = addMethod(gsKnotVector_double_Definition, ...
    "static bool gismo::gsKnotVector<double>::isConsistent(gismo::gsKnotVector<double>::knotContainer const & repKnots,gismo::gsKnotVector<double>::multContainer const & multSums)", ...
    "MATLABName", "isConsistent", ...
    "Description", "isConsistent Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(isConsistentDefinition, "repKnots", "clib.array.gismo.Double");
defineArgument(isConsistentDefinition, "multSums", "clib.array.gismo.Int");
defineOutput(isConsistentDefinition, "RetVal", "logical");
validate(isConsistentDefinition);

%% C++ class constructor for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::gsKnotVector(int degree)

gsKnotVector_double_Constructor3Definition = addConstructor(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::gsKnotVector(int degree)", ...
    "Description", "clib.gismo.gismo.gsKnotVector_double_ Constructor of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(gsKnotVector_double_Constructor3Definition, "degree", "int32");
validate(gsKnotVector_double_Constructor3Definition);

%% C++ class constructor for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::gsKnotVector(double first,double last,unsigned int interior,gismo::gsKnotVector<double>::mult_t mult_ends,gismo::gsKnotVector<double>::mult_t mult_interior,int degree)

gsKnotVector_double_Constructor4Definition = addConstructor(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::gsKnotVector(double first,double last,unsigned int interior,gismo::gsKnotVector<double>::mult_t mult_ends,gismo::gsKnotVector<double>::mult_t mult_interior,int degree)", ...
    "Description", "clib.gismo.gismo.gsKnotVector_double_ Constructor of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(gsKnotVector_double_Constructor4Definition, "first", "double");
defineArgument(gsKnotVector_double_Constructor4Definition, "last", "double");
defineArgument(gsKnotVector_double_Constructor4Definition, "interior", "uint32");
defineArgument(gsKnotVector_double_Constructor4Definition, "mult_ends", "int32");
defineArgument(gsKnotVector_double_Constructor4Definition, "mult_interior", "int32");
defineArgument(gsKnotVector_double_Constructor4Definition, "degree", "int32");
validate(gsKnotVector_double_Constructor4Definition);

%% C++ class method |initUniform| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::initUniform(double first,double last,unsigned int interior,unsigned int mult_ends,unsigned int mult_interior,int degree)

initUniformDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::initUniform(double first,double last,unsigned int interior,unsigned int mult_ends,unsigned int mult_interior,int degree)", ...
    "MATLABName", "initUniform", ...
    "Description", "initUniform Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(initUniformDefinition, "first", "double");
defineArgument(initUniformDefinition, "last", "double");
defineArgument(initUniformDefinition, "interior", "uint32");
defineArgument(initUniformDefinition, "mult_ends", "uint32");
defineArgument(initUniformDefinition, "mult_interior", "uint32");
defineArgument(initUniformDefinition, "degree", "int32");
validate(initUniformDefinition);

%% C++ class method |initUniform| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::initUniform(unsigned int numKnots,unsigned int mult_ends,unsigned int mult_interior,int degree)

initUniformDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::initUniform(unsigned int numKnots,unsigned int mult_ends,unsigned int mult_interior,int degree)", ...
    "MATLABName", "initUniform", ...
    "Description", "initUniform Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(initUniformDefinition, "numKnots", "uint32");
defineArgument(initUniformDefinition, "mult_ends", "uint32");
defineArgument(initUniformDefinition, "mult_interior", "uint32");
defineArgument(initUniformDefinition, "degree", "int32");
validate(initUniformDefinition);

%% C++ class method |initGraded| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::initGraded(unsigned int numKnots,int degree,double grading,unsigned int mult_interior)

initGradedDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::initGraded(unsigned int numKnots,int degree,double grading,unsigned int mult_interior)", ...
    "MATLABName", "initGraded", ...
    "Description", "initGraded Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(initGradedDefinition, "numKnots", "uint32");
defineArgument(initGradedDefinition, "degree", "int32");
defineArgument(initGradedDefinition, "grading", "double");
defineArgument(initGradedDefinition, "mult_interior", "uint32");
validate(initGradedDefinition);

%% C++ class method |initGraded| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::initGraded(double u0,double u1,unsigned int interior,int degree,double grading,unsigned int mult_interior)

initGradedDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::initGraded(double u0,double u1,unsigned int interior,int degree,double grading,unsigned int mult_interior)", ...
    "MATLABName", "initGraded", ...
    "Description", "initGraded Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(initGradedDefinition, "u0", "double");
defineArgument(initGradedDefinition, "u1", "double");
defineArgument(initGradedDefinition, "interior", "uint32");
defineArgument(initGradedDefinition, "degree", "int32");
defineArgument(initGradedDefinition, "grading", "double");
defineArgument(initGradedDefinition, "mult_interior", "uint32");
validate(initGradedDefinition);

%% C++ class method |includes| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: bool gismo::gsKnotVector<double>::includes(gismo::gsKnotVector<double> const & other) const

includesDefinition = addMethod(gsKnotVector_double_Definition, ...
    "bool gismo::gsKnotVector<double>::includes(gismo::gsKnotVector<double> const & other) const", ...
    "MATLABName", "includes", ...
    "Description", "includes Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(includesDefinition, "other", "clib.gismo.gismo.gsKnotVector_double_", "input");
defineOutput(includesDefinition, "RetVal", "logical");
validate(includesDefinition);

%% C++ class method |difference| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::difference(gismo::gsKnotVector<double> const & other,std::vector<double, std::allocator<double>> & result) const

differenceDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::difference(gismo::gsKnotVector<double> const & other,std::vector<double, std::allocator<double>> & result) const", ...
    "MATLABName", "difference", ...
    "Description", "difference Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(differenceDefinition, "other", "clib.gismo.gismo.gsKnotVector_double_", "input");
defineArgument(differenceDefinition, "result", "clib.array.gismo.Double");
validate(differenceDefinition);

%% C++ class method |symDifference| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::symDifference(gismo::gsKnotVector<double> const & other,std::vector<double, std::allocator<double>> & result) const

symDifferenceDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::symDifference(gismo::gsKnotVector<double> const & other,std::vector<double, std::allocator<double>> & result) const", ...
    "MATLABName", "symDifference", ...
    "Description", "symDifference Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(symDifferenceDefinition, "other", "clib.gismo.gismo.gsKnotVector_double_", "input");
defineArgument(symDifferenceDefinition, "result", "clib.array.gismo.Double");
validate(symDifferenceDefinition);

%% C++ class method |knotUnion| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double> gismo::gsKnotVector<double>::knotUnion(gismo::gsKnotVector<double> const & b) const

knotUnionDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double> gismo::gsKnotVector<double>::knotUnion(gismo::gsKnotVector<double> const & b) const", ...
    "MATLABName", "knotUnion", ...
    "Description", "knotUnion Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(knotUnionDefinition, "b", "clib.gismo.gismo.gsKnotVector_double_", "input");
defineOutput(knotUnionDefinition, "RetVal", "clib.gismo.gismo.gsKnotVector_double_");
validate(knotUnionDefinition);

%% C++ class method |knotIntersection| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double> gismo::gsKnotVector<double>::knotIntersection(gismo::gsKnotVector<double> const & b) const

knotIntersectionDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double> gismo::gsKnotVector<double>::knotIntersection(gismo::gsKnotVector<double> const & b) const", ...
    "MATLABName", "knotIntersection", ...
    "Description", "knotIntersection Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(knotIntersectionDefinition, "b", "clib.gismo.gismo.gsKnotVector_double_", "input");
defineOutput(knotIntersectionDefinition, "RetVal", "clib.gismo.gismo.gsKnotVector_double_");
validate(knotIntersectionDefinition);

%% C++ class method |reverse| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::reverse()

reverseDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::reverse()", ...
    "MATLABName", "reverse", ...
    "Description", "reverse Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
validate(reverseDefinition);

%% C++ class method |insert| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::insert(gismo::gsKnotVector<double>::knotContainer const & knots,int mult)

insertDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::insert(gismo::gsKnotVector<double>::knotContainer const & knots,int mult)", ...
    "MATLABName", "insert", ...
    "Description", "insert Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(insertDefinition, "knots", "clib.array.gismo.Double");
defineArgument(insertDefinition, "mult", "int32");
validate(insertDefinition);

%% C++ class method |set_degree| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::set_degree(int p)

set_degreeDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::set_degree(int p)", ...
    "MATLABName", "set_degree", ...
    "Description", "set_degree Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(set_degreeDefinition, "p", "int32");
validate(set_degreeDefinition);

%% C++ class method |uniformRefine| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::uniformRefine(gismo::gsKnotVector<double>::mult_t numKnots,gismo::gsKnotVector<double>::mult_t mult)

uniformRefineDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::uniformRefine(gismo::gsKnotVector<double>::mult_t numKnots,gismo::gsKnotVector<double>::mult_t mult)", ...
    "MATLABName", "uniformRefine", ...
    "Description", "uniformRefine Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(uniformRefineDefinition, "numKnots", "int32");
defineArgument(uniformRefineDefinition, "mult", "int32");
validate(uniformRefineDefinition);

%% C++ class method |addConstant| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::addConstant(double amount)

addConstantDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::addConstant(double amount)", ...
    "MATLABName", "addConstant", ...
    "Description", "addConstant Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(addConstantDefinition, "amount", "double");
validate(addConstantDefinition);

%% C++ class method |addConstant| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::addConstant(double start,double amount)

addConstantDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::addConstant(double start,double amount)", ...
    "MATLABName", "addConstant", ...
    "Description", "addConstant Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(addConstantDefinition, "start", "double");
defineArgument(addConstantDefinition, "amount", "double");
validate(addConstantDefinition);

%% C++ class constructor for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::gsKnotVector(gismo::gsKnotVector<double>::knotContainer const & uKnots,int degree,int regularity)

gsKnotVector_double_Constructor5Definition = addConstructor(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::gsKnotVector(gismo::gsKnotVector<double>::knotContainer const & uKnots,int degree,int regularity)", ...
    "Description", "clib.gismo.gismo.gsKnotVector_double_ Constructor of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(gsKnotVector_double_Constructor5Definition, "uKnots", "clib.array.gismo.Double");
defineArgument(gsKnotVector_double_Constructor5Definition, "degree", "int32");
defineArgument(gsKnotVector_double_Constructor5Definition, "regularity", "int32");
validate(gsKnotVector_double_Constructor5Definition);

%% C++ class method |initClamped| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::initClamped(int degree,unsigned int numKnots,unsigned int mult_interior)

initClampedDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::initClamped(int degree,unsigned int numKnots,unsigned int mult_interior)", ...
    "MATLABName", "initClamped", ...
    "Description", "initClamped Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(initClampedDefinition, "degree", "int32");
defineArgument(initClampedDefinition, "numKnots", "uint32");
defineArgument(initClampedDefinition, "mult_interior", "uint32");
validate(initClampedDefinition);

%% C++ class method |initClamped| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::initClamped(double u0,double u1,int degree,unsigned int interior,unsigned int mult_interior)

initClampedDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::initClamped(double u0,double u1,int degree,unsigned int interior,unsigned int mult_interior)", ...
    "MATLABName", "initClamped", ...
    "Description", "initClamped Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(initClampedDefinition, "u0", "double");
defineArgument(initClampedDefinition, "u1", "double");
defineArgument(initClampedDefinition, "degree", "int32");
defineArgument(initClampedDefinition, "interior", "uint32");
defineArgument(initClampedDefinition, "mult_interior", "uint32");
validate(initClampedDefinition);

%% C++ class method |degree| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: int gismo::gsKnotVector<double>::degree() const

degreeDefinition = addMethod(gsKnotVector_double_Definition, ...
    "int gismo::gsKnotVector<double>::degree() const", ...
    "MATLABName", "degree", ...
    "Description", "degree Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(degreeDefinition, "RetVal", "int32");
validate(degreeDefinition);

%% C++ class method |deduceDegree| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: int gismo::gsKnotVector<double>::deduceDegree() const

deduceDegreeDefinition = addMethod(gsKnotVector_double_Definition, ...
    "int gismo::gsKnotVector<double>::deduceDegree() const", ...
    "MATLABName", "deduceDegree", ...
    "Description", "deduceDegree Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(deduceDegreeDefinition, "RetVal", "int32");
validate(deduceDegreeDefinition);

%% C++ class method |greville| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: double gismo::gsKnotVector<double>::greville(int i) const

grevilleDefinition = addMethod(gsKnotVector_double_Definition, ...
    "double gismo::gsKnotVector<double>::greville(int i) const", ...
    "MATLABName", "greville", ...
    "Description", "greville Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(grevilleDefinition, "i", "int32");
defineOutput(grevilleDefinition, "RetVal", "double");
validate(grevilleDefinition);

%% C++ class method |first| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: double gismo::gsKnotVector<double>::first() const

firstDefinition = addMethod(gsKnotVector_double_Definition, ...
    "double gismo::gsKnotVector<double>::first() const", ...
    "MATLABName", "first", ...
    "Description", "first Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(firstDefinition, "RetVal", "double");
validate(firstDefinition);

%% C++ class method |last| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: double gismo::gsKnotVector<double>::last() const

lastDefinition = addMethod(gsKnotVector_double_Definition, ...
    "double gismo::gsKnotVector<double>::last() const", ...
    "MATLABName", "last", ...
    "Description", "last Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(lastDefinition, "RetVal", "double");
validate(lastDefinition);

%% C++ class method |transform| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::transform(double c,double d)

transformDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::transform(double c,double d)", ...
    "MATLABName", "transform", ...
    "Description", "transform Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(transformDefinition, "c", "double");
defineArgument(transformDefinition, "d", "double");
validate(transformDefinition);

%% C++ class method |refineSpans| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::refineSpans(std::vector<unsigned int, std::allocator<unsigned int>> const & spanIndices,gismo::gsKnotVector<double>::mult_t knotsPerSpan)

refineSpansDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::refineSpans(std::vector<unsigned int, std::allocator<unsigned int>> const & spanIndices,gismo::gsKnotVector<double>::mult_t knotsPerSpan)", ...
    "MATLABName", "refineSpans", ...
    "Description", "refineSpans Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(refineSpansDefinition, "spanIndices", "clib.array.gismo.UnsignedInt");
defineArgument(refineSpansDefinition, "knotsPerSpan", "int32");
validate(refineSpansDefinition);

%% C++ class method |refineSpans| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::refineSpans(gismo::gsKnotVector<double>::multContainer const & spanIndices,gismo::gsKnotVector<double>::mult_t knotsPerSpan)

refineSpansDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::refineSpans(gismo::gsKnotVector<double>::multContainer const & spanIndices,gismo::gsKnotVector<double>::mult_t knotsPerSpan)", ...
    "MATLABName", "refineSpans", ...
    "Description", "refineSpans Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(refineSpansDefinition, "spanIndices", "clib.array.gismo.Int");
defineArgument(refineSpansDefinition, "knotsPerSpan", "int32");
validate(refineSpansDefinition);

%% C++ class method |degreeIncrease| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::degreeIncrease(int const & i)

degreeIncreaseDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::degreeIncrease(int const & i)", ...
    "MATLABName", "degreeIncrease", ...
    "Description", "degreeIncrease Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(degreeIncreaseDefinition, "i", "int32", "input");
validate(degreeIncreaseDefinition);

%% C++ class method |degreeDecrease| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::degreeDecrease(int const & i)

degreeDecreaseDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::degreeDecrease(int const & i)", ...
    "MATLABName", "degreeDecrease", ...
    "Description", "degreeDecrease Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(degreeDecreaseDefinition, "i", "int32", "input");
validate(degreeDecreaseDefinition);

%% C++ class method |increaseMultiplicity| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::increaseMultiplicity(gismo::gsKnotVector<double>::mult_t const i,bool boundary)

increaseMultiplicityDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::increaseMultiplicity(gismo::gsKnotVector<double>::mult_t const i,bool boundary)", ...
    "MATLABName", "increaseMultiplicity", ...
    "Description", "increaseMultiplicity Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(increaseMultiplicityDefinition, "i", "int32");
defineArgument(increaseMultiplicityDefinition, "boundary", "logical");
validate(increaseMultiplicityDefinition);

%% C++ class method |reduceMultiplicity| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::reduceMultiplicity(gismo::gsKnotVector<double>::mult_t const i,bool boundary)

reduceMultiplicityDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::reduceMultiplicity(gismo::gsKnotVector<double>::mult_t const i,bool boundary)", ...
    "MATLABName", "reduceMultiplicity", ...
    "Description", "reduceMultiplicity Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(reduceMultiplicityDefinition, "i", "int32");
defineArgument(reduceMultiplicityDefinition, "boundary", "logical");
validate(reduceMultiplicityDefinition);

%% C++ class method |get| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::knotContainer const & gismo::gsKnotVector<double>::get() const

getDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::knotContainer const & gismo::gsKnotVector<double>::get() const", ...
    "MATLABName", "get", ...
    "Description", "get Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(getDefinition, "RetVal", "clib.array.gismo.Double");
validate(getDefinition);

%% C++ class method |has| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: bool gismo::gsKnotVector<double>::has(double knot) const

hasDefinition = addMethod(gsKnotVector_double_Definition, ...
    "bool gismo::gsKnotVector<double>::has(double knot) const", ...
    "MATLABName", "has", ...
    "Description", "has Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(hasDefinition, "knot", "double");
defineOutput(hasDefinition, "RetVal", "logical");
validate(hasDefinition);

%% C++ class method |u_multiplicityIndex| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: unsigned int gismo::gsKnotVector<double>::u_multiplicityIndex(size_t const & i) const

u_multiplicityIndexDefinition = addMethod(gsKnotVector_double_Definition, ...
    "unsigned int gismo::gsKnotVector<double>::u_multiplicityIndex(size_t const & i) const", ...
    "MATLABName", "u_multiplicityIndex", ...
    "Description", "u_multiplicityIndex Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(u_multiplicityIndexDefinition, "i", "uint64", "input");
defineOutput(u_multiplicityIndexDefinition, "RetVal", "uint32");
validate(u_multiplicityIndexDefinition);

%% C++ class method |firstKnotIndex| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: unsigned int gismo::gsKnotVector<double>::firstKnotIndex(size_t const & i) const

firstKnotIndexDefinition = addMethod(gsKnotVector_double_Definition, ...
    "unsigned int gismo::gsKnotVector<double>::firstKnotIndex(size_t const & i) const", ...
    "MATLABName", "firstKnotIndex", ...
    "Description", "firstKnotIndex Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(firstKnotIndexDefinition, "i", "uint64", "input");
defineOutput(firstKnotIndexDefinition, "RetVal", "uint32");
validate(firstKnotIndexDefinition);

%% C++ class method |lastKnotIndex| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: unsigned int gismo::gsKnotVector<double>::lastKnotIndex(size_t const & i) const

lastKnotIndexDefinition = addMethod(gsKnotVector_double_Definition, ...
    "unsigned int gismo::gsKnotVector<double>::lastKnotIndex(size_t const & i) const", ...
    "MATLABName", "lastKnotIndex", ...
    "Description", "lastKnotIndex Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(lastKnotIndexDefinition, "i", "uint64", "input");
defineOutput(lastKnotIndexDefinition, "RetVal", "uint32");
validate(lastKnotIndexDefinition);

%% C++ class method |knotsUntilSpan| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: unsigned int gismo::gsKnotVector<double>::knotsUntilSpan(size_t const & i) const

knotsUntilSpanDefinition = addMethod(gsKnotVector_double_Definition, ...
    "unsigned int gismo::gsKnotVector<double>::knotsUntilSpan(size_t const & i) const", ...
    "MATLABName", "knotsUntilSpan", ...
    "Description", "knotsUntilSpan Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(knotsUntilSpanDefinition, "i", "uint64", "input");
defineOutput(knotsUntilSpanDefinition, "RetVal", "uint32");
validate(knotsUntilSpanDefinition);

%% C++ class method |at| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: double gismo::gsKnotVector<double>::at(size_t const & i) const

atDefinition = addMethod(gsKnotVector_double_Definition, ...
    "double gismo::gsKnotVector<double>::at(size_t const & i) const", ...
    "MATLABName", "at", ...
    "Description", "at Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(atDefinition, "i", "uint64", "input");
defineOutput(atDefinition, "RetVal", "double");
validate(atDefinition);

%% C++ class method |isUniform| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: bool gismo::gsKnotVector<double>::isUniform(double tol) const

isUniformDefinition = addMethod(gsKnotVector_double_Definition, ...
    "bool gismo::gsKnotVector<double>::isUniform(double tol) const", ...
    "MATLABName", "isUniform", ...
    "Description", "isUniform Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(isUniformDefinition, "tol", "double");
defineOutput(isUniformDefinition, "RetVal", "logical");
validate(isUniformDefinition);

%% C++ class method |isOpen| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: bool gismo::gsKnotVector<double>::isOpen() const

isOpenDefinition = addMethod(gsKnotVector_double_Definition, ...
    "bool gismo::gsKnotVector<double>::isOpen() const", ...
    "MATLABName", "isOpen", ...
    "Description", "isOpen Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(isOpenDefinition, "RetVal", "logical");
validate(isOpenDefinition);

%% C++ class method |breaks| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::knotContainer gismo::gsKnotVector<double>::breaks() const

breaksDefinition = addMethod(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::knotContainer gismo::gsKnotVector<double>::breaks() const", ...
    "MATLABName", "breaks", ...
    "Description", "breaks Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(breaksDefinition, "RetVal", "clib.array.gismo.Double");
validate(breaksDefinition);

%% C++ class method |getUniformRefinementKnots| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::getUniformRefinementKnots(gismo::gsKnotVector<double>::mult_t knotsPerSpan,gismo::gsKnotVector<double>::knotContainer & result,gismo::gsKnotVector<double>::mult_t mult) const

getUniformRefinementKnotsDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::getUniformRefinementKnots(gismo::gsKnotVector<double>::mult_t knotsPerSpan,gismo::gsKnotVector<double>::knotContainer & result,gismo::gsKnotVector<double>::mult_t mult) const", ...
    "MATLABName", "getUniformRefinementKnots", ...
    "Description", "getUniformRefinementKnots Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(getUniformRefinementKnotsDefinition, "knotsPerSpan", "int32");
defineArgument(getUniformRefinementKnotsDefinition, "result", "clib.array.gismo.Double");
defineArgument(getUniformRefinementKnotsDefinition, "mult", "int32");
validate(getUniformRefinementKnotsDefinition);

%% C++ class method |detail| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: std::string gismo::gsKnotVector<double>::detail() const

detailDefinition = addMethod(gsKnotVector_double_Definition, ...
    "std::string gismo::gsKnotVector<double>::detail() const", ...
    "MATLABName", "detail", ...
    "Description", "detail Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(detailDefinition, "RetVal", "string");
validate(detailDefinition);

%% C++ class method |maxIntervalLength| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: double gismo::gsKnotVector<double>::maxIntervalLength() const

maxIntervalLengthDefinition = addMethod(gsKnotVector_double_Definition, ...
    "double gismo::gsKnotVector<double>::maxIntervalLength() const", ...
    "MATLABName", "maxIntervalLength", ...
    "Description", "maxIntervalLength Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(maxIntervalLengthDefinition, "RetVal", "double");
validate(maxIntervalLengthDefinition);

%% C++ class method |minIntervalLength| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: double gismo::gsKnotVector<double>::minIntervalLength() const

minIntervalLengthDefinition = addMethod(gsKnotVector_double_Definition, ...
    "double gismo::gsKnotVector<double>::minIntervalLength() const", ...
    "MATLABName", "minIntervalLength", ...
    "Description", "minIntervalLength Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineOutput(minIntervalLengthDefinition, "RetVal", "double");
validate(minIntervalLengthDefinition);

%% C++ class method |degreeElevate| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::degreeElevate(int const & i)

degreeElevateDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::degreeElevate(int const & i)", ...
    "MATLABName", "degreeElevate", ...
    "Description", "degreeElevate Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(degreeElevateDefinition, "i", "int32", "input");
validate(degreeElevateDefinition);

%% C++ class method |degreeReduce| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: void gismo::gsKnotVector<double>::degreeReduce(int const & i)

degreeReduceDefinition = addMethod(gsKnotVector_double_Definition, ...
    "void gismo::gsKnotVector<double>::degreeReduce(int const & i)", ...
    "MATLABName", "degreeReduce", ...
    "Description", "degreeReduce Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(degreeReduceDefinition, "i", "int32", "input");
validate(degreeReduceDefinition);

%% C++ class method |coarsen| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: std::vector<double, std::allocator<double>> gismo::gsKnotVector<double>::coarsen(int knotRemove,int knotSkip,gismo::gsKnotVector<double>::mult_t mul)

coarsenDefinition = addMethod(gsKnotVector_double_Definition, ...
    "std::vector<double, std::allocator<double>> gismo::gsKnotVector<double>::coarsen(int knotRemove,int knotSkip,gismo::gsKnotVector<double>::mult_t mul)", ...
    "MATLABName", "coarsen", ...
    "Description", "coarsen Method of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(coarsenDefinition, "knotRemove", "int32");
defineArgument(coarsenDefinition, "knotSkip", "int32");
defineArgument(coarsenDefinition, "mul", "int32");
defineOutput(coarsenDefinition, "RetVal", "clib.array.gismo.Double");
validate(coarsenDefinition);

%% C++ class constructor for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: gismo::gsKnotVector<double>::gsKnotVector(gismo::gsKnotVector<double> const & input1)

gsKnotVector_double_Constructor6Definition = addConstructor(gsKnotVector_double_Definition, ...
    "gismo::gsKnotVector<double>::gsKnotVector(gismo::gsKnotVector<double> const & input1)", ...
    "Description", "clib.gismo.gismo.gsKnotVector_double_ Constructor of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.
defineArgument(gsKnotVector_double_Constructor6Definition, "input1", "clib.gismo.gismo.gsKnotVector_double_", "input");
validate(gsKnotVector_double_Constructor6Definition);

%% C++ class public data member |m_deg| for C++ class |gismo::gsKnotVector<double>| 
% C++ Signature: int gismo::gsKnotVector<double>::m_deg

addProperty(gsKnotVector_double_Definition, "m_deg", "int32", ...
    "Description", "int32    Data member of C++ class gismo::gsKnotVector<double>."); % Modify help description values as needed.

%% Validate the library definition
validate(libDef);

end
