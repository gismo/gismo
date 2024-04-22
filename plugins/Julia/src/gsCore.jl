########################################################################
# gsBasis
########################################################################

export
    Basis,
    domainDim,
    targetDim,
    component,
    degree,
    numElements,
    size,
    uniformRefine,
    refineElements,
    actives,
    eval,
    Geometry,
    basis,
    coefs,
    normal,
    closest,
    invertPoints,
    MultiPatch,
    addPatch,
    patch



mutable struct Basis
    ptr::Ptr{gsCBasis}

    function Basis(basis::Ptr{gsCBasis})
        b = new(basis)
        finalizer(destroy,b)
        return b
    end

    function Basis(filename::String)
        b = new(ccall((:gsCReadFile,path_to_lib),Ptr{gsCBasis},(Cstring,),filename) )
        finalizer(destroy, b)
        return b
    end

    function destroy(b::Basis)
        ccall((:gsFunctionSet_delete,path_to_lib),Cvoid,(Ptr{gsCFunctionSet},),b.ptr)
    end
end

function domainDim(object::Basis)::Cint
    return ccall((:gsFunctionSet_domainDim,path_to_lib),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

function targetDim(object::Basis)::Cint
    return ccall((:gsFunctionSet_targetDim,path_to_lib),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

Base.show(io::IO, obj::Basis) = ccall((:gsFunctionSet_print,path_to_lib),Cvoid,(Ptr{gsCFunctionSet},),obj.ptr)

function component(obj::Basis,i::Cint)::Basis
    b = ccall((:gsBasis_component,path_to_lib),Ptr{gsCBasis},(Ptr{gsCBasis},Cint),obj.ptr,i)
    return Basis(b)
end

function degree(obj::Basis,i::Cint)::Cint
    return ccall((:gsBasis_degree,path_to_lib),Cint,(Ptr{gsCBasis},Cint),obj.ptr,i)
end

function numElements(obj::Basis)::Cint
    return ccall((:gsBasis_numElements,path_to_lib),Cint,(Ptr{gsCBasis},),obj.ptr)
end

function size(obj::Basis)::Cint
    return ccall((:gsBasis_size,path_to_lib),Cint,(Ptr{gsCBasis},),obj.ptr)
end

function uniformRefine(obj::Basis,numKnots::Cint=Int32(1),mul::Cint=Int32(1),dir::Cint=Int32(-1))::Nothing
    ccall((:gsBasis_uniformRefine,path_to_lib),Cvoid,
            (Ptr{gsCBasis},Cint,Cint,Cint),obj.ptr,numKnots,mul,dir)
end

function refineElements(obj::Basis,boxes::Vector{Cint})::Nothing
    @assert mod(size(boxes,1),2*domainDim(basis)+1)==0 "Boxes should have size 2*domainDim+1"
    ccall((:gsBasis_refineElements,path_to_lib),Ptr{gsCBasis},
            (Ptr{gsCBasis},Ptr{Cint},Cint),obj.ptr,refinement,length(refinement))
end


function actives(obj::Basis,u::Matrix{Cdouble})::EigenMatrixInt
    @assert Base.size(u,1)==domainDim(obj) "Domain dimension should be equal to the number of rows of the points"
    uu = EigenMatrix(Base.size(u,1), Base.size(u,2), pointer(u) )
    result = EigenMatrixInt()
    ccall((:gsBasis_active_into,path_to_lib),Cvoid,
      (Ptr{gsCBasis},Ptr{gsCMatrix},Ptr{gsCMatrixInt},),
      obj.ptr,uu.ptr,result.ptr)
    return result;
end

function eval(obj::Basis,u::Matrix{Cdouble})::EigenMatrix
    @assert Base.size(u,1)==domainDim(obj) "Domain dimension should be equal to the number of rows of the points"
    uu = EigenMatrix(Base.size(u,1), Base.size(u,2), pointer(u) )
    result = EigenMatrix()
    ccall((:gsFunctionSet_eval_into,path_to_lib),Cvoid,
      (Ptr{gsCBasis},Ptr{gsCMatrix},Ptr{gsCMatrix},),
      obj.ptr,uu.ptr,result.ptr)
    return result;
end

########################################################################
# gsGeometry
########################################################################

mutable struct Geometry
    ptr::Ptr{gsCGeometry}

    function Geometry(geom::Ptr{gsCGeometry})
        g = new(geom)
        finalizer(destroy,g)
        return g
    end

    function Geometry(filename::String)
        g = new(ccall((:gsCReadFile,path_to_lib),Ptr{gsCGeometry},(Cstring,),filename) )
        finalizer(destroy, g)
        return g
    end

    function destroy(g::Geometry)
        ccall((:gsFunctionSet_delete,path_to_lib),Cvoid,(Ptr{gsCFunctionSet},),g.ptr)
    end
end

function domainDim(object::Geometry)::Cint
    return ccall((:gsFunctionSet_domainDim,path_to_lib),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

function targetDim(object::Geometry)::Cint
    return ccall((:gsFunctionSet_targetDim,path_to_lib),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

Base.show(io::IO, obj::Geometry) = ccall((:gsFunctionSet_print,path_to_lib),Cvoid,(Ptr{gsCFunctionSet},),obj.ptr)

function basis(obj::Geometry)::Basis
    b = ccall((:gsBasis_basis,path_to_lib),Ptr{gsCBasis},(Ptr{gsCGeometry},),obj.ptr)
    return Basis(b)
end

function coefs(obj::Geometry)::EigenMatrix
    result = EigenMatrix()
    ccall((:gsGeometry_coefs_into,path_to_lib),Cvoid,(Ptr{gsCGeometry},Ptr{gsCMatrix},),obj.ptr,result.ptr)
    return result;
end

function eval(obj::Geometry,u::Matrix{Cdouble})::EigenMatrix
    @assert Base.size(u,1)==domainDim(obj) "Domain dimension should be equal to the number of rows of the points"
    uu = EigenMatrix(Base.size(u,1), Base.size(u,2), pointer(u) )
    result = EigenMatrix()
    ccall((:gsFunctionSet_eval_into,path_to_lib),Cvoid,
      (Ptr{gsCGeometry},Ptr{gsCMatrix},Ptr{gsCMatrix},),
      obj.ptr,uu.ptr,result.ptr)
    return result;
end

function normal(obj::Geometry,u::Matrix{Cdouble})::EigenMatrix
    @assert Base.size(u,1)==domainDim(obj) "Domain dimension should be equal to the number of rows of the points"
    uu = EigenMatrix(Base.size(u,1), Base.size(u,2), pointer(u) )
    result = EigenMatrix()
    ccall((:gsGeometry_normal_into,path_to_lib),Cvoid,
      (Ptr{gsCGeometry},Ptr{gsCMatrix},Ptr{gsCMatrix},),
      obj.ptr,uu.ptr,result.ptr)
    return result;
end

function closest(obj::Geometry,x::Vector{Cdouble},accuracy::Cdouble=1e-6)::Tuple{Cdouble,EigenMatrix}
    xx = EigenMatrix(Base.size(x,1), Base.size(x,2), pointer(x) )
    result = EigenMatrix()
    dist = ccall((:gsGeometry_closestPointTo,path_to_lib),Cdouble,
      (Ptr{gsCGeometry},Ptr{gsCMatrix},Ptr{gsCMatrix},Cdouble,),
      obj.ptr,xx.ptr,result.ptr,accuracy)
    # display(result)
    # display(result)
    # display(result)
    return dist,result;
end

function invertPoints(obj::Geometry,x::Matrix{Cdouble},accuracy::Cdouble=1e-6)::EigenMatrix
    xx = EigenMatrix(Base.size(x,1), 1, pointer(x) )
    result = EigenMatrix()
    ccall((:gsGeometry_invertPoints,path_to_lib),Cvoid,
      (Ptr{gsCGeometry},Ptr{gsCMatrix},Ptr{gsCMatrix},Cdouble,),
      obj.ptr,xx.ptr,result.ptr,accuracy)
    return result;
end

########################################################################
# gsMultiPatch
########################################################################

mutable struct MultiPatch
    ptr::Ptr{gsCMultiPatch}

    # function MultiPatch(filename::String)
    #     g = new(ccall((:gsCReadFile,path_to_lib),Ptr{gsCMultiPatch},(Cstring,),filename) )
    #     finalizer(destroy, g)
    #     return g
    # end

    function MultiPatch()
        m = new(ccall((:gsMultiPatch_create,path_to_lib),Ptr{gsCMultiPatch},(),) )
        finalizer(destroy, m)
        return m
    end

    function destroy(m::MultiPatch)
        ccall((:gsFunctionSet_delete,path_to_lib),Cvoid,(Ptr{gsCFunctionSet},),m.ptr)
    end
end

function domainDim(object::MultiPatch)::Cint
    return ccall((:gsFunctionSet_domainDim,path_to_lib),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

function targetDim(object::MultiPatch)::Cint
    return ccall((:gsFunctionSet_targetDim,path_to_lib),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

Base.show(io::IO, obj::MultiPatch) = ccall((:gsFunctionSet_print,path_to_lib),Cvoid,(Ptr{gsCFunctionSet},),obj.ptr)

function addPatch(obj::MultiPatch,geom::Geometry)::Nothing
    ccall((:gsMultiPatch_addPatch,path_to_lib),Cvoid,(Ptr{gsCMultiPatch},Ptr{gsCGeometry},),obj.ptr,geom.ptr)
end

function basis(obj::MultiPatch,i::Int)::Basis
    b = ccall((:gsMultiPatch_basis,path_to_lib),Ptr{gsCBasis},(Ptr{gsCFunctionSet},Cint),obj.ptr,i)
    return Basis(b)
end

function patch(obj::MultiPatch,i::Int)::Geometry
    g = ccall((:gsMultiPatch_patch,path_to_lib),Ptr{gsCGeometry},(Ptr{gsCMultiPatch},Cint),obj.ptr,i)
    return Geometry(g)
end

########################################################################
# gsMultiBasis
########################################################################
