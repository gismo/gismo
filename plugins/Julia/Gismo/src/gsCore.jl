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


"""
Makes a Gismo Basis
"""
mutable struct Basis
    ptr::Ptr{gsCBasis}

    """
    Makes a Gismo Basis from a pointer to a G+Smo basis.
    ...
    # Arguments
    - `basis::Ptr{gsCBasis}`: pointer to a G+Smo basis
    - `delete::Bool=true`: if true, julia will delete the pointer
    ...

    """
    function Basis(basis::Ptr{gsCBasis},delete::Bool=true)
        b = new(basis)
        if (delete)
            finalizer(destroy,b)
        end
        return b
    end

    """
    Makes a Gismo Basis from a file.
    ...
    # Arguments
    - `filename::String`: the name of the file
    ...
    """
    function Basis(filename::String)
        b = new(ccall((:gsCReadFile,libgismo),Ptr{gsCBasis},(Cstring,),filename) )
        finalizer(destroy, b)
        return b
    end

    """
    Deletes a Gismo Basis
    ...
    # Arguments
    - `b::Basis`: a Gismo Basis
    ...
    """
    function destroy(b::Basis)
        ccall((:gsFunctionSet_delete,libgismo),Cvoid,(Ptr{gsCFunctionSet},),b.ptr)
    end
end

"""
Returns the domain dimension of a basis
...
# Arguments
- `object::Basis`: a Gismo Basis
...
"""
function domainDim(object::Basis)::Cint
    return ccall((:gsFunctionSet_domainDim,libgismo),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

"""
Returns the target dimension of a basis
...
# Arguments
- `object::Basis`: a Gismo Basis
...
"""
function targetDim(object::Basis)::Cint
    return ccall((:gsFunctionSet_targetDim,libgismo),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

Base.show(io::IO, obj::Basis) = ccall((:gsFunctionSet_print,libgismo),Cvoid,(Ptr{gsCFunctionSet},),obj.ptr)

"""
Returns the component of a basis
...
# Arguments
- `obj::Basis`: a Gismo Basis
- `i::Cint`: the index of the component
...
"""
function component(obj::Basis,i::Cint)::Basis
    b = ccall((:gsBasis_component,libgismo),Ptr{gsCBasis},(Ptr{gsCBasis},Cint),obj.ptr,i)
    return Basis(b)
end

"""
Returns the degree of a basis
...
# Arguments
- `obj::Basis`: a Gismo Basis
- `i::Cint`: the index of the component
...
"""
function degree(obj::Basis,i::Cint)::Cint
    return ccall((:gsBasis_degree,libgismo),Cint,(Ptr{gsCBasis},Cint),obj.ptr,i)
end

"""
Returns the number of elements of a basis
...
# Arguments
- `obj::Basis`: a Gismo Basis
...
"""
function numElements(obj::Basis)::Cint
    return ccall((:gsBasis_numElements,libgismo),Cint,(Ptr{gsCBasis},),obj.ptr)
end

"""
Returns the size of a basis
...
# Arguments
- `obj::Basis`: a Gismo Basis
...
"""
function size(obj::Basis)::Cint
    return ccall((:gsBasis_size,libgismo),Cint,(Ptr{gsCBasis},),obj.ptr)
end

"""
Refines a basis
...
# Arguments
- `obj::Basis`: a Gismo Basis
- `numKnots::Cint=Int32(1)`: the number of knots to add
- `mul::Cint=Int32(1)`: the multiplicity of the knots
- `dir::Cint=Int32(-1)`: the direction of the refinement
...
"""
function uniformRefine(obj::Basis,numKnots::Cint=Int32(1),mul::Cint=Int32(1),dir::Cint=Int32(-1))::Nothing
    ccall((:gsBasis_uniformRefine,libgismo),Cvoid,
            (Ptr{gsCBasis},Cint,Cint,Cint),obj.ptr,numKnots,mul,dir)
end

"""
Refines a basis
...
# Arguments
- `obj::Basis`: a Gismo Basis
- `boxes::Vector{Cint}`: the boxes to refine
...
"""
function refineElements(obj::Basis,boxes::Vector{Cint})::Nothing
    @assert mod(size(boxes,1),2*domainDim(basis)+1)==0 "Boxes should have size 2*domainDim+1"
    ccall((:gsBasis_refineElements,libgismo),Ptr{gsCBasis},
            (Ptr{gsCBasis},Ptr{Cint},Cint),obj.ptr,refinement,length(refinement))
end

"""
Returns the actives of a basis
...
# Arguments
- `obj::Basis`: a Gismo Basis
- `u::Matrix{Cdouble}`: a matrix of points
...
"""
function actives(obj::Basis,u::Matrix{Cdouble})::EigenMatrixInt
    @assert Base.size(u,1)==domainDim(obj) "Domain dimension should be equal to the number of rows of the points"
    uu = EigenMatrix(Base.size(u,1), Base.size(u,2), pointer(u) )
    result = EigenMatrixInt()
    ccall((:gsBasis_active_into,libgismo),Cvoid,
      (Ptr{gsCBasis},Ptr{gsCMatrix},Ptr{gsCMatrixInt},),
      obj.ptr,uu.ptr,result.ptr)
    return result;
end

"""
Returns the evaluation of a basis
...
# Arguments
- `obj::Basis`: a Gismo Basis
- `u::Matrix{Cdouble}`: a matrix of points
...
"""
function eval(obj::Basis,u::Matrix{Cdouble})::EigenMatrix
    @assert Base.size(u,1)==domainDim(obj) "Domain dimension should be equal to the number of rows of the points"
    uu = EigenMatrix(Base.size(u,1), Base.size(u,2), pointer(u) )
    result = EigenMatrix()
    ccall((:gsFunctionSet_eval_into,libgismo),Cvoid,
      (Ptr{gsCBasis},Ptr{gsCMatrix},Ptr{gsCMatrix},),
      obj.ptr,uu.ptr,result.ptr)
    return result;
end

########################################################################
# gsGeometry
########################################################################

"""
Makes a Gismo Geometry
"""
mutable struct Geometry
    ptr::Ptr{gsCGeometry}

    # delete==true: julia will delete the pointer
    """
    Makes a Gismo Geometry from a pointer to a G+Smo geometry.
    ...
    # Arguments
    - `geom::Ptr{gsCGeometry}`: pointer to a G+Smo geometry
    - `delete::Bool=true`: if true, julia will delete the pointer
    ...

    """
    function Geometry(geom::Ptr{gsCGeometry},delete::Bool=true)
        g = new(geom)

        if (delete)
            finalizer(destroy,g)
        end
        return g
    end

    function Geometry(filename::String)
        g = new(ccall((:gsCReadFile,libgismo),Ptr{gsCGeometry},(Cstring,),filename) )
        finalizer(destroy, g)
        return g
    end

    function destroy(g::Geometry)
        ccall((:gsFunctionSet_delete,libgismo),Cvoid,(Ptr{gsCFunctionSet},),g.ptr)
    end
end

"""
Return the domain dimension of a geometry
...
# Arguments
- `object::Geometry`: a Gismo Geometry
...
"""
function domainDim(object::Geometry)::Cint
    return ccall((:gsFunctionSet_domainDim,libgismo),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

"""
Returns the target dimension of a geometry
...
# Arguments
- `object::Geometry`: a Gismo Geometry
...
"""
function targetDim(object::Geometry)::Cint
    return ccall((:gsFunctionSet_targetDim,libgismo),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

Base.show(io::IO, obj::Geometry) = ccall((:gsFunctionSet_print,libgismo),Cvoid,(Ptr{gsCFunctionSet},),obj.ptr)

"""
Returns the basis of a geometry
...
# Arguments
- `obj::Geometry`: a Gismo Geometry
...
"""
function basis(obj::Geometry)::Basis
    b = ccall((:gsBasis_basis,libgismo),Ptr{gsCBasis},(Ptr{gsCGeometry},),obj.ptr)
    return Basis(b)
end

"""
Returns the coefficients of a geometry
...
# Arguments
- `obj::Geometry`: a Gismo Geometry
...
"""
function coefs(obj::Geometry)::EigenMatrix
    result = EigenMatrix()
    ccall((:gsGeometry_coefs_into,libgismo),Cvoid,(Ptr{gsCGeometry},Ptr{gsCMatrix},),obj.ptr,result.ptr)
    return result;
end

"""
Returns the evaluation of a geometry
...
# Arguments
- `obj::Geometry`: a Gismo Geometry
- `u::Matrix{Cdouble}`: a matrix of points
...
"""
function eval(obj::Geometry,u::Matrix{Cdouble})::EigenMatrix
    @assert Base.size(u,1)==domainDim(obj) "Domain dimension should be equal to the number of rows of the points"
    uu = EigenMatrix(Base.size(u,1), Base.size(u,2), pointer(u) )
    result = EigenMatrix()
    ccall((:gsFunctionSet_eval_into,libgismo),Cvoid,
      (Ptr{gsCGeometry},Ptr{gsCMatrix},Ptr{gsCMatrix},),
      obj.ptr,uu.ptr,result.ptr)
    return result;
end

"""
Returns the normal of a geometry
...
# Arguments
- `obj::Geometry`: a Gismo Geometry
- `u::Matrix{Cdouble}`: a matrix of points
...
"""
function normal(obj::Geometry,u::Matrix{Cdouble})::EigenMatrix
    @assert Base.size(u,1)==domainDim(obj) "Domain dimension should be equal to the number of rows of the points"
    uu = EigenMatrix(Base.size(u,1), Base.size(u,2), pointer(u) )
    result = EigenMatrix()
    ccall((:gsGeometry_normal_into,libgismo),Cvoid,
      (Ptr{gsCGeometry},Ptr{gsCMatrix},Ptr{gsCMatrix},),
      obj.ptr,uu.ptr,result.ptr)
    return result;
end

"""
Returns the closest point of a geometry
...
# Arguments
- `obj::Geometry`: a Gismo Geometry
- `x::Vector{Cdouble}`: a vector of points
- `accuracy::Cdouble=1e-6`: the accuracy of the computation
...
"""
function closest(obj::Geometry,x::Vector{Cdouble},accuracy::Cdouble=1e-6)::Tuple{Cdouble,EigenMatrix}
    xx = EigenMatrix(Base.size(x,1), Base.size(x,2), pointer(x) )
    result = EigenMatrix()
    dist = ccall((:gsGeometry_closestPointTo,libgismo),Cdouble,
      (Ptr{gsCGeometry},Ptr{gsCMatrix},Ptr{gsCMatrix},Cdouble,),
      obj.ptr,xx.ptr,result.ptr,accuracy)
    # display(result)
    # display(result)
    # display(result)
    return dist,result;
end

"""
Inverts a set of points
...
# Arguments
- `obj::Geometry`: a Gismo Geometry
- `x::Matrix{Cdouble}`: a matrix of points
- `accuracy::Cdouble=1e-6`: the accuracy of the computation
...
"""
function invertPoints(obj::Geometry,x::Matrix{Cdouble},accuracy::Cdouble=1e-6)::EigenMatrix
    xx = EigenMatrix(Base.size(x,1), 1, pointer(x) )
    result = EigenMatrix()
    ccall((:gsGeometry_invertPoints,libgismo),Cvoid,
      (Ptr{gsCGeometry},Ptr{gsCMatrix},Ptr{gsCMatrix},Cdouble,),
      obj.ptr,xx.ptr,result.ptr,accuracy)
    return result;
end

########################################################################
# gsMultiPatch
########################################################################

"""
Makes a Gismo MultiPatch
"""
mutable struct MultiPatch
    ptr::Ptr{gsCMultiPatch}

    # function MultiPatch(filename::String)
    #     g = new(ccall((:gsCReadFile,libgismo),Ptr{gsCMultiPatch},(Cstring,),filename) )
    #     finalizer(destroy, g)
    #     return g
    # end

    function MultiPatch()
        m = new(ccall((:gsMultiPatch_create,libgismo),Ptr{gsCMultiPatch},(),) )
        finalizer(destroy, m)
        return m
    end

    function destroy(m::MultiPatch)
        ccall((:gsMultiPatch_delete,libgismo),Cvoid,(Ptr{gsCFunctionSet},),m.ptr)
    end
end

"""
Adds a patch to a MultiPatch
...
# Arguments
- `obj::MultiPatch`: a Gismo MultiPatch
- `geom::Geometry`: a Gismo Geometry
...
"""
function domainDim(object::MultiPatch)::Cint
    return ccall((:gsFunctionSet_domainDim,libgismo),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

"""
Returns the target dimension of a MultiPatch
...
# Arguments
- `object::MultiPatch`: a Gismo MultiPatch
...
"""
function targetDim(object::MultiPatch)::Cint
    return ccall((:gsFunctionSet_targetDim,libgismo),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

Base.show(io::IO, obj::MultiPatch) = ccall((:gsFunctionSet_print,libgismo),Cvoid,(Ptr{gsCFunctionSet},),obj.ptr)

"""
Returns the basis of a MultiPatch
...
# Arguments
- `obj::MultiPatch`: a Gismo MultiPatch
...
"""
function addPatch(obj::MultiPatch,geom::Geometry)::Nothing
    ccall((:gsMultiPatch_addPatch,libgismo),Cvoid,(Ptr{gsCMultiPatch},Ptr{gsCGeometry},),obj.ptr,geom.ptr)
end

"""
Returns the coefficients of a MultiPatch
...
# Arguments
- `obj::MultiPatch`: a Gismo MultiPatch
...
"""
function basis(obj::MultiPatch,i::Int)::Basis
    b = ccall((:gsMultiPatch_basis,libgismo),Ptr{gsCBasis},(Ptr{gsCFunctionSet},Cint),obj.ptr,i)
    return Basis(b,false)
end

"""
Returns the coefficients of a MultiPatch
...
# Arguments
- `obj::MultiPatch`: a Gismo MultiPatch
...
"""
function patch(obj::MultiPatch,i::Int)::Geometry
    g = ccall((:gsMultiPatch_patch,libgismo),Ptr{gsCGeometry},(Ptr{gsCMultiPatch},Cint),obj.ptr,i)
    return Geometry(g,false)
end

########################################################################
# gsMultiBasis
########################################################################
