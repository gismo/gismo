module Gismo

path_to_lib = "build/lib/libgismo"
import Base.print

########################################################################
# CTypes
########################################################################

mutable struct gsCMatrix end
mutable struct gsCMatrixInt end
mutable struct gsCVector end
mutable struct gsCFunctionSet end
mutable struct gsCMultiPatch end
mutable struct gsCBasis end
mutable struct gsCGeometry end
mutable struct gsCKnotVector end

########################################################################
# gsMatrix
########################################################################

#@ccall library.function_name(argvalue1::argtype1, ...)::returntype

#@inline function wrap_gsCMatrix(m::Ptr{gsCMatrix})
#    return unsafe_wrap(Array{Float64}, gsMatrix_data(m), (gsMatrix_rows(m), gsMatrix_cols(m)))
#end

mutable struct EigenMatrix
    ptr::Ptr{gsCMatrix}

    function EigenMatrix()
        m = new( ccall((:gsMatrix_create,path_to_lib),Ptr{gsCMatrix},()) )
        finalizer(destroy, m)
        return m
    end

    function EigenMatrix(r::Int64,c::Int64)
        m = new(ccall((:gsMatrix_create_rc,path_to_lib),Ptr{gsCMatrix},
                     (Cint,Cint), r, c) )
        finalizer(destroy, m)
        return m
    end

    function EigenMatrix(r::Int,c::Int, data::Ptr{Cdouble})
        m = new(ccall((:gsMatrix_create_rcd,path_to_lib),Ptr{gsCMatrix},
                     (Cint,Cint,Ptr{Cdouble},), r, c, data) )
        finalizer(destroy, m)
        return m
    end

    function destroy(m::EigenMatrix)
        ccall((:gsMatrix_delete,path_to_lib),Cvoid,(Ptr{gsCMatrix},),m.ptr)
    end
end

function rows(object::EigenMatrix)::Int64
    return ccall((:gsMatrix_rows,path_to_lib),Cint,(Ptr{gsCMatrix},),object.ptr)
end

function cols(object::EigenMatrix)::Int64
    return ccall((:gsMatrix_cols,path_to_lib),Cint,(Ptr{gsCMatrix},),object.ptr)
end

function data(object::EigenMatrix)::Ptr{Cdouble}
    return ccall((:gsMatrix_data,path_to_lib),Ptr{Cdouble},(Ptr{gsCMatrix},),object.ptr)
end

function asMatrix(object::EigenMatrix)::Matrix{Cdouble}
    return unsafe_wrap(Array, data(object), (rows(object),cols(object)); own = false)
end


Base.deepcopy(obj::EigenMatrix) = EigenMatrix(rows(obj),cols(obj),data(obj))

# Base.show(io::IO, obj::EigenMatrix) = asMatrix(obj)
Base.show(io::IO, obj::EigenMatrix) = ccall((:gsMatrix_print,path_to_lib),Cvoid,(Ptr{gsCMatrix},),obj.ptr)

function setZero(object::EigenMatrix)::Nothing
    ccall((:gsMatrix_setZero,path_to_lib),Cvoid,(Ptr{gsCMatrix},),object.ptr)
end

########################################################################
# gsMatrixInt
########################################################################

#@ccall library.function_name(argvalue1::argtype1, ...)::returntype

#@inline function wrap_gsCMatrix(m::Ptr{gsCMatrix})
#    return unsafe_wrap(Array{Float64}, gsMatrix_data(m), (gsMatrix_rows(m), gsMatrix_cols(m)))
#end

mutable struct EigenMatrixInt
    ptr::Ptr{gsCMatrixInt}

    function EigenMatrixInt()
        m = new( ccall((:gsMatrixInt_create,path_to_lib),Ptr{gsCMatrixInt},()) )
        finalizer(destroy, m)
        return m
    end

    function EigenMatrixInt(r::Int64,c::Int64)
        m = new(ccall((:gsMatrixInt_create_rc,path_to_lib),Ptr{gsCMatrixInt},
                     (Cint,Cint), r, c) )
        finalizer(destroy, m)
        return m
    end

    function EigenMatrixInt(r::Int,c::Int, data::Ptr{Cint})
        m = new(ccall((:gsMatrixInt_create_rcd,path_to_lib),Ptr{gsCMatrixInt},
                     (Cint,Cint,Ptr{Cint},), r, c, data) )
        finalizer(destroy, m)
        return m
    end

    function destroy(m::EigenMatrixInt)
        ccall((:gsMatrixInt_delete,path_to_lib),Cvoid,(Ptr{gsCMatrixInt},),m.ptr)
    end
end

function rows(object::EigenMatrixInt)::Int64
    return ccall((:gsMatrixInt_rows,path_to_lib),Cint,(Ptr{gsCMatrixInt},),object.ptr)
end

function cols(object::EigenMatrixInt)::Int64
    return ccall((:gsMatrixInt_cols,path_to_lib),Cint,(Ptr{gsCMatrixInt},),object.ptr)
end

function data(object::EigenMatrixInt)::Ptr{Cint}
    return ccall((:gsMatrixInt_data,path_to_lib),Ptr{Cint},(Ptr{gsCMatrixInt},),object.ptr)
end

function asMatrixInt(object::EigenMatrixInt)::MatrixInt{Cint}
    return unsafe_wrap(Array, data(object), (rows(object),cols(object)); own = false)
end


Base.deepcopy(obj::EigenMatrixInt) = EigenMatrixInt(rows(obj),cols(obj),data(obj))

# Base.show(io::IO, obj::EigenMatrixInt) = asMatrixInt(obj)
Base.show(io::IO, obj::EigenMatrixInt) = ccall((:gsMatrixInt_print,path_to_lib),Cvoid,(Ptr{gsCMatrixInt},),obj.ptr)

function setZero(object::EigenMatrixInt)::Nothing
    ccall((:gsMatrixInt_setZero,path_to_lib),Cvoid,(Ptr{gsCMatrixInt},),object.ptr)
end

########################################################################
# gsBasis
########################################################################

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


########################################################################
# gsKnotVector
########################################################################

mutable struct KnotVector
    ptr::Ptr{gsCKnotVector}

    function KnotVector(filename::String)
        g = new(ccall((:gsCReadFile,path_to_lib),Ptr{gsCKnotVector},(Cstring,),filename) )
        finalizer(destroy, g)
        return g
    end

    function KnotVector(a::Vector{Float64})
        kv = new(ccall((:gsKnotVector_create,path_to_lib),Ptr{gsCKnotVector},(Ptr{Cdouble},Cint),a, length(a)) )
        finalizer(destroy, kv)
        return kv
    end

    function destroy(kv::KnotVector)
        ccall((:gsKnotVector_delete,path_to_lib),Cvoid,(Ptr{gsCKnotVector},),kv.ptr)
    end
end

Base.show(io::IO, obj::KnotVector) = ccall((:gsKnotVector_print,path_to_lib),Cvoid,(Ptr{gsCKnotVector},),obj.ptr)


########################################################################
# gsBSplineBasis
########################################################################

function BSplineBasis(kv::KnotVector)::Basis
    b = ccall((:gsBSplineBasis_create,path_to_lib),Ptr{gsCBasis},(Ptr{gsCKnotVector},),kv.ptr)
    return Basis(b)
end

########################################################################
# gsBSpline
########################################################################

function BSpline(basis::Basis,coefs::Matrix{Cdouble})::Geometry
    @assert Base.size(coefs,1)==Gismo.size(basis) "Number of rows of the coefficients should be equal to the number of elements of the basis"
    cc = EigenMatrix(Base.size(coefs,1), Base.size(coefs,2), pointer(coefs) )
    g = ccall((:gsBSpline_create,path_to_lib),Ptr{gsCGeometry},
              (Ptr{gsCBasis},Ptr{gsCMatrix},),
              basis.ptr,cc.ptr)
    return Geometry(g)
end

########################################################################
# gsTensorBSplineBasis
########################################################################

function TensorBSplineBasis(kv1::KnotVector,kv2::KnotVector)::Basis
    b = ccall((:gsTensorBSplineBasis2_create,path_to_lib),Ptr{gsCBasis},(Ptr{gsCKnotVector},Ptr{gsCKnotVector},),kv1.ptr,kv2.ptr)
    return Basis(b)
end

function TensorBSplineBasis(kv1::KnotVector,kv2::KnotVector,kv3::KnotVector)::Basis
    b = ccall((:gsTensorBSplineBasis3_create,path_to_lib),Ptr{gsCBasis},(Ptr{gsCKnotVector},
                                                                         Ptr{gsCKnotVector},
                                                                         Ptr{gsCKnotVector},),
                                                                        kv1.ptr,
                                                                        kv2.ptr,
                                                                        kv3.ptr)
    return Basis(b)
end

function TensorBSplineBasis(kv1::KnotVector,kv2::KnotVector,kv3::KnotVector,kv4::KnotVector)::Basis
    b = ccall((:gsTensorBSplineBasis3_create,path_to_lib),Ptr{gsCBasis},(Ptr{gsCKnotVector},
                                                                         Ptr{gsCKnotVector},
                                                                         Ptr{gsCKnotVector},
                                                                         Ptr{gsCKnotVector},),
                                                                        kv1.ptr,
                                                                        kv2.ptr,
                                                                        kv3.ptr,
                                                                        kv4.ptr)
    return Basis(b)
end


########################################################################
# gsTensorBSpline
########################################################################

function TensorBSpline(basis::Basis,coefs::Matrix{Cdouble})::Geometry
    @assert Base.size(coefs,1)==Gismo.size(basis) "Number of rows of the coefficients should be equal to the number of elements of the basis"
    cc = EigenMatrix(Base.size(coefs,1), Base.size(coefs,2), pointer(coefs) )
    if (domainDim(basis)==2)
        g = ccall((:gsTensorBSpline2_create,path_to_lib),Ptr{gsCGeometry},
              (Ptr{gsCBasis},Ptr{gsCMatrix},),
              basis.ptr,cc.ptr)
    elseif (domainDim(basis)==3)
        g = ccall((:gsTensorBSpline3_create,path_to_lib),Ptr{gsCGeometry},
              (Ptr{gsCBasis},Ptr{gsCMatrix},),
              basis.ptr,cc.ptr)
    elseif (domainDim(basis)==4)
        g = ccall((:gsTensorBSpline4_create,path_to_lib),Ptr{gsCGeometry},
              (Ptr{gsCBasis},Ptr{gsCMatrix},),
              basis.ptr,cc.ptr)
    else
        error("TensorBSpline not implemented for this dimension")
    end
    return Geometry(g)
end

########################################################################
# gsTHBSplineBasis
########################################################################

function THBSplineBasis(basis::Basis)::Basis
    if (domainDim(basis)==1)
        b = ccall((:gsTHBSplineBasis1_create,path_to_lib),Ptr{gsCBasis},(Ptr{gsCBasis},),basis.ptr)
    elseif (domainDim(basis)==2)
        b = ccall((:gsTHBSplineBasis2_create,path_to_lib),Ptr{gsCBasis},(Ptr{gsCBasis},),basis.ptr)
    elseif (domainDim(basis)==3)
        b = ccall((:gsTHBSplineBasis3_create,path_to_lib),Ptr{gsCBasis},(Ptr{gsCBasis},),basis.ptr)
    elseif (domainDim(basis)==4)
        b = ccall((:gsTHBSplineBasis4_create,path_to_lib),Ptr{gsCBasis},(Ptr{gsCBasis},),basis.ptr)
    else
        error("THBSplineBasis not implemented for this dimension")
    end
    return Basis(b)
end

end #module
