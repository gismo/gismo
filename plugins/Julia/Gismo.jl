module Gismo

path_to_lib = "build/lib/libgismo"
import Base.print

########################################################################
# CTypes
########################################################################

mutable struct gsCMatrix end
mutable struct gsCFunctionSet end
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

Base.deepcopy(obj::EigenMatrix) = EigenMatrix(rows(obj),cols(obj),data(obj))

Base.show(io::IO, obj::EigenMatrix) = ccall((:gsMatrix_print,path_to_lib),Cvoid,(Ptr{gsCMatrix},),obj.ptr)

function setZero(object::EigenMatrix)::Nothing
    ccall((:gsMatrix_setZero,path_to_lib),Cvoid,(Ptr{gsCMatrix},),object.ptr)
end

########################################################################
# gsGeometry
########################################################################

mutable struct Geometry
    ptr::Ptr{gsCGeometry}

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
    return ccall((:domainDim,path_to_lib),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

function targetDim(object::Geometry)::Cint
    return ccall((:targetDim,path_to_lib),Cint,(Ptr{gsCFunctionSet},),object.ptr)
end

Base.show(io::IO, obj::Geometry) = ccall((:gsFunctionSet_print,path_to_lib),Cvoid,(Ptr{gsCFunctionSet},),obj.ptr)

function eval(obj::Geometry,u::Matrix{Cdouble})::EigenMatrix
    uu = EigenMatrix(size(u,1), size(u,2), pointer(u) )
    result = EigenMatrix()
    ccall((:eval_into,path_to_lib),Cvoid,
	  (Ptr{gsCGeometry},Ptr{gsCMatrix},Ptr{gsCMatrix},),
	  obj.ptr,uu.ptr,result.ptr)
    return result;
end

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

end #module
