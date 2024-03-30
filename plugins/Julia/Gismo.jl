module Gismo

path_to_lib = "build/lib/libgismo"

########################################################################
# CTypes
########################################################################

mutable struct gsCMatrix end
mutable struct gsCFunctionSet end
mutable struct gsCGeometry end

########################################################################
# gsMatrix
########################################################################

@inline function wrap_gsCMatrix(m::Ptr{gsCMatrix})
    M = unsafe_load(m)
    @assert M.size2==M.tda "Cannot unsafe_wrap gsl_matrix with tda != size2."     
    return unsafe_wrap(Array{Float64}, gsMatrix_data(m), (gsMatrix_rows(m), gsMatrix_cols(m)))
end

mutable struct EigenMatrix
    ptr::Ptr{gsCMatrix}

    function EigenMatrix()
        m = new( ccall((:gsMatrix_create,path_to_lib),Ptr{gsCMatrix},()) )
        finalizer(destroy, m)
        return m
    end

    function EigenMatrix(r::Int64,c::Int64)
        m = new(ccall((:gsMatrix_create_rc,path_to_lib),Ptr{Cdouble},
                     (Cint,Cint), r, c) )
        finalizer(destroy, m)
        return m
    end

    function destroy(m::EigenMatrix)
        ccall((:gsMatrix_delete,path_to_lib),Cvoid,(Ptr{gsCMatrix},),m.ptr)
    end
end

function print(object::EigenMatrix)::Nothing
    ccall((:gsMatrix_print,path_to_lib),Cvoid,(Ptr{gsCMatrix},),object.ptr)
    #printf("\n")
    return;
end

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

function print(object::Geometry)::Nothing
	ccall((:gsFunctionSet_print,path_to_lib),Cvoid,(Ptr{gsCFunctionSet},),object.ptr)
	return;
end

function eval(obj::Geometry,u::Matrix{Cdouble})::EigenMatrix
    uu = ccall((:gsMatrix_create_rcd,path_to_lib),Ptr{Cdouble},
                     (Cint,Cint,Ptr{Cdouble}), size(u,1), size(u,2), u)
    result = EigenMatrix()
    ccall((:eval_into,path_to_lib),Ptr{gsCMatrix},
	  ( Ptr{gsCFunctionSet},Ptr{gsCMatrix},Ptr{gsCMatrix}),
	  obj.ptr,uu,result.ptr)
    return result;
end

########################################################################
# gsReadFile
########################################################################

end #module

