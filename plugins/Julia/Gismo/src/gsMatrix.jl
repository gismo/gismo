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
        m = new( ccall((:gsMatrix_create,libgismo),Ptr{gsCMatrix},()) )
        finalizer(destroy, m)
        return m
    end

    function EigenMatrix(r::Int64,c::Int64)
        m = new(ccall((:gsMatrix_create_rc,libgismo),Ptr{gsCMatrix},
                     (Cint,Cint), r, c) )
        finalizer(destroy, m)
        return m
    end

    function EigenMatrix(r::Int,c::Int, data::Ptr{Cdouble})
        m = new(ccall((:gsMatrix_create_rcd,libgismo),Ptr{gsCMatrix},
                     (Cint,Cint,Ptr{Cdouble},), r, c, data) )
        finalizer(destroy, m)
        return m
    end

    function destroy(m::EigenMatrix)
        ccall((:gsMatrix_delete,libgismo),Cvoid,(Ptr{gsCMatrix},),m.ptr)
    end
end

function rows(object::EigenMatrix)::Int64
    return ccall((:gsMatrix_rows,libgismo),Cint,(Ptr{gsCMatrix},),object.ptr)
end

function cols(object::EigenMatrix)::Int64
    return ccall((:gsMatrix_cols,libgismo),Cint,(Ptr{gsCMatrix},),object.ptr)
end

function data(object::EigenMatrix)::Ptr{Cdouble}
    return ccall((:gsMatrix_data,libgismo),Ptr{Cdouble},(Ptr{gsCMatrix},),object.ptr)
end

function asMatrix(object::EigenMatrix)::Matrix{Cdouble}
    return unsafe_wrap(Array, data(object), (rows(object),cols(object)); own = false)
end


Base.deepcopy(obj::EigenMatrix) = EigenMatrix(rows(obj),cols(obj),data(obj))

# Base.show(io::IO, obj::EigenMatrix) = asMatrix(obj)
Base.show(io::IO, obj::EigenMatrix) = ccall((:gsMatrix_print,libgismo),Cvoid,(Ptr{gsCMatrix},),obj.ptr)

function setZero(object::EigenMatrix)::Nothing
    ccall((:gsMatrix_setZero,libgismo),Cvoid,(Ptr{gsCMatrix},),object.ptr)
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
        m = new( ccall((:gsMatrixInt_create,libgismo),Ptr{gsCMatrixInt},()) )
        finalizer(destroy, m)
        return m
    end

    function EigenMatrixInt(r::Int64,c::Int64)
        m = new(ccall((:gsMatrixInt_create_rc,libgismo),Ptr{gsCMatrixInt},
                     (Cint,Cint), r, c) )
        finalizer(destroy, m)
        return m
    end

    function EigenMatrixInt(r::Int,c::Int, data::Ptr{Cint})
        m = new(ccall((:gsMatrixInt_create_rcd,libgismo),Ptr{gsCMatrixInt},
                     (Cint,Cint,Ptr{Cint},), r, c, data) )
        finalizer(destroy, m)
        return m
    end

    function destroy(m::EigenMatrixInt)
        ccall((:gsMatrixInt_delete,libgismo),Cvoid,(Ptr{gsCMatrixInt},),m.ptr)
    end
end

function rows(object::EigenMatrixInt)::Int64
    return ccall((:gsMatrixInt_rows,libgismo),Cint,(Ptr{gsCMatrixInt},),object.ptr)
end

function cols(object::EigenMatrixInt)::Int64
    return ccall((:gsMatrixInt_cols,libgismo),Cint,(Ptr{gsCMatrixInt},),object.ptr)
end

function data(object::EigenMatrixInt)::Ptr{Cint}
    return ccall((:gsMatrixInt_data,libgismo),Ptr{Cint},(Ptr{gsCMatrixInt},),object.ptr)
end

function asMatrixInt(object::EigenMatrixInt)::MatrixInt{Cint}
    return unsafe_wrap(Array, data(object), (rows(object),cols(object)); own = false)
end


Base.deepcopy(obj::EigenMatrixInt) = EigenMatrixInt(rows(obj),cols(obj),data(obj))

# Base.show(io::IO, obj::EigenMatrixInt) = asMatrixInt(obj)
Base.show(io::IO, obj::EigenMatrixInt) = ccall((:gsMatrixInt_print,libgismo),Cvoid,(Ptr{gsCMatrixInt},),obj.ptr)

function setZero(object::EigenMatrixInt)::Nothing
    ccall((:gsMatrixInt_setZero,libgismo),Cvoid,(Ptr{gsCMatrixInt},),object.ptr)
end
