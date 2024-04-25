########################################################################
# gsKnotVector
########################################################################

mutable struct KnotVector
    ptr::Ptr{gsCKnotVector}

    function KnotVector(filename::String)
        g = new(ccall((:gsCReadFile,libgismo),Ptr{gsCKnotVector},(Cstring,),filename) )
        finalizer(destroy, g)
        return g
    end

    function KnotVector(a::Vector{Float64})
        kv = new(ccall((:gsKnotVector_create,libgismo),Ptr{gsCKnotVector},(Ptr{Cdouble},Cint),a, length(a)) )
        finalizer(destroy, kv)
        return kv
    end

    function destroy(kv::KnotVector)
        ccall((:gsKnotVector_delete,libgismo),Cvoid,(Ptr{gsCKnotVector},),kv.ptr)
    end
end

Base.show(io::IO, obj::KnotVector) = ccall((:gsKnotVector_print,libgismo),Cvoid,(Ptr{gsCKnotVector},),obj.ptr)


########################################################################
# gsBSplineBasis
########################################################################

function BSplineBasis(kv::KnotVector)::Basis
    b = ccall((:gsBSplineBasis_create,libgismo),Ptr{gsCBasis},(Ptr{gsCKnotVector},),kv.ptr)
    return Basis(b)
end

########################################################################
# gsBSpline
########################################################################

function BSpline(basis::Basis,coefs::Matrix{Cdouble})::Geometry
    @assert Base.size(coefs,1)==Gismo.size(basis) "Number of rows of the coefficients should be equal to the number of elements of the basis"
    cc = EigenMatrix(Base.size(coefs,1), Base.size(coefs,2), pointer(coefs) )
    g = ccall((:gsBSpline_create,libgismo),Ptr{gsCGeometry},
              (Ptr{gsCBasis},Ptr{gsCMatrix},),
              basis.ptr,cc.ptr)
    return Geometry(g)
end

########################################################################
# gsTensorBSplineBasis
########################################################################

function TensorBSplineBasis(kv1::KnotVector,kv2::KnotVector)::Basis
    b = ccall((:gsTensorBSplineBasis2_create,libgismo),Ptr{gsCBasis},(Ptr{gsCKnotVector},Ptr{gsCKnotVector},),kv1.ptr,kv2.ptr)
    return Basis(b)
end

function TensorBSplineBasis(kv1::KnotVector,kv2::KnotVector,kv3::KnotVector)::Basis
    b = ccall((:gsTensorBSplineBasis3_create,libgismo),Ptr{gsCBasis},(Ptr{gsCKnotVector},
                                                                         Ptr{gsCKnotVector},
                                                                         Ptr{gsCKnotVector},),
                                                                        kv1.ptr,
                                                                        kv2.ptr,
                                                                        kv3.ptr)
    return Basis(b)
end

function TensorBSplineBasis(kv1::KnotVector,kv2::KnotVector,kv3::KnotVector,kv4::KnotVector)::Basis
    b = ccall((:gsTensorBSplineBasis3_create,libgismo),Ptr{gsCBasis},(Ptr{gsCKnotVector},
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
        g = ccall((:gsTensorBSpline2_create,libgismo),Ptr{gsCGeometry},
              (Ptr{gsCBasis},Ptr{gsCMatrix},),
              basis.ptr,cc.ptr)
    elseif (domainDim(basis)==3)
        g = ccall((:gsTensorBSpline3_create,libgismo),Ptr{gsCGeometry},
              (Ptr{gsCBasis},Ptr{gsCMatrix},),
              basis.ptr,cc.ptr)
    elseif (domainDim(basis)==4)
        g = ccall((:gsTensorBSpline4_create,libgismo),Ptr{gsCGeometry},
              (Ptr{gsCBasis},Ptr{gsCMatrix},),
              basis.ptr,cc.ptr)
    else
        error("TensorBSpline not implemented for this dimension")
    end
    return Geometry(g)
end
