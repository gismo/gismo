# Inspired by https://github.com/precice/PreCICE.jl/blob/develop/src/PreCICE.jl

module Gismo

path_to_lib::String = "build/lib/libgismo"

function print(object::Ptr{Cvoid})::Nothing
	ccall((:print,path_to_lib),Cvoid,(Ptr{Cvoid},),object)
	return;
end

function read(filename::String)::Ptr{Cvoid}
	return ccall((:read,path_to_lib),Ptr{Cvoid},(Cstring,),filename)
end

function domainDim(object::Ptr{Cvoid})::Cint
	return ccall((:domainDim,path_to_lib),Cint,(Ptr{Cvoid},),object)
end

function targetDim(object::Ptr{Cvoid})::Cint
	return ccall((:targetDim,path_to_lib),Cint,(Ptr{Cvoid},),object)
end

function eval(object::Ptr{Cvoid},points::Array{Cdouble})::Array{Cdouble}
	@assert ndims(points)==2
	_rows = size(points,1)
	_cols = size(points,2)
	_size = _cols * targetDim(object);
	result = Array{Cdouble,1}(undef, _size)

	ccall((:eval_into,path_to_lib),Cvoid,
		  ( Ptr{Cvoid},
		  	Ptr{Cdouble},
		  	Cint,
		  	Cint,
		  	Ref{Cdouble},
		  	Cint,
		  	Ref{Cint},
		  	Ref{Cint}),
			object,
			points,
			_rows,
			_cols,
			result,
			_size,
			_rows,
			_cols)
	return reshape(result,_rows,_cols);
end

function knotVector(knots::Array{Cdouble,1})::Ptr{Cvoid}
	_size = size(knots,1)
	return ccall((:knotVector,path_to_lib),Ptr{Cvoid},(Ptr{Cdouble},Cint,),knots,_size)
end

function bsplineBasis(KV::Ptr{Cvoid})::Ptr{Cvoid}
	return ccall((:bsplineBasis,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},),KV)
end

function tensorBSplineBasis(KV1::Ptr{Cvoid},KV2::Ptr{Cvoid})::Ptr{Cvoid}
	return ccall((:tensorBSplineBasis2,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cvoid},),KV1,KV2)
end

function tensorBSplineBasis(KV1::Ptr{Cvoid},KV2::Ptr{Cvoid},KV3::Ptr{Cvoid})::Ptr{Cvoid}
	return ccall((:tensorBSplineBasis3,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cvoid},Ptr{Cvoid},),KV1,KV2,KV3)
end

function tensorBSplineBasis(KV1::Ptr{Cvoid},KV2::Ptr{Cvoid},KV3::Ptr{Cvoid},KV4::Ptr{Cvoid})::Ptr{Cvoid}
	return ccall((:tensorBSplineBasis4,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cvoid},Ptr{Cvoid},Ptr{Cvoid},),KV1,KV2,KV3,KV4)
end

function tensorBSplineBasis(KV1::Ptr{Cvoid},KV2::Ptr{Cvoid},KV3::Ptr{Cvoid})::Ptr{Cvoid}
	return ccall((:tensorBSplineBasis3,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cvoid},Ptr{Cvoid},),KV1,KV2,KV3)
end

function tensorBSplineBasis(KV1::Ptr{Cvoid},KV2::Ptr{Cvoid},KV3::Ptr{Cvoid},KV4::Ptr{Cvoid})::Ptr{Cvoid}
	return ccall((:tensorBSplineBasis4,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cvoid},Ptr{Cvoid},Ptr{Cvoid},),KV1,KV2,KV3,KV4)
end

function tensorBSpline(basis::Ptr{Cvoid},coefs::Array{Cdouble})::Ptr{Cvoid}
	@assert ndims(coefs)==2
	_rows = size(coefs,1)
	_cols = size(coefs,2)
	if     domainDim(basis) == 2
		return ccall((:tensorBSpline2,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cdouble},Cint,Cint,),basis,coefs,_rows,_cols)
	elseif domainDim(basis) == 3
		return ccall((:tensorBSpline3,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cdouble},Cint,Cint,),basis,coefs,_rows,_cols)
	elseif domainDim(basis) == 4
		return ccall((:tensorBSpline4,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cdouble},Cint,Cint,),basis,coefs,_rows,_cols)
	end
end

function THBSplineBasis(basis::Ptr{Cvoid})::Ptr{Cvoid}
	if     domainDim(basis) == 1
		return ccall((:THBSplineBasis1,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},),basis)
	elseif domainDim(basis) == 2
		return ccall((:THBSplineBasis2,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},),basis)
	elseif domainDim(basis) == 3
		return ccall((:THBSplineBasis3,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},),basis)
	elseif domainDim(basis) == 4
		return ccall((:THBSplineBasis4,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},),basis)
	end
end

function THBSpline(basis::Ptr{Cvoid},coefs::Array{Cdouble})::Ptr{Cvoid}
	@assert ndims(coefs)==2
	_rows = size(coefs,1)
	_cols = size(coefs,2)
	if     domainDim(basis) == 1
		return ccall((:THBSpline1,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cdouble},Cint,Cint,),basis,coefs,_rows,_cols)
	elseif domainDim(basis) == 2
		return ccall((:THBSpline2,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cdouble},Cint,Cint,),basis,coefs,_rows,_cols)
	elseif domainDim(basis) == 3
		return ccall((:THBSpline3,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cdouble},Cint,Cint,),basis,coefs,_rows,_cols)
	elseif domainDim(basis) == 4
		return ccall((:THBSpline4,path_to_lib),Ptr{Cvoid},(Ptr{Cvoid},Ptr{Cdouble},Cint,Cint,),basis,coefs,_rows,_cols)
	end
end

# function eval(pointer::Ptr{Cvoid})::Arrayxxxx
# 	ccall((:eval,path_to_lib),Cvoid,(Ptr{Cvoid},),geom,points)
# end

# function deriv(pointer::Ptr{Cvoid})::Arrayxxxx
# 	ccall((:eval,path_to_lib),Cvoid,(Ptr{Cvoid},),geom,points)
# end

# function dim(pointer::Ptr{Cvoid})::Arrayxxxx
# 	ccall((:dim,path_to_lib),Cvoid,(Ptr{Cvoid},),geom)
# end

# ##
# ## @brief      Computes second derivative
# ##
# ## @param      pointer object to compute on
# ## @param      points  points to compute
# ##
# ## @return     { description_of_the_return_value }
# ##
# function deriv2(pointer::Ptr{Cvoid})::Arrayxxxx
# 	ccall((:eval,path_to_lib),Cvoid,(Ptr{Cvoid},),geom,points)
# end

end
