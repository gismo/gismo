include("./Gismo.jl")
import .Gismo

m = Gismo.EigenMatrix(3,3)
Gismo.setZero(m)
display(m)

geom = Gismo.Geometry( "filedata/surfaces/simple.xml" )
println(geom)
println("Domain dim: ", Gismo.domainDim(geom) )
println("Target dim: ", Gismo.targetDim(geom) )

pts = zeros((2, 3))
pts[:,2] .= 0.5;
pts[:,3] .= 0.99;
println(pts)
result = Gismo.eval(geom,pts)
println("Rows: ", Gismo.rows(result) )
println("Cols: ", Gismo.cols(result) )
display(result)

knots = [0.0,0.0,0.5,1.0,1.0]
kv = Gismo.KnotVector(knots)
println(kv)

#basis = Gismo.TensorBSplineBasis(kv,kv)
#println(basis)

#basis = Gismo.THBSplineBasis(basis)
#println(basis)

#coefs = zeros((4, 3))
#geom = Gismo.THBSpline(basis,coefs)
#println(geom)
