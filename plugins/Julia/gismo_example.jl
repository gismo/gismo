include("./Gismo.jl")
import .Gismo

m = Gismo.EigenMatrix(3,3)
Gismo.setZero(m)
display(Gismo.asMatrix(m))

geom = Gismo.Geometry( "filedata/surfaces/simple.xml" )
println("Geometry:")
println(geom)
println("Domain dim: ", Gismo.domainDim(geom) )
println("Target dim: ", Gismo.targetDim(geom) )

pts = zeros((2, 3))
pts[:,2] .= 0.5;
pts[:,3] .= 0.99;
println("Evaluation points:")
display(pts)

println("Evaluation result:")
values = Gismo.eval(geom,pts)
println("Rows: ", Gismo.rows(values) )
println("Cols: ", Gismo.cols(values) )

vals = Gismo.asMatrix(values)
display(vals)
pts2 = Gismo.invertPoints(geom,vals)
display(Gismo.asMatrix(pts2))

normals = Gismo.normal(geom,pts)
println("Normals at points:")
display(Gismo.asMatrix(normals))

display(Gismo.asMatrix(values)[:,1])
dist,par = Gismo.closest(geom,Gismo.asMatrix(values)[:,1])
display(Gismo.asMatrix(par))


# basis = Gismo.Basis( "filedata/bspbasis/tpBSpline2_01.xml" )
# println("Basis:")
# println(basis)

# values = Gismo.eval(basis,pts)
# println("Rows: ", Gismo.rows(values) )
# println("Cols: ", Gismo.cols(values) )
# vals = Gismo.asMatrix(values)
# display(vals)

kv = Gismo.KnotVector([0.,0.,0.,0.,0.5,1.,1.,1.,1.])
basis = Gismo.TensorBSplineBasis(kv,kv)
println("Basis:")
println(basis)

values = Gismo.eval(basis,pts)
println("Rows: ", Gismo.rows(values) )
println("Cols: ", Gismo.cols(values) )
vals = Gismo.asMatrix(values)
display(vals)

coefs = rand(Gismo.size(basis),3)
geom = Gismo.TensorBSpline(basis,coefs)
values = Gismo.eval(geom,pts)
println("Rows: ", Gismo.rows(values) )
println("Cols: ", Gismo.cols(values) )
vals = Gismo.asMatrix(values)
display(vals)

Gismo.uniformRefine(basis)
println("Basis (refined):")
println(basis)

display(Gismo.actives(basis,pts))

display(Gismo.asMatrix(Gismo.coefs(geom)))

Gismo.THBSplineBasis(basis)
values = Gismo.eval(basis,pts)
println("Rows: ", Gismo.rows(values) )
println("Cols: ", Gismo.cols(values) )
vals = Gismo.asMatrix(values)
display(vals)

mp = Gismo.MultiPatch()
Gismo.addPatch(mp,geom)
display(mp)
display(Gismo.patch(mp,0))
# display(dist)

# knots = [0.0,0.0,0.5,1.0,1.0]
# kv = Gismo.KnotVector(knots)
# println("KnotVector:")
# println(kv)
# println("\n")

#basis = Gismo.TensorBSplineBasis(kv,kv)
#println(basis)

#basis = Gismo.THBSplineBasis(basis)
#println(basis)

#coefs = zeros((4, 3))
#geom = Gismo.THBSpline(basis,coefs)
#println(geom)
