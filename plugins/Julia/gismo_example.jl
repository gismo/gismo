include("./Gismo.jl")
import .Gismo

path_to_lib::String = "build/lib/libgismo"
geom = Gismo.read("filedata/domain2d/lake.xml")
Gismo.print(geom)

pts = zeros((2, 3))
display(pts)

result = Gismo.eval(geom,pts)

knots = [0.0,0.0,1.0,1.0]
display(result)
# Print the geometry

kv = Gismo.knotVector(knots)
# basis = Gismo.bsplineBasis(kv)
basis = Gismo.tensorBSplineBasis(kv,kv)
basis = Gismo.THBSplineBasis(basis)
Gismo.print(basis)

coefs = zeros((4, 3))
geom = Gismo.THBSpline(basis,coefs)

Gismo.print(geom)

