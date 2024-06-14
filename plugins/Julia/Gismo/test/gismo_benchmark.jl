include("../src/Gismo.jl")
using BenchmarkTools
using BenchmarkPlots

using .Gismo

geom = Geometry( "surfaces/simple.xml" )

# display(@benchmark eval(geom,rand(Float64, (2, 1))))

display(@benchmark normal(geom,rand(Float64, (2, 1))))



