include("./Gismo.jl")
import .Gismo
using BenchmarkTools
using BenchmarkPlots, StatsPlots

path_to_lib::String = "build/lib/libgismo"
geom = Gismo.Geometry( "filedata/surfaces/simple.xml" )


display(@benchmark Gismo.eval(geom,rand(Float64, (2, 1))))

display(@benchmark Gismo.normal(geom,rand(Float64, (2, 1))))

# display(@benchmark Gismo.closest(geom,rand(Float64, (3))))



