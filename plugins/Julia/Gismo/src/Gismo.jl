module Gismo

# Get the path to the library
include("Gismo_path.jl")

# Forward declaration of structs
include("Declarations.jl")

# Load gsMatrix
include("gsMatrix.jl")

# Load gsMultiPatch
include("gsCore.jl")

# Load gsNurbs
include("gsNurbs.jl")

# Load gsHSplines
include("gsHSplines.jl")


end #module
