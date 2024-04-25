module Gismo

# Get the path to the library
include("src/Gismo_path.jl")

# Forward declaration of structs
include("src/Declarations.jl")

# Load gsMatrix
include("src/gsMatrix.jl")

# Load gsMultiPatch
include("src/gsCore.jl")

# Load gsNurbs
include("src/gsNurbs.jl")

# Load gsHSplines
include("src/gsHSplines.jl")


end #module
