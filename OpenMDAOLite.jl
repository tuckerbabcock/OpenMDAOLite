"""
# OpenMDAOLite.jl

OpenMDAOLite.jl is a package that implements a minimal interface similar to
OpenMDAO, written in Julia.
"""

module OpenMDAOLite

# export Vector

include("vector.jl")
include("component.jl")
include("problem.jl")

end
