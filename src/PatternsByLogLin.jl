module PatternsByLogLin

using LogLinearModels
import MultivariateStats

export bestset, toppatterns, patterncoords

include("bestset.jl")
include("toppatterns.jl")
include("patterncoords.jl")

end # module
