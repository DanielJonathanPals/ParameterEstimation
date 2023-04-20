module ParameterEstimation

export Trajectory
export push!
export capacity
export space_dim
export current_length
export linear_parameter_estimates

include("LinearEstimator.jl")
include("Trajectory.jl")

end
