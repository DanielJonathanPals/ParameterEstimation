module ParameterEstimation

include("LinearEstimator.jl")
include("Trajectory.jl")

export(Trajectory)
export(push!)
export(capacity)
export(space_dim)
export(current_length)
export(linear_parameter_estimates)

end
