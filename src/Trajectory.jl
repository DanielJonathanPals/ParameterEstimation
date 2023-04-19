using CircularArrayBuffers
import Base: push!


"""
    Trajectory(data,Δt,capacity,space_dim)

# Fields
- 'data::CircularArrayBuffer{Float64}': The datapoints of a recorded trajectory
- 'Δt::Number': The time intervall between consecutive datapoints

The data can for instanciation of the struct also be given in the form of an Array{Float64}. In this case the capacity of the CircularArrayBuffer must be specified in an additional argument.
"""
struct Trajectory
    data::CircularArrayBuffer{Float64}
    Δt::Number
end


function Trajectory(data::Array{Float64},Δt::Number,capacity::Int)
    buff = CircularArrayBuffer{Float64}(size(data)[1],capacity)
    for i in 1:size(data)[2]
        push!(buff,data[:,i])
    end
    return Trajectory(buff,Δt)
end


function push!(traj::Trajectory,data)
    if size(data)[1] != traj.space_dim
        error("The dimension of the data that is attemted to be pushed into the Trajectory does not have the correct dimension")
    end
    for i in 1:size(data)[2]
        push!(traj.data,data[:,i])
    end
end


"""
    capacity(traj::Trajectory)

# Arguments
- 'traj::Trajectory'

# Returns
- The capacity of the CircularArrayBuffer containing the trajectory
"""
capacity(traj::Trajectory) = CircularArrayBuffers.capacity(traj.data)


"""
    space_dim(traj::Trajectory)

# Arguments
- 'traj::Trajectory'

# Returns
- The spacial dimension corresponding to the trajectory
"""
space_dim(traj::Trajectory) = size(traj.data)[1]