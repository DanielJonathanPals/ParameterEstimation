using LinearAlgebra

include("Trajectory.jl")
include("Converter.jl")



function exp_λ(traj::Trajectory) # Returns the estimate of exp(-λ⋅Δt)
    if size(traj.data)[2] < 2
        error("The trajectory does not contain enought datapoints for the analysis")
    end
    n = current_length(traj)-1
    m_0 = sum(traj.data[:,1:end-1],dims=2) ./ n
    m_1 = sum(traj.data[:,2:end],dims=2) ./ n
    d_0 = traj.data[:,1:end-1] .- m_0
    d_1 = traj.data[:,2:end] .- m_1
    s_1 = sum([d_1[:,i]*d_0[:,i]' for i in 1:n])
    s_2 = sum([d_0[:,i]*d_0[:,i]' for i in 1:n])
    return s_1/s_2
end


λ(traj::Trajectory) = -log(exp_λ(traj)) ./ traj.Δt


function μ(traj::Trajectory; exp_λ_p::Union{Array,Nothing}=nothing)
    if exp_λ_p === nothing
        exp_λ_p = exp_λ(traj)
    end
    if size(exp_λ_p)[1] != size(exp_λ_p)[2]
        error("The matrix exp_λ_p is not a sqare matrix")
    end
    if space_dim(traj) != size(exp_λ_p)[1]
        error("The dimensions of exp_λ_p and the trajectory mismatch")
    end
    n = current_length(traj)-1
    d = space_dim(traj)
    v = sum(traj.data[:,2:end] - exp_λ_p*traj.data[:,1:end-1], dims=2) ./ n
    return (Matrix{Float64}(I,d,d) - exp_λ_p) \ v
end


function B(traj; exp_λ_p::Union{Array,Nothing}=nothing, μ_p::Union{Array,Nothing}=nothing)
    if exp_λ_p === nothing
        exp_λ_p = exp_λ(traj)
    end
    if μ_p === nothing
        μ_p = μ(traj,exp_λ_p = exp_λ_p)
    end
    n = current_length(traj)-1
    d = space_dim(traj)
    u = traj.data[:,2:end] - exp_λ_p*traj.data[:,1:end-1] .- (Matrix{Float64}(I,d,d) - exp_λ_p)*μ_p
    m = sum([u[:,i]*u[:,i]' for i in 1:n])
    if det(m) == 0
        error("The matrix under consideration is not invertible")
    end
    return n .* m^(-1)
end


function σ(traj; exp_λ_p::Union{Array,Nothing}=nothing, μ_p::Union{Array,Nothing}=nothing, B_p::Union{Array,Nothing}=nothing)    # Returns the Matrix σσᵀ
    if exp_λ_p === nothing
        exp_λ_p = exp_λ(traj)
    end
    if μ_p === nothing
        μ_p = μ(traj,exp_λ_p = exp_λ_p)
    end
    if B_p === nothing
        B_p = B(traj, exp_λ_p = exp_λ_p, μ_p = μ_p)
    end
    d = space_dim(traj)
    x = matrix_to_vec(B_p^(-1))
    λ = -log(exp_λ(traj)) ./ traj.Δt
    A = λ_to_A(λ)
    if det(Matrix{Float64}(I,d^2,d^2) - exp(-A .* traj.Δt)) == 0
        error("(1 - e⁻ᴬᵗ) is not invertible")
    end
    b = (Matrix{Float64}(I,d^2,d^2) - exp(-A  .* traj.Δt)) \ A * x
    return vec_to_matrix(b)
end


"""
    linear_parameter_estimates(traj::Trajectory)

This function returns the most likely values for the matracies λ and σσᵀ as well as for the vector μ 
    assuming that the measured trajectory follows an SODE of the form dx = -λx⋅dt-σdW

# Arguments
- 'traj::Trajectory': Trajectory to be analysed

# Returns
- 'λ': Estimate for λ
- 'μ': Estimate for μ
-  'σσᵀ': Estimate for σσᵀ
"""
function linear_parameter_estimates(traj::Trajectory)
    exp_λ_p = exp_λ(traj)
    μ_p = μ(traj; exp_λ_p = exp_λ_p)
    σ_p = σ(traj; exp_λ_p = exp_λ_p, μ_p = μ_p)
    λ_p = λ(traj)
    return λ_p, μ_p, σ_p
