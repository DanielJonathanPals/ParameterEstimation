using BlockArrays
using LinearAlgebra


function vec_to_matrix(vec::Array)
    vec = reshape(vec,length(vec),1)
    try
        d = Int64(sqrt(length(vec)))
    catch
        error("The number of entries of the vector must be a square number")
    end
    d = convert(Int64,d)
    return reshape(vec,d,d)
end


function matrix_to_vec(matrix::Array)
    if size(matrix)[1] != size(matrix)[2]
        error("The matrix in the argument of the function must be a square matrix")
    end
    return reshape(matrix,size(matrix)[1]^2,1)
end


function λ_to_A(λ)
    if size(λ)[1] != size(λ)[2]
        error("λ must be quadratic")
    end
    d = size(λ)[1]
    x = [d for i in 1:d]
    A = BlockArray{Complex{Float64}}(zeros(d^2,d^2), x, x)
    for i in 1:d, j in 1:d
        A[Block(i),Block(j)] = (λ[j,i]) .* Matrix{Float64}(I,d,d)
        if i == j
            A[Block(i),Block(j)] += λ
        end
    end
    return Array(A)
end