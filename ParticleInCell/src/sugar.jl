import LinearAlgebra
LinearAlgebra.dot(λ::Number, A::AbstractArray) = λ*A;
LinearAlgebra.dot(A::AbstractArray, λ::Number) = λ*A;
