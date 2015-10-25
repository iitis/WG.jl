function conjugate_partition(λ::Vector{Int})
  a = zeros(maximum(λ), length(λ))
  for i=1:length(λ)
    a[1:λ[i], i] = 1
  end
  vec(sum(a, 2))
end

function character_at_id(λ::Vector{Int})
  conj_λ = conjugate_partition(λ)
  d = 1
  for i=1:length(λ), j=1:λ[i]
    d *= (λ[i] - j + conj_λ[j] - i + 1)
  end
  factorial(sum(λ)) // d
end

function schur_poly(λ::Vector{Int}, d::Int)
  ret = Rational{Int}(1)
  λ0 = hcat(λ, zeros(d-length(λ)))
  for j=1:d, i=1:j-1
    ret *= (λ[i] - λ[j] +j - i)//(j - i)
  end
  ret
end

character_symmetric_group(λ::Vector{Int}, perm_type::Vector{Int}) = mninner(λ, perm_type)
character_symmetric_group(λ::Vector{Int} = character_at_id(λ)


