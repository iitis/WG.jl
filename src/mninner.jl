function mninner(R::Vector{Int}, v::Vector{Int}, t::Int=1)
  if t>=length(v)
    return 1
  end
  Χ = 0
  σ = 1
  pow = length(find(x->x==0, R))
  σ *= -1^pow
  for i=1:length(R)-v[t]+1
    if R[i] != R[i+v[t]-1] σ *= -1; end
    if i+v[t] <= length(R) && R[i] == 1 && R[i+v[t]] == 0
      R[i], R[i+v[t]] = R[i+v[t]], R[i]
      Χ += σ * mninner(R, v, t+1)
      R[i], R[i+v[t]] = R[i+v[t]], R[i]
    end
  end
  Χ
end

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
  for i=1:length(λ)
    for j=1:λ[i]
      d *= (λ[i] - j + conj_λ[j] - i + 1)
    end
  end
  factorial(sum(λ)) // d
end


