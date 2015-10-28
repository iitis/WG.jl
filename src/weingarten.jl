using Iterators
import Base.Rational
include("utils.jl")

function mninner(R::Vector{Int}, v::Vector{Int}, t::Int=1)
  if t>=length(v)
    return 1
  end
  Χ = 0
  σ = 1
  pow = length(find(x->x==0, R[1:min(v[t], length(R))]))
  σ *= (-1)^pow
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

function wg(perm_type::Vector{Int}, d::Int)
  sort!(perm_type, rev=true)
  n = sum(perm_type)
  r = Rational(0)
  for λ in partitions_up_to_length(n, d)
    r += character_symmetric_group(λ)^2 * character_symmetric_group(λ, perm_type) / schur_poly(λ, d)
  end
  r = 1//factorial(n)^2 * r
end
