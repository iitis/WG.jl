using Iterators
include("utils.jl")

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

# function wg(perm_type::Vector{Int}, d::UInt)
#   n = sum(perm_type)
#   a
#   for plen=1:d
#     for λ in partitions(n, plen)
#       mapreduce
#     end
#   end
# end
