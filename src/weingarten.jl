using Iterators
include("mninner.jl")

character_symmetric_group(λ::Vector{Int}, perm_type::Vector{Int}) = mninner(λ, perm_type)
character_symmetric_group(λ::Vector{Int} = character_at_id(λ)

# function wg(perm_type::Vector{Int}, d::UInt)
#   n = sum(perm_type)
#   a
#   for plen=1:d
#     for λ in partitions(n, plen)
#       mapreduce
#     end
#   end
# end
