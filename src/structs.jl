struct Mutation
    location::Int64
    ancestral_value::Float64
    new_value::Float64
    fitness_diff::Float64
end

struct Genotype
    grn::Array{Float64,2}
    fixed_mutations::Array{Mutation,1}
end
