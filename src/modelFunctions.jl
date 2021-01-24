include("dictionaries.jl")

function generate_genotype(model::GeneModel, instruction::RandomGenome,
                           id::Int)
    """
    returns a geneotype, based on the specified size, connectivity, and
    interaction type.

    interaction_types:
    binary - all interacting genes have interaction of 1 or -1
    normal - normal distribution mean 0, sd 1
    uniform - uniform distribution over [-1,1]
    """
    genes = instruction.genomesize
    grn = rand(interaction_dict[instruction.interactiontype],genes,genes)

    num_unconnected = convert(Int64, floor((1-instruction.connectivity)*genes^2))
    locs = sample(1:genes^2,num_unconnected,replace=false,ordered=true)
    for loc in locs
        grn[id] = 0.0
    end
    model.individuals[id].grn = grn
    return 0
end
