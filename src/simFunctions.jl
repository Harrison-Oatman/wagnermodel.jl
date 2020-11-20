using ConfParser
using Random
using GLM
using Distributions

include("structs.jl")

function generate_config(conpath::String="none")
    """
    Takes in the path of a config file and produces sim_dict, a dictionary of
    experimental parameters.
    """
    if conpath == "none" # DEFAULT config is sim.cfg, currently in src
        print("using default config...")
        conpath = string(@__DIR__) * "\\sim.cfg"
    end
    conf = ConfParse(conpath)
    parse_conf!(conf)

    sim_dict = Dict{String,Any}()
    sim_dict["random_seed"] = parse(Int64,retrieve(conf, "RANDOM_SEED"))
    sim_dict["genes"] = parse(Int64,retrieve(conf, "GENES"))
    sim_dict["ancestor_connectivity"] = parse(Float64,retrieve(conf, "ANCESTOR_CONNECTIVITY"))
    sim_dict["interaction_type"] = retrieve(conf,"INTERACTION_TYPE")
    sim_dict["mu_gain"]= parse(Float64,retrieve(conf, "GAIN_MUTATION_RATE"))
    sim_dict["mu_loss"] = parse(Float64,retrieve(conf, "LOSS_MUTATION_RATE"))
    sim_dict["mu_change"] = parse(Float64,retrieve(conf, "CHANGE_MUTATION_RATE"))
    sim_dict["treatment"] = retrieve(conf,"TREATMENT")
    sim_dict["max_convergence_time"] = parse(Int64,retrieve(conf, "MAX_CONVERGENCE_TIME"))

    sim_dict["experiment_type"] = retrieve(conf,"EXPERIMENT_TYPE")
    sim_dict["selection_method"] = retrieve(conf,"SELECTION_METHOD")
    sim_dict["sigma"] = parse(Float64,retrieve(conf, "SIGMA"))

    sim_dict["population_size"] = parse(Int64,retrieve(conf, "POPULATION_SIZE"))

    sim_dict["walk_length"] = parse(Int64,retrieve(conf,"WALK_LENGTH"))
    sim_dict["reverse_signal"] = parse(Bool,retrieve(conf,"REVERSE_SIGNAL"))

    sim_dict["parent_folder"] = retrieve(conf,"PARENT_FOLDER")

    # for the purpose of testing
    sim_dict["test_param"] = parse(Bool,retrieve(conf,"TEST_PARAM"))

    return sim_dict
end

function assess_stability(test_grn::Array{Float64,2},test_S_0::Array{Float64,1},config::Dict)
    """
    Goal: This function measures whether a given GRN results in a stable final
    gene state vector for a given initial gene state vector.

    Inputs:
        1. test_grn: gene regulatory network matrix to be tested
        2. test_S_0: vector of initial gene states

    Outputs:
        1. boolean variable on whether the grn is stable
        2. the vector of final gene states
        3. the number of steps it took to reach a stable state
    """

    function filter_function(x::Float64)
        """
        Goal: This function normalizes the value of a gene's state to either
        0.0 or 1.0.
        """

        if x > 0
            return 1.0
        else
            return 0.0
        end
    end

    S_prev = deepcopy(test_S_0)
    for t=1:config["max_convergence_time"]
        S_next_temp = test_grn*S_prev
        S_next = map(filter_function,S_next_temp)

        #If stable, return stability, final gene state, and number of iterations
        #it took to reach state
        if S_next == S_prev
            return true, S_next, t
        else
            S_prev = deepcopy(S_next)
        end
    end

    #If not stable, return 0 vector for stable state
    return false,zeros(length(test_S_0)),-1
end

function generate_genotype(genes::Int,interaction_type::String,
    ancestor_connectivity::Real)
    """
    returns a geneotype, based on the specified size, connectivity, and
    interaction type.

    interaction_types:
    binary - all interacting genes have interaction of 1 or -1
    normal - normal distribution mean 0, sd 1
    uniform - uniform distribution over [-1,1]
    """

    if interaction_type == "binary"
        grn = rand([-1.0,1.0],genes,genes)
    elseif interaction_type == "normal"
        grn = rand(Normal(),genes,genes)
    elseif interaction_type == "uniform"
        grn = rand(Uniform(-1,1),genes,genes)
    else
        throw(error("unknown interaction type"))
    end

    num_unconnected = convert(Int64, (1-ancestor_connectivity)*genes^2)
    idx = sample(1:genes^2,num_unconnected,replace=false,ordered=true)
    for id in idx
        grn[id] = 0.0
    end
    return Genotype(grn, Array{Mutation,1}())
end

function generate_genotype(genes::Int,sample_distribution::Sampleable,
    ancestor_connectivity::Real)
    """
    returns a genotype based on the specified size and connectivity, using the
    speicified sampling method.
    """
    grn = rand(sample_distribution,genes,genes)

    num_unconnected = convert(Int64, (1-ancestor_connectivity)*genes^2)
    idx = sample(1:genes^2,num_unconnected,replace=false,ordered=true)
    for id in idx
        grn[id] = 0.0
    end
    return Genotype(grn, Array{Mutation,1}())
end

function generate_genotype(config::Dict)
    return generate_genotype(config["genes"],config["interaction_type"],
        config["ancestor_connectivity"])
end

function generate_genotype(config::Dict,sample_distribution::Sampleable)
    return generate_genotype(config["genes"],sample_distribution,
        config["ancestor_connectivity"])
end

function mutate_genotype!(genotype::Genotype,mutation_type::String,interaction_type::String)
    """
    Mutates a grn with using the specified mutation type
    """
    mutation_found = false
    grn = genotype.grn
    genes = size(grn,1)
    loc = -1
    oldval = 0
    newval = 0
    while mutation_found == false
        loc = sample(1:genes^2)
        if mutation_type == "loss" && grn[loc] != 0.0
            newval = 0
        elseif mutation_type == "change" && grn[loc] != 0.0
            newval = grn[loc]*-1
        elseif mutation_type == "gain" && grn[loc] == 0.0
            if interaction_type == "binary"
                newval = rand([-1.0,1.0])
            elseif interaction_Type == "normal"
                newval = rand(Normal())
            elseif interaction_type == "uniform"
                newval = rand(Uniform(-1.0,1.0))
            else
                throw(error("unsupported interaction type"))
            end
        else
            continue
        end

        oldval = grn[loc]
        mutation_found = true
    end
    grn[loc] = newval
    push!(genotype.fixed_mutations, Mutation(loc, oldval, newval))
    return genotype
end

function mutate_genotype!(genotype::Genotype,mutation_type::String,sample_distribution::Sampleable)
    """
    Mutates a grn with using the specified mutation type, with the specified
    distribution
    """
    mutation_found = false
    grn = genotype.grn
    genes = size(grn,1)
    loc = -1
    oldval = 0
    newval = 0
    while mutation_found == false
        loc = sample(1:genes^2)
        if mutation_type == "loss" && grn[loc] != 0.0
            newval = 0
        elseif mutation_type == "change" && grn[loc] != 0.0
            # currently just inverts interaction
            newval = grn[loc]*-1
        elseif mutation_type == "gain" && grn[loc] == 0.0
            newval = rand(sample_distribution)
        else
            continue
        end
        oldval = grn[loc]
        mutation_found = true
    end
    grn[loc] = newval
    push!(genotype.fixed_mutations, Mutation(loc, oldval, newval))
    return genotype
end

function rand_mutate!(genotype::Genotype,config::Dict)
    genes = size(genotype.grn,1)
    empty = count(i->(i==0.0),genotype.grn)
    println(empty)
    gains = rand(Binomial(empty,config["mu_gain"]))
    losses = rand(Binomial(genes^2 - empty,config["mu_loss"]))
    changes = rand(Binomial(genes^2 - empty,config["mu_change"]))
    println(gains,losses,changes)
    for i in 1:gains
        mutate_genotype!(genotype,"gain",config["interaction_type"])
    end
    for i in 1:losses
        mutate_genotype!(genotype,"loss",config["interaction_type"])
    end
    for i in 1:changes
        mutate_genotype!(genotype,"change",config["interaction_type"])
    end
    return gains + losses + changes
end

function generate_population(config::Dict,init_conditions::Array{Float64,1},single_ancestor=true)
    """
    Generates a population given the specifications of the config file.
    Each member of the population will be viable under the init_conditions specified
    """
    population = Population(Array{Tuple{Genotype,Int},1}())
    popsize = config["population_size"]
    unique_individuals = single_ancestor ? 1 : popsize
    for individual in 1:unique_individuals
        viable = false
        genotype = nothing
        while viable == false
            genotype = generate_genotype(config)
            viable = assess_stability(genotype.grn,init_conditions,config)[1]
        end
        push!(population.individuals,(genotype,single_ancestor ? popsize : 1))
    end
    return population
end

config = generate_config()
grn = generate_genotype(config)
init_cond = rand([-1.0,1.0],config["genes"])
oldgrn = deepcopy(grn)
rand_mutate!(grn, config)
population = generate_population(config,init_cond,false)
