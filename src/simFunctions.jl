using ConfParser
using Random
using GLM
using Distributions
using StatsBase
using Plots

include("dictionaries.jl")
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

function assess_stability(test_grn::Array{Float64,2},test_S_0::Array{Float64,1},max_convergence_time::Int)
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

        if x > 0.5
            return 1.0
        else
            return 0.0
        end
    end

    S_prev = deepcopy(test_S_0)
    for t=1:max_convergence_time
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

    grn = rand(interaction_dict[interaction_type],genes,genes)

    num_unconnected = convert(Int64, floor((1-ancestor_connectivity)*genes^2))
    idx = sample(1:genes^2,num_unconnected,replace=false,ordered=true)
    for id in idx
        grn[id] = 0.0
    end
    return Genotype(grn, Array{Mutation,1}())
end

function generate_condition(genes::Int,connectivity::Float64,deterministic::Bool=false)
    if deterministic
        n = floor(genes*connectivity)
        condition = [0.0 for _ in 1:(genes-n)]
        push!(condition,[1.0 for _ in 1:n]...)
        return sample(condition,10,replace=false)
    else
        return sample([0.0,1.0],Weights([connectivity,1-connectivity]))
    end
end

# function generate_genotype(genes::Int,sample_distribution::Sampleable,
#     ancestor_connectivity::Real)
#     """
#     returns a genotype based on the specified size and connectivity, using the
#     speicified sampling method.
#     """
#     grn = rand(sample_distribution,genes,genes)
#
#     num_unconnected = convert(Int64, (1-ancestor_connectivity)*genes^2)
#     idx = sample(1:genes^2,num_unconnected,replace=false,ordered=true)
#     for id in idx
#         grn[id] = 0.0
#     end
#     return Genotype(grn, Array{Mutation,1}())
# end

function hamming_distance(A::Array{Float64,1},B::Array{Float64,1})
    distance = 0.0
    for i in 1:length(A)
        if A[i] != B[i]
            distance += 1.0
        end
    end
    return distance
end

function gen_grn_from_states(S_0::Array{Float64,1},target_eq::Array{Float64,1},
    interaction_type::String,ancestor_connectivity::Real,max_convergence_time::Int)
    """
    generates a genotype given an initial state, that converges to the target state
    """


    genes = length(S_0)
    @assert length(target_eq) == genes

    grn = nothing
    stable = false
    current_eq = target_eq
    while stable==false
        grn = generate_genotype(genes, interaction_type, ancestor_connectivity).grn
        stable, current_eq, t = assess_stability(grn,S_0,max_convergence_time)
    end

    p_change = 0.5
    p_gain = ancestor_connectivity*0.5
    p_loss = (1.0 - ancestor_connectivity)*0.5

    grn_found = false
    distance = hamming_distance(current_eq,target_eq)
    total_steps = 0
    while grn_found == false
        total_steps += 1
        new_grn = simple_random_mutate(grn,[p_change,p_gain,p_loss],interaction_type,
        sample(1:10))
        stable, next_eq, t = assess_stability(new_grn,S_0,max_convergence_time)
        if stable
            next_distance = hamming_distance(next_eq,target_eq)
            accept = false
            if next_distance < distance
                accept = true
            elseif next_distance == distance
                accept = sample([true,false])
            else
                p_accept = (distance/next_distance)^2
                accept = sample([true,false],Weights([p_accept,1-p_accept]))
            end
            if accept
                grn = new_grn
                current_eq = next_eq
                distance = next_distance
            end
        end
        if distance == 0.0
            grn_found = true
        end
    end
    println(total_steps)
    return grn
end

function simple_random_mutate(oldgrn::Array{Float64,2},mutation_probs::Array{Float64,1},
        interaction_type::String,mutations::Int)
    """
    finds a new mutation without altering the old grn, and without updating the
    mutation record of any genotype. iterating this function will produce a
    steady-state connectivity of p_gain/p_change
    """
    grn = deepcopy(oldgrn)
    genes = size(grn,1)

    for _ in 1:mutations
        mutation_found = false
        loc = -1
        newval = 0
        while mutation_found == false
            loc = sample(1:genes^2)
            mutation_type = sample(["change","gain","loss"],Weights(mutation_probs))

            if mutation_type == "change" && grn[loc] != 0.0
                # currently only negates the interaction, for nonbinary interaction
                # should this be a reroll?
                newval = grn[loc]*-1
            elseif mutation_type == "gain" && grn[loc] == 0.0
                newval = rand(interaction_dict[interaction_type])
            elseif mutation_type == "loss" && grn[loc] != 0.0
                newval = 0
            else
                continue
            end
            mutation_found = true
        end
        grn[loc] = newval
    end
    return grn
end

function mutate_genotype!(genotype::Genotype,mutation_type::String,interaction_type::String)
    """
    makes a single mutation on a grn of the specified type, at a random
    valid location
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
            # currently only negates the interaction, for nonbinary interaction
            # should this be a reroll?
            newval = grn[loc]*-1
        elseif mutation_type == "gain" && grn[loc] == 0.0
            newval = rand(interaction_dict[interaction_type])
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

function pi(X, duplicates)
    output = Array{Float64,1}()
    for i in 1:length(X)
        push!(output, X[i])
        if i in duplicates
            push!(output, X[i])
        end
    end
    return output
end

function pi_mat(A, duplicates)
    new_A = []
    for i in 1:length(A[1,:])
        push!(new_A,pi(A[i,:],duplicates))
        if i in duplicates
            push!(new_A,pi(A[i,:],duplicates))
        end
    end
    return hcat(new_A...)
end

# function mutate_genotype!(genotype::Genotype,mutation_type::String,sample_distribution::Sampleable)
#     """
#     Mutates a grn with using the specified mutation type, with the specified
#     distribution
#     """
#     mutation_found = false
#     grn = genotype.grn
#     genes = size(grn,1)
#     loc = -1
#     oldval = 0
#     newval = 0
#     while mutation_found == false
#         loc = sample(1:genes^2)
#         if mutation_type == "loss" && grn[loc] != 0.0
#             newval = 0
#         elseif mutation_type == "change" && grn[loc] != 0.0
#             # currently just inverts interaction
#             newval = grn[loc]*-1
#         elseif mutation_type == "gain" && grn[loc] == 0.0
#             newval = rand(sample_distribution)
#         else
#             continue
#         end
#         oldval = grn[loc]
#         mutation_found = true
#     end
#     grn[loc] = newval
#     push!(genotype.fixed_mutations, Mutation(loc, oldval, newval))
#     return genotype
# end

function rand_mutate!(genotype::Genotype,config::Dict)
    """
    determines the number of mutations for an individual using the given
    probabilities of interaction loss, gain, and change
    """
    genes = size(genotype.grn,1)
    empty = count(i->(i==0.0),genotype.grn)

    gains = rand(Binomial(empty,config["mu_gain"]))
    losses = rand(Binomial(genes^2 - empty,config["mu_loss"]))
    changes = rand(Binomial(genes^2 - empty,config["mu_change"]))

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
            genotype = generate_genotype(config["genes"],config["interaction_type"],
                        config["ancestor_connectivity"])
            viable = assess_stability(genotype.grn,init_conditions,config)[1]
        end
        push!(population.individuals,(genotype,single_ancestor ? popsize : 1))
    end
    return population
end

function wagner_experiment(genes,k_values,trials,subtrials)
    connectivity = 1.0
    S_0_list = []
    S_eq_list = []
    grn_list = []
    for i in 1:trials
        push!(S_0_list,generate_condition(genes,connectivity,true))
        push!(S_eq_list,generate_condition(genes,connectivity,true))
        push!(grn_list,gen_grn_from_states(S_0_list[i],S_eq_list[i],"binary",connectivity,100))
    end

    similar = zeros(length(k_values))
    stables = zeros(length(k_values))
    for (i,k) in enumerate(k_values)
        for trial in 1:trials
            for subtrial in 1:subtrials
                duplicates = sample(1:genes,k,replace=false,ordered=true)
                # println(duplicates)
                new_grn = pi_mat(grn_list[trial],duplicates)
                stable, new_S_eq, t = assess_stability(new_grn,pi(S_0_list[trial],duplicates),100)
                # print(stable)
                if stable
                    stables[i] += 1.0
                    # println(new_S_eq,pi(S_eq_list[trial],duplicates))
                    if new_S_eq == pi(S_eq_list[trial],duplicates)
                        similar[i] += 1.0
                    end
                end
            end
        end
        similar[i] = 1 - similar[i]/stables[i]
    end
    print(stables)
    return similar
end

# code for testing
config = generate_config()

S_0 = generate_condition(10,0.5,true)
S_eq = generate_condition(10,0.5,true)
grn = gen_grn_from_states(S_0,S_eq,config["interaction_type"], config["ancestor_connectivity"],20)
# grn = generate_genotype(config["genes"],config["interaction_type"],
#     config["ancestor_connectivity"])
# init_cond = rand([-1.0,1.0],config["genes"])
# oldgrn = deepcopy(grn)
# rand_mutate!(grn, config)
# population = generate_population(config,init_cond,false)
m =  pi_mat(grn,[1,2,3])


results = wagner_experiment(10,[0,1,2,3,4,5,6,7,8,9,10],100,100)
plot([0,1,2,3,4,5,6,7,8,9,10],results)
