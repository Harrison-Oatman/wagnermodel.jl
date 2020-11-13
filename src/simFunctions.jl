using ConfParser

function gen_config(conpath::String="none")
    """
    Takes in the path of a config file and produces sim_dict, a dictionary of
    experimental parameters.
    """
    if conpath == "none"
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
    sim_dict["walk_length"] = parse(Int64,retrieve(conf,"WALK_LENGTH"))
    sim_dict["mu_gain"]= parse(Float64,retrieve(conf, "GAIN_MUTATION_RATE"))
    sim_dict["mu_loss"] = parse(Float64,retrieve(conf, "LOSS_MUTATION_RATE"))
    sim_dict["mu_change"] = parse(Float64,retrieve(conf, "CHANGE_MUTATION_RATE"))
    sim_dict["sigma"] = parse(Float64,retrieve(conf, "SIGMA"))
    sim_dict["treatment"] = retrieve(conf,"TREATMENT")
    sim_dict["max_convergence_time"] = parse(Int64,retrieve(conf, "MAX_CONVERGENCE_TIME"))

    sim_dict["experiment_type"] = retrieve(conf,"EXPERIMENT_TYPE")
    sim_dict["selection_method"] = retrieve(conf,"SELECTION_METHOD")

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

config = gen_config()
