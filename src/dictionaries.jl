using Distributions

interaction_dict = Dict()

interaction_dict["binary"] = [-1.0,1.0]
interaction_dict["normal"] = Normal()
interaction_dict["uniform"] = Uniform(-1.0,1.0)
