include("agents.jl")

"""
Model declarations
"""

mutable struct GeneModel{A <: AbstractIndividual}
    individuals::Dict{Int, A}
    # maxid::Base.RefValue{Int64}
end
