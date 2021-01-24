"""
Individual Types
"""

"""
Any individual is a subtype of AbstractIndividual.
The individual should be a mutable struct with fields that store
data that need to be passed through steps of an experiment
"""
abstract type AbstractIndividual end


"""
A BasicIndividual is a minimal concrete type of the abstract type
AbstractIndividual, and it includes only the fields that are necessary
to complete an experiment.
"""
mutable struct BasicIndividual <: AbstractIndividual
    grn::Array{AbstractFloat,2}
end
