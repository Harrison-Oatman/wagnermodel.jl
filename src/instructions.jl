using Parameters
include("modelFunctions.jl")
"""
An instruction is an object used to tell the experiment manager when and how
to run a step of the experiment. all instructions are subtypes of this
abstract type
"""
abstract type AbstractInstruction end

"""
An individual instruction tells the experiment manager to run a step of the
experiment that acts on one or more individuals.
--
An individual instruction can include the field "tag," which takes a Symbol.
When the individual instruction
"""
abstract type IndividualInstruction <: AbstractInstruction end
abstract type ModelInstruction <: AbstractInstruction end

@with_kw struct RandomGenome <: IndividualInstruction
    genomesize::Int
    interactiontype::String
    connectivity::Real
    # filter parameters
    tag::Union{Symbol,Nothing} = nothing
    value::Any = nothing
    exclude::Bool = false
    # associated function
    name::function = generate_genotype
end
