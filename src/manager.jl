include("instructions.jl")

function manager(instruction::T, model::GeneModel) where {T <: IndividualInstruction}
    if hasfield(T, :tag)
        tag = instruction.tag
        for id in collect(keys(model.individuals))
            if tag = nothing
                instruction.name(model, instruction, id)
            end
        end
    end
end

instruction = RandomGenome(2,"binary",1.0, nothing, nothing, false, generate_genotype)

individuals = Dict(1,BasicIndividual(zeros(2,2)))

model = GeneModel(individuals)

manager(instruction, model)
