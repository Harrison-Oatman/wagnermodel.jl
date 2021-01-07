using Parameters

"""
Types and Structs
"""
abstract type Instruction end

@with_kw struct SetPlayers <: Instruction
    name::String = "players"
    writeResults::Bool = true
    playerone::Function
    playertwo::Function
end

@with_kw struct Repeater <: Instruction
    itinerary::Vector{Instruction}
    iterates::Int
end

@with_kw struct PlayGame <: Instruction
    name::String = "playedGame"
    writeResults::Bool = true
end

@with_kw struct RandomBoard <: Instruction
    name::String = "board"
    choices::Int = 2
    writeResults::Bool = true
    symmetric::Bool = false
end

@with_kw struct CustomBoard <: Instruction
    board::Array
    name::String = "board"
    writeResults::Bool = true
end

@with_kw struct Trials <: Instruction
    itineraries::Vector{Vector{Instruction}}
end

"""
Instruction Handlers
"""
function manager(instruction::Instruction, inputDict::Dict, config::Dict)
    """
    'hub' function, which directs the call to the function corresponding
    to the given function
    """
    if haskey(inputDict,"dir") == false
        println("here")
        currentdir = string(@__DIR__)
        if isdir(currentdir*"/tests") == false
            mkdir(currentdir*"/tests")
        end
        inputDict["dir"] = currentdir*"/tests"
    end

    functionDict = Dict(
        Repeater => repeater,
        PlayGame => playgame,
        RandomBoard => randomboard,
        CustomBoard => customboard,
        Trials => trials,
        SetPlayers => setplayers
        )

    instructionType = typeof(instruction)
    outputDict = functionDict[instructionType](instruction, inputDict, config)
    if hasfield(instructionType,:writeResults)
        if instruction.writeResults
            save(outputDict,outputDict["dir"],instruction.name)
        end
    end
    return outputDict
end

function repeater(instruction::Instruction,inputDict::Dict, config::Dict)
    outputDict = inputDict
    basedir = inputDict["dir"]
    for i in 1:instruction.iterates
        for step in instruction.itinerary
            outputDict["dir"] = basedir*"/iterate_"*string(i)
            if isdir(outputDict["dir"]) == false
                mkdir(outputDict["dir"])
            end
            outputDict = manager(step,outputDict,config)
        end
    end
    return outputDict
end

function playgame(instruction::Instruction,inputDict::Dict,config::Dict)
    """
    requires:
        board: a matrix of choices
        player_one, player_two: two functions which accept the dictionary of
            available information, and return a choice
    """
    # if a KeyError occurs, the previous instruction may be incompatible
    # with the currect instruction.
    board = inputDict["board"]

    player_one = inputDict["player_one"]
    player_two = inputDict["player_two"]

    prev_choice_one = 0
    prev_choice_two = 0

    try
        prev_choice_one = inputDict["prev_choice_one"]
        prev_choice_two = inputDict["prev_choice_two"]
    catch KeyError
    end

    playerOneDict = Dict(
    "board" => board,
    "my_prev_choice" => prev_choice_one,
    "their_prev_choice" => prev_choice_two
    )
    playerTwoDict = Dict(
    "board" => [board[1],board[2]],
    "my_prev_choice" => prev_choice_two,
    "their_prev_choice" => prev_choice_one
    )

    i = player_one(playerOneDict)
    j = player_two(playerTwoDict)

    outputDict = inputDict
    outputDict["prev_choice_one"] = i
    outputDict["prev_choice_two"] = j

    println("Player One Score: ",board[1][i,j],". Player Two Score: ",board[2][j,i],".")

    return outputDict
end

function randomboard(instruction::Instruction,inputDict::Dict,config::Dict)
    outputDict = inputDict
    n = instruction.choices
    board = [zeros(n,n),zeros(n,n)]
    for i in 1:length(board[1])
        board[1][i] = rand(1:10)
    end
    if instruction.symmetric
        board[2] = deepcopy(board[1])
    else
        for i in 1:length(board[1])
            board[2][i] = rand(1:10)
        end
    end
    outputDict["board"] = board
    outputDict["choices"] = length(board[1][1,:])
    return outputDict
end

function customboard(instruction::Instruction,inputDict::Dict,config::Dict)
    outputDict = inputDict
    outputDict["board"] = instruction.board
    outputDict["choices"] = length(outputDict["board"][1][1,:])
    return outputDict
end

function setplayers(instruction::Instruction,inputDict::Dict,config::Dict)
    outputDict = inputDict
    outputDict["player_one"] = instruction.playerone
    outputDict["player_two"] = instruction.playertwo
    outputDict["prev_choice_one"] = 0
    outputDict["prev_choice_two"] = 0
    return outputDict
end

function trials(instruction::Instruction,inputDict::Dict,config::Dict)
end

"""
Helper Functions
"""
function save(dict::Dict, dir::String, filename::String)
    f = open(dir*"/"*filename*".dat","w")
    for (i,key) in enumerate(sort(collect(keys(dict))))
        if key != "dir"
            if i != 1
                write(f, '\n')
            end
            write(f,string(key)*' '*string(dict[key]))
        end
    end
    close(f)
end

"""
Strategy Functions
"""
function A100(infoDict)
    return 1
end

function B100(infoDict)
    return 2
end

function greedy(infoDict)
    idx = argmax(infoDict["board"][1])
    return idx[1]
end


config = Dict()
inputDict = Dict()

players = SetPlayers(playerone = greedy, playertwo = B100)
genboard = RandomBoard(symmetric = true)
rungame = PlayGame()

trial = Repeater([players,genboard,rungame],100)

outputDict = manager(trial,inputDict,config)
