using WagnerModel
using Test

@testset "WagnerModel.jl" begin
    config = gen_config()
    @test config["test_param"]

    grn = [[1.0,0.0],[0.0,1.0]]
    grn = hcat(grn...)
    S_0 = [1.0,0.0]
    @test assess_stability(grn, S_0, config)[1]
    @test assess_stability(grn, S_0, config)[2] == [1.0,0.0]
    @test assess_stability(grn, S_0, config)[3] == 1
end
