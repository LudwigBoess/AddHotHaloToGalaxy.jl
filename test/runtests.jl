using AddHotHaloToGalaxy, Test

@testset "Parameter File" begin
    parfile = "../example_parameter_file.yml"

    par = AddHotHaloToGalaxy.read_parameter_file(parfile)

    @test par["G"] ≈ 6.672e-8
end