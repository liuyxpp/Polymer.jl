@testset "config.jl: load_config" begin
    configfile = joinpath(@__DIR__, "ABS.yml")
    config_abs = load_config(configfile, top="system")
    @test config_abs isa PolymerSystemConfig
end

@testset "config.jl: save_config" begin
    configfile = joinpath(@__DIR__, "ABS.yml")
    config_abs = load_config(configfile, top="system")
    configfile = joinpath(@__DIR__, "ABS_save.yml")
    save_config(configfile, config_abs)
    # The config file is saved without top level key.
    config_abs_reload = load_config(configfile)
    @test config_abs_reload isa PolymerSystemConfig
end

nothing