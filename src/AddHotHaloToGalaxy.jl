module AddHotHaloToGalaxy

    include("io.jl")
    include("construct_halo.jl")

    export add_hot_halo_to_galaxy

    function add_hot_halo_to_galaxy(parameterfile::String; verbose::Bool=true)

        if verbose
            @info "Reading parameter file"
            t1 = time_ns()
        end

        par = read_parameter_file(parameterfile)

        if verbose
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end
        

        if par["verbose"]
            @info "Sampling Gas halo"
            t1 = time_ns()
        end

        pos_halo, vel_halo, u_halo, m_halo = construct_hot_halo(par)

        if par["verbose"]
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

    end

end # module