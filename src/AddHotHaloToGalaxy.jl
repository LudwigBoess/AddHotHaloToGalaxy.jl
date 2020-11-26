module AddHotHaloToGalaxy

    include("io.jl")
    include("construct_halo.jl")
    include("cut_ids_from_disk.jl")

    export add_hot_halo_to_galaxy

    function add_hot_halo_to_galaxy(parameterfile::String; verbose::Bool=true)

        if verbose
            t1_total = time_ns()
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

        pos_halo, vel_halo, u_halo, m_halo, rho_halo = construct_hot_halo(par)

        if par["verbose"]
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

        if par["verbose"]
            r_dummy = @. sqrt(pos_halo[:,1]^2 + pos_halo[:,2]^2 + pos_halo[:,3]^2 )
            @info "r_max = $(maximum(r_dummy))"
        end

        if par["verbose"]
            @info "Selecting IDs of Gas particles to be cut out"
            t1 = time_ns()
        end

        cut_ids = cut_ids_from_halo(pos_halo, rho_halo, par)

        if par["verbose"]
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

        npart_before = size(u_halo,1)

        # only keep relevant particles
        pos_halo = pos_halo[cut_ids,:]
        vel_halo = vel_halo[cut_ids,:]
        rho_halo = rho_halo[cut_ids]
        u_halo   = u_halo[cut_ids]
        m_halo   = m_halo[cut_ids]

        npart_after = size(u_halo,1)
        n_cut = size(findall(cut_ids .== false),1)

        if ( npart_before - npart_after ) != n_cut
            error("Wrrong number of particles cut!
                   nbefore = $npart_before
                   nafter  = $npart_after
                   n_cut   = $n_cut
                   npart_before - npart_after = $(npart_before - npart_after)")
        end 

        

        if par["verbose"]
            @info "Writing to file"
            t1 = time_ns()
        end

        write_to_file(pos_halo, vel_halo, u_halo, m_halo, par)

        if par["verbose"]
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

        # Output total time
        if verbose
            t2_total = time_ns()
            println()
            @info "Total runtime $(output_time(t1_total,t2_total)) s"
        end

    end

end # module