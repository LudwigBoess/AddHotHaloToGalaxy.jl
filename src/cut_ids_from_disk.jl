using TriangularShapedCloudInterpolation
using ProgressMeter
using Base.Threads 

function calculate_COM(x)
    com = zeros(3)
    for i = 1:3
        com[i] = sum(x[:,i]) / length(x[:,i])
    end
    for i = 1:3
        x[:,i] .-= com[i]
    end
    return x
end

function run_tsc(pos_disk, rho, resx, resy, resz, verbose)

    pos_disk_tsc = zeros(length(rho),3)

    minx = minimum(pos_disk[:,1])
    maxx = maximum(pos_disk[:,1]) .* (1.0+1.e-6)
    dx   = -(minx - maxx) / resx
    pos_disk_tsc[:,1] = ( pos_disk[:,1] .- minx ) ./ dx 

    miny = minimum(pos_disk[:,2])
    maxy = maximum(pos_disk[:,2]) .* (1.0+1.e-6)
    dy   = -(miny - maxy) / resy
    pos_disk_tsc[:,2] = ( pos_disk[:,2] .- miny ) ./ dy 

    minz = minimum(pos_disk[:,3])
    maxz = maximum(pos_disk[:,3]) .* (1.0+1.e-6)
    dz   = -(minz - maxz) / resz
    pos_disk_tsc[:,3] = ( pos_disk[:,3] .- minz ) ./ dz

    if verbose
        @info "    Running TSC"
        t1 = time_ns()
    end

    tsc = TSCInterpolation(rho, pos_disk_tsc[:,1], resx, 
                                pos_disk_tsc[:,2], resy, 
                                pos_disk_tsc[:,3], resz, 
                            average=true)

    if verbose
        t2 = time_ns()
        @info "    Done!"
        @info "    Took $(output_time(t1,t2)) s"
    end

    return tsc, minx, miny, minz, dx, dy, dz

end


@inline function find_ids_from_density_criterion(   keep_id, pos_halo, rho_halo, rho_tsc, 
                                                    resx, resy, resz,
                                                    dx, dy, dz,
                                                    minx, miny, minz
                                                )

    @showprogress "Filtering IDs..." for ix = 1:resx, iy = 1:resy, iz = 1:resz

        lx = (ix-1)*dx + minx
        rx = lx + dx
        ly = (iy-1)*dy + miny
        ry = ly + dy
        lz = (iz-1)*dz + minz
        rz = lz + dz

        # @inbounds id = findall( (lx .< pos_halo[:,1] .< rx) .&
        #                         (ly .< pos_halo[:,2] .< ry) .&
        #                         (lz .< pos_halo[:,3] .< rz) )

        
        # idid = findall( rho_halo[id] .< 0.1*rho_tsc[ix, iy, iz] )

        # keep_id[idid] .= false

        #@threads 
        # for i = 1:length(id)
        #     if rho_halo[id[i]] < 0.1*rho_tsc[ix, iy, iz]
        #         @inbounds keep_id[id[i]] = false
        #     end
        # end

        @threads for i = 1:length(rho_halo)
            @inbounds if (  (lx < pos_halo[i,2] < rx) &
                            (ly < pos_halo[i,1] < ry) &
                            (lz < pos_halo[i,3] < rz) )

                if rho_halo[i] < 0.1*rho_tsc[ix, iy, iz]
                    keep_id[i] = false
                end
            end
        end
    end

    return keep_id

end

function cut_ids_from_halo(pos_halo, rho_halo, par)

    if par["verbose"]
        @info "  Reading galaxy data"
        t1 = time_ns()
    end

    # read in data from galaxy
    pos_disk, rho_disk, hsml_disk = read_pos_rho_hsml_galaxy(par["input_snap"], par["ic_format"])

    if par["verbose"]
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

    if par["verbose"]
        @info "  Shifting to center of mass."
        t1 = time_ns()
    end
    # move to center of mass
    x_halo = calculate_COM(pos_halo)
    x_disk = calculate_COM(pos_halo)

    if par["verbose"]
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

    resx = par["interp_resolution"]
    resy = par["interp_resolution"]
    resz = par["interp_resolution"]

    if par["verbose"]
        @info "  TSC Interpolation"
        t1 = time_ns()
    end

    rho_tsc, minx, miny, minz, dx, dy, dz = run_tsc(pos_disk, rho_disk, resx, resy, resz, par["verbose"])

    if par["verbose"]
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end    

    if par["verbose"]
        @info "  Filtering IDs by density criterion"
        t1 = time_ns()
    end

    keep_id = trues(length(rho_halo))
    keep_id = find_ids_from_density_criterion(  keep_id, pos_halo, rho_halo, rho_tsc, 
                                                resx, resy, resz,
                                                dx, dy, dz,
                                                minx, miny, minz )

    if par["verbose"]
        t2 = time_ns()
        @info "  Done! We need to delete $(length(findall(keep_id .== false) )) particles."
        @info "  Took $(output_time(t1,t2)) s."
    end

    return keep_id
end