using GadgetIO
import YAML

# helper function to get elapsed time in seconds
function output_time(t1, t2)
    return Float64((t2-t1))*1.e-9
end

function add_constants(par)

    par["G"]  = 6.672e-8
    par["mp"] = 1.6e-24
    par["mu"] = 4.0/( 8.0 - 5.0*(1.0-par["xH"]))
    par["kB"] = 1.3806e-16
    
    par["r200"] = par["v200"] * par["UnitLength_in_cm"]
    par["rs"]   = par["r200"] / par["cc"]
    par["rc"]   = par["rc_frac"] * par["rs"]
    par["a"]    = par["rs"] * sqrt( 2.0 * ( log(1.0 + par["cc"] ) - par["cc"] / (1.0 + par["cc"]) ) )

    par["Mdm"]  *= par["UnitMass_in_g"]
    #par["Mgas"] *= par["UnitMass_in_g"]

    par["UnitEnergy_in_cgs"]  = par["UnitMass_in_g"] * par["UnitVelocity_in_cm_s"]^2
    par["SpecEnergyUnit"]     = par["UnitEnergy_in_cgs"] / par["UnitMass_in_g"]
    par["UnitDensity_in_cgs"] = par["UnitMass_in_g"] / par["UnitLength_in_cm"]^3

    return par
end

function read_parameter_file(parameter_file::String)
    # read parameters into dictionary
    par = YAML.load(open(parameter_file))
    # add calculated parameters and return
    return add_constants(par)
end

function read_pos_rho_hsml_galaxy(fi::String, ic_format::String="float")

    if ic_format == "float"
        dtype = Float32
    elseif ic_format == "double"
        dtype = Float64
    end

    # define info for readin
    pos_info  = InfoLine("POS",  dtype, 3, [1, 1, 1, 1, 1, 1])
    rho_info  = InfoLine("RHO",  dtype, 1, [1, 1, 1, 1, 1, 1])
    hsml_info = InfoLine("HSML", dtype, 1, [1, 1, 1, 1, 1, 1])
    
    # read the blocks
    pos  = read_block(fi, "POS",  info=pos_info,  parttype=0)
    rho  = read_block(fi, "RHO",  info=rho_info,  parttype=0) 
    hsml = read_block(fi, "HSML", info=hsml_info, parttype=0)   

    return pos, rho, hsml 
end

function read_data_gas(fi::String, ic_format::String="float")

    if ic_format == "float"
        dtype = Float32
    elseif ic_format == "double"
        dtype = Float64
    end

    info_3d = InfoLine("POS", dtype, 3, [1, 1, 1, 1, 1, 1])
    info_1d = InfoLine("POS", dtype, 1, [1, 1, 1, 1, 1, 1])

    pos  = read_block(fi, "POS",  info=info_3d, parttype=0)
    vel  = read_block(fi, "VEL",  info=info_3d, parttype=0) 
    u    = read_block(fi, "U",    info=info_1d, parttype=0)
    m    = read_block(fi, "MASS", info=info_1d, parttype=0)

    calculate_COM!(pos)

    return Float32.(pos), Float32.(vel), Float32.(u), Float32.(m)
end

function read_data_collisionless(fi::String, parttype::Int64, ic_format::String="float")

    if ic_format == "float"
        dtype = Float32
    elseif ic_format == "double"
        dtype = Float64
    end

    info_3d = InfoLine("POS", dtype, 3, [1, 1, 1, 1, 1, 1])
    info_1d = InfoLine("POS", dtype, 1, [1, 1, 1, 1, 1, 1])

    pos  = read_block(fi, "POS",  info=info_3d, parttype=parttype)
    vel  = read_block(fi, "VEL",  info=info_3d, parttype=parttype)
    m    = read_block(fi, "MASS", info=info_1d, parttype=parttype) 

    calculate_COM!(pos)

    return Float32.(pos), Float32.(vel), Float32.(m)
end

function write_to_file(pos_halo, vel_halo, u_halo, m_halo, par)

    if par["verbose"]
        @info "  Reading an merging snapshot and halo"
        t1 = time_ns()
    end

    # store number of particles in halo 
    Nhalo  = size(u_halo,1)

    if par["verbose"]
        @info "  Nhalo = $Nhalo"
        t1 = time_ns()
    end

    # read the header of the input snap
    header = head_to_obj(par["input_snap"])

    # store total number of particles
    Ntotal = sum(header.npart) + Nhalo
    
    # set up id array
    ids = UInt32.(collect(1:Ntotal))

    # allocate empty arrays for particles
    pos = zeros(Float32, Ntotal, 3 )
    vel = zeros(Float32, Ntotal, 3 )
    u   = zeros(Float32, header.npart[1]+Nhalo )
    m   = zeros(Float32, Ntotal )

    # assign halo particles
    pos[1:Nhalo,:] = Float32.(pos_halo)
    vel[1:Nhalo,:] = Float32.(vel_halo)
    u[1:Nhalo]     = Float32.(u_halo)
    m[1:Nhalo]     = Float32.(m_halo)
 
    # read old gas particles
    Nstart = Nhalo + 1
    Npart  = header.npart[1]
    pos[Nstart:Nstart+Npart-1,:], vel[Nstart:Nstart+Npart-1,:], u[Nstart:Nstart+Npart-1], m[Nstart:Nstart+Npart-1] = read_data_gas(par["input_snap"], par["ic_format"])
    Nstart += Npart
    
    # read collisionless particles
    @inbounds for parttype = 1:5
        Npart = header.npart[parttype+1]
        if Npart > 0
            pos[Nstart:Nstart+Npart-1,:], vel[Nstart:Nstart+Npart-1,:], m[Nstart:Nstart+Npart-1] = read_data_collisionless(par["input_snap"], parttype, par["ic_format"])
            Nstart += Npart
        end
    end

    # update number of gas particles in snapshot
    header.npart[1] += Nhalo
    header.nall[1]  += Nhalo
    header.massarr  .= 0.0

    B = Float32.(par["B0"] .* ones(Ntotal, 3) )

    if par["verbose"]
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

    if par["verbose"]
        @info "  Writing IC file"
        t1 = time_ns()
    end

    # write to file 
    f = open(par["output_file"], "w")
    write_header(f, header)
    write_block(f, pos, "POS")
    write_block(f, vel, "VEL")
    write_block(f, ids, "ID")
    write_block(f, m, "MASS")
    write_block(f, u, "U")
    write_block(f, B, "BFLD")
    close(f)

    if par["verbose"]
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

end