using GadgetIO
import YAML

# helper function to get elapsed time in seconds
function output_time(t1, t2)
    return Float64((t2-t1))*1.e-9
end

function calculate_COM!(x)
    com = zeros(3)
    
    @inbounds for i = 1:3
        com[i] = sum(x[i,:]) / size(x,2)
    end
    @inbounds for N = 1:size(x,2), i = 1:3
        x[i,N] -= com[i]
    end
    return x
end

function find_COM(x)
    com = zeros(3)
    
    @inbounds for i = 1:3
        com[i] = sum(x[i,:]) / size(x,2)
    end

    return com
end

function add_constants(par)

    par["G"]  = 6.672e-8
    par["mp"] = 1.6726e-24
    par["mu"] = 4.0/( 8.0 - 5.0*(1.0-par["xH"]))
    # yhelium = ( 1.0 - par["xH"] ) / ( 4.0 * par["xH"] )
    # par["mu"] = (1.0 + 4.0 * yhelium) / (1.0 + 3.0 * yhelium + 1.0)
    par["kB"] = 1.3806e-16
    
    par["r200"] = par["v200"] * par["UnitLength_in_cm"]
    par["rs"]   = par["r200"] / par["cc"]
    par["rc"]   = par["rc_frac"] * par["rs"]
    par["a"]    = par["rs"] * sqrt( 2.0 * ( log(1.0 + par["cc"] ) - par["cc"] / (1.0 + par["cc"]) ) )

    par["UnitEnergy_in_cgs"]  = par["UnitMass_in_g"] * par["UnitVelocity_in_cm_s"]^2
    par["SpecEnergyUnit"]     = par["UnitEnergy_in_cgs"] / par["UnitMass_in_g"]
    par["UnitDensity_in_cgs"] = par["UnitMass_in_g"] / par["UnitLength_in_cm"]^3

    return par
end

function calculate_pars_from_galaxy(par, galaxy)

    # mass of the dm halo
    par["Mdm"] = sum(galaxy["PartType1"]["MASS"]) * par["UnitMass_in_g"]

    # if number of gas particles is set to zero we use the same number as in the gas disk.
    if iszero(par["npartgas"])
        par["npartgas"] = size(galaxy["PartType0"]["MASS"],1)
    end

    # Mgas is set to 0 calculate it from m_gas + m_star + m_bulge
    if iszero(par["Mgas"])
        parttypes = ["PartType0", "PartType2", "PartType3"]
        
        Mgas = 0.0
        for part ∈ parttypes
            if haskey(galaxy, part)
                Mgas += sum(galaxy[part]["MASS"])
            end
        end
        # save in cgs units
        par["Mgas"] = Mgas

    end

    # calculate number of gas particles from total gas mass divided by disk gas particle mass
    # this makes sure that SPH particles have the same mass
    if par["npartgas"] == -1
        par["npartgas"] = floor(Int64, par["Mgas"] / galaxy["PartType0"]["MASS"][1,1] )
    end

    par["Mgas"] *= par["UnitMass_in_g"]

    return par
end

function read_parameter_file(parameter_file::String)
    # read parameters into dictionary
    par = YAML.load(open(parameter_file))
    # add calculated parameters and return
    return add_constants(par)
end

# function read_pos_rho_hsml_galaxy(fi::String, ic_format::String="float")

#     if ic_format == "float"
#         dtype = Float32
#     elseif ic_format == "double"
#         dtype = Float64
#     end

#     # define info for readin
#     pos_info  = InfoLine("POS",  dtype, 3, [1, 1, 1, 1, 1, 1])
#     rho_info  = InfoLine("RHO",  dtype, 1, [1, 1, 1, 1, 1, 1])
#     hsml_info = InfoLine("HSML", dtype, 1, [1, 1, 1, 1, 1, 1])
    
#     # read the blocks
#     pos  = read_block(fi, "POS",  info=pos_info,  parttype=0)
#     rho  = read_block(fi, "RHO",  info=rho_info,  parttype=0) 
#     hsml = read_block(fi, "HSML", info=hsml_info, parttype=0)   

#     return pos, rho, hsml 
# end

# function read_data_gas(fi::String, ic_format::String="float")

#     if ic_format == "float"
#         dtype = Float32
#     elseif ic_format == "double"
#         dtype = Float64
#     end

#     info_3d = InfoLine("POS", dtype, 3, [1, 1, 1, 1, 1, 1])
#     info_1d = InfoLine("POS", dtype, 1, [1, 1, 1, 1, 1, 1])

#     pos  = read_block(fi, "POS",  info=info_3d, parttype=0)
#     vel  = read_block(fi, "VEL",  info=info_3d, parttype=0) 
#     u    = read_block(fi, "U",    info=info_1d, parttype=0)
#     m    = read_block(fi, "MASS", info=info_1d, parttype=0)

#     calculate_COM!(pos)

#     return Float32.(pos), Float32.(vel), Float32.(u), Float32.(m)
# end

# function read_data_collisionless(fi::String, parttype::Int64, ic_format::String="float")

#     if ic_format == "float"
#         dtype = Float32
#     elseif ic_format == "double"
#         dtype = Float64
#     end

#     info_3d = InfoLine("POS", dtype, 3, [1, 1, 1, 1, 1, 1])
#     info_1d = InfoLine("POS", dtype, 1, [1, 1, 1, 1, 1, 1])

#     pos  = read_block(fi, "POS",  info=info_3d, parttype=parttype)
#     vel  = read_block(fi, "VEL",  info=info_3d, parttype=parttype)
#     m    = read_block(fi, "MASS", info=info_1d, parttype=parttype) 

#     calculate_COM!(pos)

#     return Float32.(pos), Float32.(vel), Float32.(m)
# end

function read_galaxy_data(filename::String, ic_format::String="float")

    if ic_format == "float"
        dtype = Float32
    elseif ic_format == "double"
        dtype = Float64
    end
    
    data = Dict()

    h = read_header(filename)

    # read dummy particles
    info_3d = InfoLine("POS", dtype, 3, [1, 1, 1, 1, 1, 1])
    info_1d = InfoLine("POS", dtype, 1, [1, 1, 1, 1, 1, 1])

    # read positios of dm halo
    dm_pos  = read_block(filename, "POS", info=info_3d, parttype=1)

    # calculate center of mass of DM halo
    com = find_COM(dm_pos)

    parttype_arr = ["PartType0", "PartType1", "PartType2", "PartType3", "PartType4", "PartType5"]

    for p = 1:size(parttype_arr,1)

        if h.npart[p] > 0
            # read data
            data[parttype_arr[p]] = Dict()
            data[parttype_arr[p]]["POS"] = read_block(filename, "POS",  info=info_3d, parttype=(p-1))
            # move to center of mass
            for i = 1:h.npart[p], dim = 1:3
                data[parttype_arr[p]]["POS"][dim,i] -= com[dim]
            end

            data[parttype_arr[p]]["VEL"]  = read_block(filename, "VEL",  info=info_3d, parttype=(p-1))
            data[parttype_arr[p]]["MASS"] = read_block(filename, "MASS", info=info_1d, parttype=(p-1))
        
        end
    end

    data["PartType0"]["U"]   = read_block(filename, "U",  info=info_1d, parttype=0)
    data["PartType0"]["RHO"] = read_block(filename, "RHO",  info=info_1d, parttype=0)

    return data

end

function write_to_file(galaxy, pos_halo, vel_halo, u_halo, rho_halo, m_halo, par)

    if par["verbose"]
        @info "  Merging snapshot and halo"
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

     # update number of gas particles in snapshot
    header.npart[1] += Nhalo
    header.nall[1]  += Nhalo
    header.massarr  .= 0.0

    if par["black_hole"]
        header.npart[end] = 1
        header.nall[end]  = 1
    end

    # store total number of particles
    Ntotal = sum(header.npart)
    
    # set up id array
    ids = UInt32.(collect(1:Ntotal))

    # allocate empty arrays for particles
    pos = Matrix{Float32}(undef, 3, Ntotal )
    vel = Matrix{Float32}(undef, 3, Ntotal )
    u   = zeros(Float32, header.npart[1] )
    rho = zeros(Float32, header.npart[1] )
    m   = zeros(Float32, Ntotal )

    # assign halo particles
    pos[:,1:Nhalo] = Float32.( pos_halo )
    vel[:,1:Nhalo] = Float32.( vel_halo )
    u[1:Nhalo]     = Float32.( u_halo )
    rho[1:Nhalo]   = Float32.( rho_halo )
    m[1:Nhalo]     = Float32.( m_halo )
 
    # assign old gas particles
    Nstart = Nhalo + 1
    Npart  = header.npart[1] - Nhalo
    pos[:,Nstart:Nstart+Npart-1] = galaxy["PartType0"]["POS"]
    vel[:,Nstart:Nstart+Npart-1] = galaxy["PartType0"]["VEL"]
    u[Nstart:Nstart+Npart-1]     = galaxy["PartType0"]["U"]
    rho[Nstart:Nstart+Npart-1]   = galaxy["PartType0"]["RHO"]
    m[Nstart:Nstart+Npart-1]     = galaxy["PartType0"]["MASS"]
    Nstart += Npart
    
    # assign collisionless particles
    @inbounds for parttype = 1:5
        Npart = header.npart[parttype+1]
        if Npart > 0
            pos[:,Nstart:Nstart+Npart-1] = galaxy["PartType$parttype"]["POS"]
            vel[:,Nstart:Nstart+Npart-1] = galaxy["PartType$parttype"]["VEL"]
            m[Nstart:Nstart+Npart-1]     = galaxy["PartType$parttype"]["MASS"]
            Nstart += Npart
        end
    end

    if par["black_hole"]

        pos[:,end] = Float32[0.0, 0.0, 0.0]
        vel[:,end] = Float32[0.0, 0.0, 0.0]
        m[end]     = Float32(par["M_bh"] * 1.989e+33 / par["UnitMass_in_g"]) # Convert to Gadget units

    end

    B = Float32.(par["B0"] .* ones( 3, header.npart[1]) )

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
    write_block(f, m,   "MASS")
    write_block(f, u,   "U")
    write_block(f, B,   "BFLD")

    if par["black_hole"]
        Mbh = par["M_bh"] * 1.989e+33 / par["UnitMass_in_g"] # Convert to Gadget units
        write_block(f, Float32(Mbh), "BHMA")
        write_block(f, Float32(0.0), "BHMD")
        write_block(f, Int32(1), "BHPC")
    end

    close(f)

    if par["verbose"]
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

end