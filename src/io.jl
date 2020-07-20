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
    pos_info  = Info_Line("POS",  dtype, 3, [1, 1, 1, 1, 1, 1])
    rho_info  = Info_Line("RHO",  dtype, 1, [1, 1, 1, 1, 1, 1])
    hsml_info = Info_Line("HSML", dtype, 1, [1, 1, 1, 1, 1, 1])
    
    # read the blocks
    pos  = read_block_by_name(fi, "POS",  info=pos_info,  parttype=0)
    rho  = read_block_by_name(fi, "RHO",  info=rho_info,  parttype=0) 
    hsml = read_block_by_name(fi, "HSML", info=hsml_info, parttype=0)   

    return pos, rho, hsml 
end

function read_data_gas(filename::String)

    info_3d = Info_Line("POS", dtype, 3, [1, 1, 1, 1, 1, 1])
    info_1d = Info_Line("POS", dtype, 1, [1, 1, 1, 1, 1, 1])

end


function write_to_file(pos_halo, vel_halo, rho_halo, u_halo, m_halo, par)

    ref_header = head_to_obj(par["input_snap"])

end