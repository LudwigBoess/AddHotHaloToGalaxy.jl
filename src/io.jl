using GadgetIO
import YAML

# helper function to get elapsed time in seconds
function output_time(t1, t2)
    return Float16((t2-t1)*1.e-9)
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

    par["Mdm"] *= par["UnitMass_in_g"]

    par["UnitEnergy_in_cgs"] = par["UnitMass_in_g"] * par["UnitVelocity_in_cm_s"]^2
    par["SpecEnergyUnit"]    = par["UnitEnergy_in_cgs"] / par["UnitMass_in_g"]
    par["DensityUnit"]       = par["UnitMass_in_g"] / par["UnitLength_in_cm"]^3

    return par
end

function read_parameter_file(parameter_file::String)
    # read parameters into dictionary
    par = YAML.load(open(parameter_file))
    # add calculated parameters and return
    return add_constants(par)
end

function read_old_pos(fi)
    
end