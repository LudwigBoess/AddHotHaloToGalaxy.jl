using Roots
using GadgetIO
using Base.Threads
using ProgressMeter

# which r does fullfill the following equation?
# Mgas(<r) / Mgas_total = q
# q(r) = 4*!PI*rc^3*rho_0 / Mgas * ( r/rc - atan( r/rc ) )
# therefore the root of the following equation is found by newton approach
@inline function rootfinding_for_r(r, rc::Float64, rho0::Float64, Mgas::Float64, q::Float64)
    @fastmath 4π * rc^3 * rho0 / Mgas * ( r/rc - atan( r/rc ) ) - q
end

function find_rmax(rc::Float64, rc_frac::Float64, rho0::Float64, Mgas::Float64)

    q = 1.0
    start = rc/rc_frac

    # define helper function because root finding allows only one input parameter
    r_helper(r) = rootfinding_for_r(r, rc, rho0, Mgas, q)
    r_max = find_zero(r_helper, (start, 100start)) # use bisection, since it works

    return r_max
end

function construct_large_glass(pos::Array{Float64,2}, n::Int32, ntot::Int64, 
                               nx::Int64, ny::Int64, nz::Int64)

    pos_out = zeros(ntot, 3)

    count = 0

    @inbounds @fastmath for ix = 0:(nx-1), iy = 0:(ny-1), iz = 0:(nz-1)
        
        @inbounds @fastmath @simd for i = 1:n
            pos_out[count*n+i, 1] = pos[i,1] + ix 
            pos_out[count*n+i, 2] = pos[i,2] + iy  
            pos_out[count*n+i, 3] = pos[i,3] + iz
        end

        count += 1
    end

    # @inbounds @fastmath @simd for dim = 1:3
    #     pos_out[:,dim] ./= nn[dim]
    # end

    pos_out[:,1] ./= nx
    pos_out[:,2] ./= ny
    pos_out[:,3] ./= nz

    x = 2.0 .* pos_out[:,1] .- 1.0
    y = 2.0 .* pos_out[:,2] .- 1.0
    z = 2.0 .* pos_out[:,3] .- 1.0

    return x, y, z
end

@inline function xy_to_phi(x::Vector{Float64}, y::Vector{Float64})

    phi = zeros(size(x,1))
    
    @threads for i = 1:size(x,1)
        if ( x[i] > 0.0 )
            phi[i] = atan( y[i]/ x[i] )
        elseif ( x[i] == 0.0 )
            phi[i] = y[i] / abs(y[i]) * 0.5π
        elseif ( x[i] < 0.0 && y[i] >= 0.0 )
            phi[i] = atan( y[i]/ x[i])  + π
        elseif ( x[i] < 0.0 && y[i] < 0.0)
            phi[i] = atan( y[i]/ x[i])  - π
        end
    end

    return phi

end

function construct_spherical_coordinates(x::Array{Float64,1}, y::Array{Float64,1}, z::Array{Float64,1},
                                         nparticles::Integer)

    n_samples = size(x,1)

    phi       = zeros(n_samples)
    theta     = zeros(n_samples)
    r         = zeros(n_samples)

    for i = 1:n_samples
        r[i]     = sqrt( x[i]^2 + y[i]^2 + z[i]^2 )
        #phi[i]   = atan(y[i], x[i])
        theta[i] = atan( sqrt(x[i]^2 + y[i]^2 )/ z[i] )
    end
    phi   = xy_to_phi(x, y)

    r = r.^3

    id = findall( r .< 1.0)
    count_in = size(id,1)

    # cut out sphere so that particle number equals about nparticles
    id = findall( r .<= nparticles/count_in )

    # rescale sphere
    r .*= count_in / nparticles

    return r[id], phi[id], theta[id]

end

function glass_spherical_dist(nparticles::Int64, glass_file::String)

    h = head_to_obj(glass_file)

    nnn = nparticles / h.npart[1] * ( 8.0 / ( 4.0/3.0 * π ) )

    # number of stacks in x, y, z direction
    nx = ny = nz = floor(Int64, nnn^(1.0/3.0)) + 1

    pos_info = InfoLine("POS", Float32, 3, [ 1, 0, 0, 0, 0, 0 ] )
    pos = read_block(glass_file, "POS", info=pos_info, parttype=0) ./ h.boxsize
    pos = Float64.(pos) # convert from Float32 to Float64

    n =  h.npart[1]
    ntot = n * nx * ny * nz
    
    x, y, z = construct_large_glass(pos, n, ntot, nx, ny, nz)

    return x, y, z

    r, theta, phi = construct_spherical_coordinates(x, y, z, nparticles)

    return [r phi theta]

end

@inline function calc_rho_rmax(rho0::Float64, rmax::Float64, rc::Float64, 
                               UnitMass_in_g::Float64, UnitLength_in_cm::Float64)

    rho_rmax  = rho0 / ( 1.0 + (rmax/rc)^2 )
    rho_rmax /= (UnitMass_in_g / UnitLength_in_cm^3)
    return rho_rmax
end

function rescale_r!(r_in::Array{Float64,1}, 
                   rc::Float64, rc_frac::Float64, rho0::Float64, Mgas::Float64)

    start = rc/rc_frac

    @inbounds for i = 1:size(r_in,1)
        q = r_in[i]
        # define helper function because root finding allows only one input parameter
        r_helper(r) = rootfinding_for_r(r, rc, rho0, Mgas, q)
        r           = find_zero(r_helper, start)
        r_in[i]     = r
        start       = r
    end

    nothing
end

function gashalo_beta(par)

    # convert and compute parameters
    Mgas = par["Mgas"] * par["UnitMass_in_g"]

    n    = par["npartgas"]
    
    # find maximum sampling radius of gas halo
    rmax = find_rmax(par["rc"], par["rc_frac"], par["rho0"], Mgas)

    if par["verbose"]
        @info "    Maximum sampling radius: rmax = $(rmax/par["UnitLength_in_cm"])   rmax/rc = $(rmax/par["rc"])"
    end

    # find density at maximum sampling radius
    rho_rmax = calc_rho_rmax(par["rho0"], rmax, par["rc"], par["UnitMass_in_g"], par["UnitLength_in_cm"] )

    if par["verbose"]
        @info "    rho_rmax = $rho_rmax"
    end

    # if no glass file specified sample halo from random distribution
    if par["glass_file"] == ""
        if par["verbose"]
            @info "    Sampling random positions."
        end
        r_phi_theta = rand(n, 3)
        r_phi_theta[:,3] = @. acos( 2.0 * r_phi_theta[:,3] - 1.0)
        r_phi_theta[:,2] = @. 2π * r_phi_theta[:,2]

    # otherwise construct it from stacking a glass file
    else
        if par["verbose"]
            @info "    Sampling from glass file."
        end
        r_phi_theta = glass_spherical_dist(n, par["glass_file"])
        n           = size(r_phi_theta,1)
    end

    if par["verbose"]
        
        @info "    maximum r = $(maximum(r_phi_theta[:,1]))"
        @info "    rescaling r, rmax = $(rmax/par["UnitLength_in_cm"])"

    end

    sort!(r_phi_theta[:,1])

    rescale_r!(r_phi_theta[:,1], par["rc"], par["rc_frac"], par["rho0"], Mgas)


    if par["verbose"]
        @info "    Converting back into code units."
    end
    if par["verbose"]
        @info "    maximum r = $(maximum(r_phi_theta[:,1]))"
    end
    max_sampled = maximum(r_phi_theta[:,1])
    r_phi_theta[:,1] .*= (rmax / max_sampled)

    r_phi_theta[:,1] ./= par["UnitLength_in_cm"]

    if par["verbose"]
        @info "    maximum r = $(maximum(r_phi_theta[:,1]))"
        #error("done")
        @info "    Converting spherical to cartesian coordinates."
    end

    # convert back to cartesian coordinates
    x = r_phi_theta[:,1] .* cos.(r_phi_theta[:,2]) .* sin.(r_phi_theta[:,3])
    y = r_phi_theta[:,1] .* sin.(r_phi_theta[:,2]) .* sin.(r_phi_theta[:,3])
    z = r_phi_theta[:,1]                           .* cos.(r_phi_theta[:,3]) 

    return r_phi_theta[:,1], [ x y z ]
end

function gashalo_vphi(r_in::Vector{Float64}, par)

    r = r_in .* par["UnitLength_in_cm"]
    
    vphi = @. sqrt( par["G"] * par["Mdm"] * r ) / ( r + par["a"] ) * par["lambda"]

    return vphi ./ par["UnitVelocity_in_cm_s"]
end

function convert_vphi_to_cart(vphi::Vector{Float64}, pos::Array{Float64,2})

    phi = xy_to_phi(pos[:,1], pos[:,2])

    vel = zeros(size(pos,1), 3)
    
    for i = 1:size(pos,1)
        vel[i,1] = -vphi[i] * sin(phi[i])
        vel[i,2] =  vphi[i] * cos(phi[i])
    end

    return vel

end

function halo_temp(r, par)

    F0 = par["rc"] / (par["a"]^2 + par["rc"]^2)^2 * 
        ( π/2. * ( par["a"]^2 - par["rc"]^2 ) + par["rc"] * 
        (par["a"]^2 + par["rc"]^2) / (par["a"]+r) 
        - (par["a"]^2 - par["rc"]^2) * atan(r/par["rc"]) 
        - par["rc"] * par["a"] * log( (par["a"]+r)^2/(r^2+par["rc"]^2) ) )

    F1 = π^2/(8*par["rc"]) - (atan(r/par["rc"]))^2/(2*par["rc"]) - atan(r/par["rc"])/r

    T = par["G"] * par["mu"] * par["mp"] / par["kB"] * ( 1.0 + r^2 / par["rc"]^2 ) * ( par["Mdm"] * F0 + 4 * π * par["rc"]^3 * par["rho0"] * F1 )

    return T
end

function halo_u(r, par)

    # compute temperature as function of radius
    T = zeros(size(r,1))
    @threads for i = 1:size(r,1)
        @inbounds T[i] = halo_temp(r[i], par)
    end

    if par["verbose"]
        Tc = halo_temp(par["rc"], par)
        @info "    Characteristic temperature of cluster: Tc = $(Float32(Tc)) K"
    end

    # convert temperature back to internal energy
    return T ./ ( par["UnitVelocity_in_cm_s"]^2 .* 2.0/3.0 .* par["mu"] ./ par["kB"] )
end

function halo_rho(r, par)
    # β-model with β = 2/3
    return @. par["rho0"] / ( 1.0 + (r * par["UnitLength_in_cm"] / par["rc"])^2) / par["UnitDensity_in_cgs"]
end

function construct_hot_halo(par)

    if par["verbose"]
        @info "  Gas halo positions"
        t1 = time_ns()
    end

    r, pos = gashalo_beta(par) 

    if par["verbose"]
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

    if par["verbose"]
        @info "  Computing particle Masses"
        t1 = time_ns()
    end
    
    # find mass of a gas particle
    particle_mass_gas = par["Mgas"] / size(pos,1)
    m = particle_mass_gas .* ones(size(pos,1))

    if par["verbose"]
        t2 = time_ns()
        @info "  Done! Mass of a single particle: Mg = $particle_mass_gas"
        @info "  Took $(output_time(t1,t2)) s"
    end

    if par["verbose"]
        @info "  Computing internal energy of particles"
        t1 = time_ns()
    end

    # halo_u missing!
    u = halo_u(r, par)

    if par["verbose"]
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

    if par["verbose"]
        @info "  Computing velocity of particles"
        t1 = time_ns()
    end

    vphi = gashalo_vphi( sqrt.( pos[:,1].^2 .+ pos[:,2].^2 ), par)
    vel  = convert_vphi_to_cart(vphi, pos)

    if par["verbose"]
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

    if par["verbose"]
        @info "  Computing ideal density"
        t1 = time_ns()
    end

    rho = halo_rho(r, par)

    if par["verbose"]
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

    # compensate for 'little h'.
    pos .*= par["h"]
    u   .*= par["h"]
    
    return pos, vel, u, m, rho

end