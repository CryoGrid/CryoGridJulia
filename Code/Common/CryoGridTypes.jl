module CryoGridTypes
export constants, parameter, statvar, temporary, forcingparameter, forcingtemporary, forcingdata, sediment

#stratigraphy, forcing, statvar, grid, out, para

struct constants #24 Array{Float64,1}
    day_sec::Array{Float64,1}
    year_sec::Array{Float64,1}
    rho_w::Array{Float64,1}
    rho_i::Array{Float64,1}
    L_f::Array{Float64,1}
    L_v::Array{Float64,1}
    c_i::Array{Float64,1}
    c_w::Array{Float64,1}
    c_o::Array{Float64,1}
    c_m::Array{Float64,1}
    c_a::Array{Float64,1}
    k_a::Array{Float64,1}
    k_i::Array{Float64,1}
    k_w::Array{Float64,1}
    k_o::Array{Float64,1}
    k_m::Array{Float64,1}
    sigma::Array{Float64,1}
    R::Array{Float64,1}
    R_a::Array{Float64,1}
    kappa::Array{Float64,1}
    Tmelt::Array{Float64,1}
    g::Array{Float64,1}
    benthicSalt::Array{Float64,1}
    tau::Array{Float64,1}
end

mutable struct parameter # 11 Array{Float64,1}
    heatFlux_lb::Array{Float64,1}
    albedo::Array{Float64,1}
    epsilon::Array{Float64,1}
    rs::Array{Float64,1}
    z::Array{Float64,1}
    z0::Array{Float64,1}

    max_timestep::Array{Float64,1}
    dE_max::Array{Float64,1}
    dt_max::Array{Float64,1}

    a::Array{Float64,1}
    b::Array{Float64,1}
end

mutable struct statvar # 19 Array{Float64,1}
    deltaT::Array{Float64,1}

    upperPos::Array{Float64,1}
    lowerPos::Array{Float64,1}

    layerThick::Array{Float64,1}

    energy::Array{Float64,1}
    T::Array{Float64,1}
    saltConc::Array{Float64,1}

    waterIce::Array{Float64,1}
    liqWater::Array{Float64,1}
    water::Array{Float64,1}
    ice::Array{Float64,1}
    mineral::Array{Float64,1}
    organic::Array{Float64,1}
    porosity::Array{Float64,1}
    air::Array{Float64,1}

    thermCond::Array{Float64,1}
    c_eff::Array{Float64,1}
    Tmelt::Array{Float64,1}
    soilType::Array{Float64,1}
end

mutable struct temporary # 4 Array{Float64,1}
    T_ub::Array{Float64,1}
    heatFlux_ub::Array{Float64,1}
    heatFlux_lb::Array{Float64,1}
    divT::Array{Float64,1}
end

struct forcingparameter # 3 Array{int64,1}, 15 Array{Float64,1}
    SL_no::Array{Int64,1}
    TF_no::Array{Int64,1}
    saltForcingSwitch::Array{Int64,1}

    T_freeze::Array{Float64,1}
    T_IceSheet::Array{Float64,1}
    IS::Array{Float64,1}
    benthicSalt::Array{Float64,1}
    latitude::Array{Float64,1}
    longitude::Array{Float64,1}
    altitude::Array{Float64,1}
    heatFlux_lb::Array{Float64,1}
    startForcing::Array{Float64,1}
    dtForcing::Array{Float64,1}
    endForcing::Array{Float64,1}

    domain_depth::Array{Float64,1}

    start_time::Array{Float64,1}
    end_time::Array{Float64,1}
end

mutable struct forcingtemporary # 3 Array{Float64,1}
    TForcing::Array{Float64,1}
    surfaceState::Array{Float64,1}
    saltConcForcing::Array{Float64,1}
end

mutable struct forcingdata # 6 Array{float64,1}
    timeForcing::Array{Float64,1}
    seaLevel::Array{Float64,1}
    airTemp::Array{Float64,1}
    TForcing::Array{Float64,1}
    surfaceState::Array{Float64,1}
    saltConcForcing::Array{Float64,1}
end

struct sediment # 4 Array{int64,1}
    soilType::Array{Float64,1}
    mineral::Array{Float64,1}
    organic::Array{Float64,1}
    salinity::Array{Float64,1}
end

end