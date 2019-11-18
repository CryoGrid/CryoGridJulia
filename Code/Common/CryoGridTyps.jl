module CryoGridTyps
export stratigraphy, temporary, forcing, statvar, grid, out, para
struct stratigraphy
    Water::Array{Float64,2}
    Mineral::Array{Float64,2}
    Organic::Array{Float64,2}
    WaterIce::Array{Float64,2}
end
struct temporary
    lat_flux::Array{Float64,2}
    SnowDepth::Array{Float64,2}
    kp::Array{Float64,2}
    cp::Array{Float64,2}
    kn::Array{Float64,1}
    ks::Array{Float64,1}
end
struct forcing
    t_span::Array{Float64,1}
    Tair::Array{Float64,1}
    snowfall::Array{Float64,1}
end
struct statvar
    T::Array{Float64,2}
    H::Array{Float64,2}
end
struct grid
    Zp::Array{Float64,1}
    Zn::Array{Float64,1}
    Zs::Array{Float64,1}
    dxo::Array{Float64,2}
    dxp::Array{Float64,1}
    dxn::Array{Float64,1}
    dxs::Array{Float64,1}
    An::Array{Float64,2}
    As::Array{Float64,2}
    Ao::Array{Float64,2}
    Vp::Array{Float64,2}
end
struct out
    #Temparatur and water content on output time step
    T::Array{Float64,3}
    Water::Array{Float64,3}
    WaterIce::Array{Float64,3}

    #average and min max soil temperature profiles
    T_av::Array{Float64,3}
    T_min::Array{Float64,3}
    T_max::Array{Float64,3}

    #average and min max liquid water content profiles
    W_av::Array{Float64,3}
    W_min::Array{Float64,3}
    W_max::Array{Float64,3}

    #lateral heat fluxes of tiles
    Q_lat::Array{Float64,3}

    #freezing and thawing degree days
    FDD::Array{Float64,3}
    TDD::Array{Float64,3}
    FrostDays::Array{Int64,3}
end
struct para
    SnowDensity::Array{Float64,1}
    WaterDensity::Array{Float64,1}
    WaterDepth::Array{Float64,2}
    Qgeo::Array{Float64,1}
    SnowCoverMax::Array{Float64,2}
end


end
