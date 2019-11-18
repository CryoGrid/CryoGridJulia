module cryoGridTechnical
    using MAT
    #using Dates

    function LayerIndex(vec)
        #clear all
        #vec=[0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

        A=collect(1:length(vec));
        C=A[vec];

        if isempty(C)
            ind_ub=[];
            ind_lb=[];
        else
            ind_ub=C[1];
            ind_lb=C[end];
        end
        return ind_lb, ind_ub
    end

    function makeGrids(PARA)
        GRID=Dict()
        GRID["general"]=Dict()
        GRID["air"]=Dict()
        GRID["snow"]=Dict()
        GRID["soil"]=Dict()
        GRID["lake"]=Dict()
        GRID["lake"]["water"]=Dict()
        GRID["lake"]["ice"]=Dict()

        GRID["snow"]["snowCellSize"]=(PARA["technical"]["SWEperCell"]/(PARA["snow"]["rho_snow"]/PARA["constants"]["rho_w"]))*100.0/100.0;
        a = -(PARA["technical"]["maxSWE"]./(PARA["technical"]["SWEperCell"])+2.0).*GRID["snow"]["snowCellSize"];
        a = trunc(a,10)
        GRID["snow"]["snowGrid"]=collect(a:GRID["snow"]["snowCellSize"]:-GRID["snow"]["snowCellSize"]);
        #tsvd added
        GRID["lake"]["water"]["waterGrid"]=collect(0.0:0.02:PARA["water"]["depth"]-0.02);  #with water (it is recommended to use a grid resolution of about 2cm)
        #GRID["soil"]["soilGrid"]=PARA["technical"]["subsurfaceGrid"]; #no water
        if !isempty(GRID["lake"]["water"]["waterGrid"])
            GRID["soil"]["soilGrid"]=PARA["technical"]["subsurfaceGrid"] + (2.0*GRID["lake"]["water"]["waterGrid"][end]-GRID["lake"]["water"]["waterGrid"][end-1]);
        else
            GRID["soil"]["soilGrid"]=PARA["technical"]["subsurfaceGrid"];
        end

        #tsvd  #K_grid =[GRID["snow"]["snowGrid"]; GRID["soil"]["soilGrid"]]; #grid on which the conductivty information lives (edges of grid cells)
        K_grid =[GRID["snow"]["snowGrid"]; GRID["lake"]["water"]["waterGrid"]; GRID["soil"]["soilGrid"]]; #grid on which the conductivty information lives (edges of grid cells)
        cT_grid=(K_grid[1:end-1]+K_grid[2:end])/2; #grid on which heat capacity and temperature information lives (midpoints of grid cells]
        cT_delta=-cT_grid[1:end-1]+cT_grid[2:end];
        K_delta=-K_grid[1:end-1]+K_grid[2:end]

        GRID["general"]["cT_grid"]=cT_grid;
        GRID["general"]["K_grid"]=K_grid;
        GRID["general"]["cT_delta"]=cT_delta;
        GRID["general"]["K_delta"]=K_delta;

        #set air grid
        GRID["air"]["cT_domain"]= convert(Array{Bool},Base.falses(size(GRID["general"]["cT_grid"])));
        GRID["air"]["cT_domain"][1:length(GRID["snow"]["snowGrid"])]=true;
        GRID["air"]["K_domain"]= convert(Array{Bool},Base.falses(size(GRID["general"]["K_grid"])));
        GRID["air"]["K_domain"][1:length(GRID["snow"]["snowGrid"])]=true;
        GRID["air"]["cT_domain_lb"], GRID["air"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["air"]["cT_domain"]);
        GRID["air"]["K_domain_lb"], GRID["air"]["K_domain_ub"]   = cryoGridTechnical.LayerIndex(GRID["air"]["K_domain"]);

        #set soil grid
        GRID["soil"]["cT_domain"]= convert(Array{Bool},Base.falses(size(GRID["general"]["cT_grid"])));
        #GRID["soil"]["cT_domain"][GRID["general"]["cT_grid"].>0]=true;
        #tsvd    GRID["soil"]["cT_domain"](GRID["general"]["cT_grid"]>0)=1;
        GRID["soil"]["cT_domain"][end-length(GRID["soil"]["soilGrid"])+2:end]=true; #required if there is a lake on top
        GRID["soil"]["cT_domain_lb"], GRID["soil"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["soil"]["cT_domain"]);
        GRID["soil"]["K_domain"]= convert(Array{Bool},Base.falses(size(GRID["general"]["K_grid"])));
        GRID["soil"]["K_domain"][GRID["soil"]["cT_domain_ub"]:end]=true;
        GRID["soil"]["K_domain_lb"], GRID["soil"]["K_domain_ub"]   = cryoGridTechnical.LayerIndex(GRID["soil"]["K_domain"]);

        #set snow grid
        GRID["snow"]["cT_domain"]=convert(Array{Bool},Base.falses(size(GRID["general"]["cT_grid"])));
        GRID["snow"]["K_domain"]=convert(Array{Bool},Base.falses(size(GRID["general"]["K_grid"])));
        GRID["snow"]["cT_domain_lb"], GRID["snow"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["snow"]["cT_domain"]);
        GRID["snow"]["K_domain_lb"], GRID["snow"]["K_domain_ub"]   = cryoGridTechnical.LayerIndex(GRID["snow"]["K_domain"]);

        # JAN: currently these are "secondary" domains which belong to the
        # "primary" soil domain. Lake cells are treated like soil cells during
        # integration of the heat conduction, but the surface properties are
        # different and mixing occurs in summer. Initially, the lake domain is
        # empty, but it is updated in the first time step.

        #set water and ice cover grid
        #GRID["lake"]["cT_domain"] = convert(Array{Bool},Base.falses(size(GRID["general"]["cT_grid"])));
        #GRID["lake"]["K_domain"]= convert(Array{Bool},Base.falses(size(GRID["general"]["K_grid"])));
        #GRID["lake"]["cT_domain_lb"], GRID["lake"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["lake"]["cT_domain"]);
        #GRID["lake"]["K_domain_lb"], GRID["lake"]["K_domain_ub"]   = cryoGridTechnical.LayerIndex(GRID["lake"]["K_domain"]);

        #GRID["lake"]["water"]["cT_domain"]= convert(Array{Bool},Base.falses(size(GRID["general"]["cT_grid"])));
        #GRID["lake"]["water"]["K_domain"]= convert(Array{Bool},Base.falses(size(GRID["general"]["K_grid"])));
        #GRID["lake"]["water"]["cT_domain_lb"], GRID["lake"]["water"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["lake"]["water"]["cT_domain"]);
        #GRID["lake"]["water"]["K_domain_lb"], GRID["lake"]["water"]["K_domain_ub"]   = cryoGridTechnical.LayerIndex(GRID["lake"]["water"]["K_domain"]);

        #GRID["lake"]["ice"]["cT_domain"]=convert(Array{Bool},Base.falses(size(GRID["general"]["cT_grid"])));
        #GRID["lake"]["ice"]["K_domain"]=convert(Array{Bool},Base.falses(size(GRID["general"]["K_grid"])));
        #GRID["lake"]["ice"]["cT_domain_lb"], GRID["lake"]["ice"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["lake"]["ice"]["cT_domain"]);
        #GRID["lake"]["ice"]["K_domain_lb"], GRID["lake"]["ice"]["K_domain_ub"]   = cryoGridTechnical.LayerIndex(GRID["lake"]["ice"]["K_domain"]);
        #GRID["lake"]["unfrozenWaterSurface"]=false;

        #tsvd
        #set water and ice cover grid
        GRID["lake"]["water"]["cT_domain"] = convert(Array{Bool},Base.falses(size(GRID["general"]["cT_grid"])));
        GRID["lake"]["water"]["cT_domain"][GRID["air"]["cT_domain_lb"]+1:GRID["soil"]["cT_domain_ub"]-1]=true;
        GRID["lake"]["water"]["K_domain"] = convert(Array{Bool},Base.falses(size(GRID["general"]["K_grid"])));
        GRID["lake"]["water"]["K_domain"][GRID["air"]["K_domain_lb"]+1:GRID["soil"]["K_domain_ub"]-1]=true;
        GRID["lake"]["water"]["cT_domain_lb"], GRID["lake"]["water"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["lake"]["water"]["cT_domain"]);
        GRID["lake"]["water"]["K_domain_lb"], GRID["lake"]["water"]["K_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["lake"]["water"]["K_domain"]);
        GRID["lake"]["ice"]["cT_domain"] = convert(Array{Bool},Base.falses(size(GRID["general"]["cT_grid"])));
        GRID["lake"]["ice"]["K_domain"] = convert(Array{Bool},Base.falses(size(GRID["general"]["K_grid"])));
        GRID["lake"]["ice"]["cT_domain_lb"], GRID["lake"]["ice"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["lake"]["ice"]["cT_domain"]);
        GRID["lake"]["ice"]["K_domain_lb"], GRID["lake"]["ice"]["K_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["lake"]["ice"]["K_domain"]);
        GRID["lake"]["ice"]["z_ice"] = 0.0;
        GRID["lake"]["ice"]["dz_dt_ice"] = 0.0;
        GRID["lake"]["ice"]["dE_dt_cond_residual"] = 0.0;
        GRID["lake"]["ice"]["melt_flag"] = false;
        return GRID
    end

    function generateTemporary(T::Array{Float64,1}, PARA)
        TEMPORARY=Dict();

        t=PARA["technical"]["starttime"];
        TEMPORARY["outputTime"]=t+PARA["technical"]["outputTimestep"];
        if !isempty(PARA["technical"]["saveInterval"])
            TEMPORARY["saveTime"]=Dates.datetime2julian(Dates.DateTime(Dates.year(Dates.julian2datetime(PARA["technical"]["starttime"]))+1,  parse(Int8,PARA["technical"]["saveDate"][4:5]), parse(Int8,PARA["technical"]["saveDate"][1:2])))-PARA["technical"]["outputTimestep"];
        else
            TEMPORARY["saveTime"]=PARA["technical"]["endtime"]+PARA["technical"]["outputTimestep"];
        end
        TEMPORARY["syncTime"]=round((t + PARA["technical"]["syncTimeStep"])./PARA["technical"]["syncTimeStep"]).*PARA["technical"]["syncTimeStep"];
        TEMPORARY["Qh_sum"]=0.0;
        TEMPORARY["Qe_sum"]=0.0;
        TEMPORARY["Qnet_sum"]=0.0;
        TEMPORARY["Qg_sum"]=0.0;

        #for EB checks
        TEMPORARY["Qsurf_sum"] = 0.0;
        TEMPORARY["dE_dt_SEB_sum"] = zeros(size(T));
        TEMPORARY["dE_dt_cond_sum"] = zeros(size(T));

        TEMPORARY["dE_soil_sens"] = 0.0;
        TEMPORARY["dE_soil_lat"] = 0.0;
        TEMPORARY["dE_soil"] = 0.0;
        TEMPORARY["dE_snow_sens"] = 0.0;
        TEMPORARY["dE_snow_lat"] = 0.0;
        TEMPORARY["dE_snow"] = 0.0;

        TEMPORARY["timestep_sum"]=0.0;

        TEMPORARY["T_sum"]=zeros(size(T));

        TEMPORARY["t_last"]=t;
        TEMPORARY["counter"]=0;
        return t, TEMPORARY
    end

    function generateOUT(GRID, PARA, TEMPORARY)
        OUT=Dict();
        OUT["snow"]=Dict();
        OUT["SEB"]=Dict();
        OUT["soil"]=Dict();
        OUT["snow"]=Dict();
        OUT["WB"]=Dict();
        OUT["EB"]=Dict();
        OUT["location"]=Dict();
        OUT["lateral"]=Dict();

        s1_interval=TEMPORARY["saveTime"]-Dates.datetime2julian(Dates.julian2datetime(TEMPORARY["saveTime"])-Dates.Year(1))
        s2_interval=TEMPORARY["saveTime"]-PARA["technical"]["starttime"]
        s3_interval=PARA["technical"]["endtime"]-Dates.datetime2julian(Dates.julian2datetime(TEMPORARY["saveTime"])-Dates.Year(1))

        s_interval=s1_interval
        if s2_interval>s1_interval
            if s2_interval/s1_interval>=2.0
                s_interval=s1_interval
            else
                s_interval=s2_interval
            end
        end
        if s3_interval<s1_interval
            s_interval=s3_interval
        end

        OUTsize=floor(Int64,s_interval/PARA["technical"]["outputTimestep"])

        OUT["snow"]["outSnow_i"]=[];
        OUT["snow"]["outSnow_a"]=[];
        OUT["snow"]["outSnow_w"]=[];

        OUT["cryoGrid3"]=ones(size(GRID["general"]["cT_grid"],1),OUTsize).*NaN;
        OUT["water"]=[];
        OUT["liquidWater"]=[]; # for distinction between total and liquid water content
        OUT["timestamp"]=ones(1,OUTsize).*NaN;
        OUT["TIMESTEP"]=ones(1,OUTsize).*NaN;

        OUT["SEB"]["Lsta"]=ones(1,OUTsize).*NaN;
        OUT["SEB"]["QE"]=ones(1,OUTsize).*NaN;
        OUT["SEB"]["QH"]=ones(1,OUTsize).*NaN;
        OUT["SEB"]["QNET"]=ones(1,OUTsize).*NaN;
        OUT["SEB"]["QG"]=ones(1,OUTsize).*NaN;
        OUT["SEB"]["Tsurf"]=ones(1,OUTsize).*NaN;
        OUT["SEB"]["albedo_stored"]=ones(1,OUTsize).*NaN;
        OUT["SEB"]["Qsurf"]=ones(1,OUTsize).*NaN;

        OUT["soil"]["topPosition"]=ones(1,OUTsize).*NaN;
        OUT["soil"]["lakeFloor"]=ones(1,OUTsize).*NaN;
        OUT["soil"]["soil"]=ones(1,OUTsize).*NaN;

        OUT["snow"]["outSnow_i"]=[];
        OUT["snow"]["outSnow_a"]=[];
        OUT["snow"]["outSnow_w"]=[];
        OUT["snow"]["topPosition"]=ones(1,OUTsize).*NaN;
        OUT["snow"]["botPosition"]=ones(1,OUTsize).*NaN;


        # derived characteristics and related to geometry
        OUT["location"]["area"] = [];
        OUT["location"]["altitude"] = [];
        OUT["location"]["soil_altitude"] = [];
        OUT["location"]["surface_altitude"] = [];
        OUT["location"]["active_layer_depth_altitude"] = [];
        OUT["location"]["water_table_altitude"] = [];

        # lateral fluxes
        OUT["lateral"]["terrain_index_snow"] = [];
        OUT["lateral"]["water_fluxes"] = [];     # vector containing water fluxes in [m/s] to the current worker
        OUT["lateral"]["snow_fluxes"] = [];      # vector containing snow fluxes in [m SWE / s] to the current worker
        OUT["lateral"]["heat_fluxes"] = [];      # vector containing depth-integrated heat fluxes in [J/m^2 s ] to the current worker


        # water balance (WB)
        # all flows are defined as positive when they go into the soil/snow column
        # cumulative values per output interval in [mm]
        # storage
        OUT["WB"]["dW_soil"] = [];
        OUT["WB"]["dW_snow"] = [];
        # precipitation
        OUT["WB"]["dp_rain"]=[];
        OUT["WB"]["dp_snow"]=[]; # SWE
        # evapotranspiration and sublimation
        OUT["WB"]["de"]=[];
        OUT["WB"]["ds"]=[];
        # runoff
        OUT["WB"]["dr_surface"]=[];
        OUT["WB"]["dr_subsurface"]=[];
        OUT["WB"]["dr_snowmelt"]=[];
        OUT["WB"]["dr_excessSnow"]=[];
        OUT["WB"]["dr_rain"]=[];  # this is only rain on frozen ground

        # energy balance (EB)
        # accumulated energy fluxes per output time in [ J / m^2 ]
        OUT["EB"]["Qg"] = [];         # ground heat flux [positive into ground)
        OUT["EB"]["Qe"] = [];         # latent heat flux [positive into ground)
        OUT["EB"]["Qh"] = [];         # sensible heat flux [positive into ground)
        OUT["EB"]["Qnet"] = [];
        OUT["EB"]["Qgeo"] = [];       # geothermal heat flux

        OUT["EB"]["dE_soil_sens"] = [];
        OUT["EB"]["dE_soil_lat"] = [];
        OUT["EB"]["dE_soil"] = [];
        OUT["EB"]["dE_snow_sens"] = [];
        OUT["EB"]["dE_snow_lat"] = [];
        OUT["EB"]["dE_snow"] = [];
        return OUT
    end

    function interpolateForcingData(t::Float64, FORCING)
        t_span=FORCING["data"]["t_span"]::Array{Float64,2};
        snowfall=FORCING["data"]["snowfall"]::Array{Float64,2};
        rainfall=FORCING["data"]["rainfall"]::Array{Float64,2};
        Lin=FORCING["data"]["Lin"]::Array{Float64,2};
        Sin=FORCING["data"]["Sin"]::Array{Float64,2};
        Tair=FORCING["data"]["Tair"]::Array{Float64,2};
        wind=FORCING["data"]["wind"]::Array{Float64,2};
        q=FORCING["data"]["q"]::Array{Float64,2};
        p=FORCING["data"]["p"]::Array{Float64,2};

        dt=(t_span[2]-t_span[1])::Float64;
        posit=floor(Int32,(t-t_span[1])./dt)+1::Int64;

        #.*size(FORCING["data"]["t_span"]]) + 1;
        # if t<FORCING["data"]["t_span"][posit]
        #     only for security, there can be problems due to numeric resolution
        #     posit=posit-1;
        # end

        FORCING["i"]["snowfall"]=snowfall[posit]+(snowfall[posit+1]-snowfall[posit])*(t-t_span[posit])./dt;
        FORCING["i"]["rainfall"]=rainfall[posit]+(rainfall[posit+1]-rainfall[posit]).*(t-t_span[posit])./dt;
        FORCING["i"]["Lin"]=Lin[posit]+(Lin[posit+1]-Lin[posit])*(t-t_span[posit])./dt;
        FORCING["i"]["Sin"]=Sin[posit]+(Sin[posit+1]-Sin[posit]).*(t-t_span[posit])./dt;
        FORCING["i"]["Tair"]=Tair[posit]+(Tair[posit+1]-Tair[posit]).*(t-t_span[posit])./dt;
        FORCING["i"]["wind"]=wind[posit]+(wind[posit+1]-wind[posit]).*(t-t_span[posit])./dt;
        FORCING["i"]["q"]=q[posit]+(q[posit+1]-q[posit]).*(t-t_span[posit])./dt;
        FORCING["i"]["p"]=p[posit]+(p[posit+1]-p[posit]).*(t-t_span[posit])./dt;
        return FORCING
    end

    function updateBALANCE(T::Array{Float64,1}, wc::Array{Float64,1}, c_cTgrid::Array{Float64,1}, lwc_cTgrid::Array{Float64,1}, BALANCE, GRID, PARA)
        c_w=PARA["constants"]["c_w"]::Float64;
        c_i=PARA["constants"]["c_i"]::Float64;
        c_m=PARA["constants"]["c_m"]::Float64;
        c_o=PARA["constants"]["c_o"]::Float64;
        rho_w=PARA["constants"]["rho_w"]::Float64;
        L_sl=PARA["constants"]["L_sl"]::Float64;

        soil_cT_domain=GRID["soil"]["cT_domain"]::Array{Bool,1};
        soil_cT_mineral=GRID["soil"]["cT_mineral"]::Array{Float64,1};
        soil_cT_organic=GRID["soil"]["cT_organic"]::Array{Float64,1};
        general_K_delta=GRID["general"]["K_delta"]::Array{Float64,1};
        snow_cT_domain=GRID["snow"]["cT_domain"]::Array{Bool,1};
        snow_Snow_i=GRID["snow"]["Snow_i"]::Array{Float64,1};
        snow_Snow_w=GRID["snow"]["Snow_w"]::Array{Float64,1};
        snow_SWEinitial=GRID["snow"]["SWEinitial"]::Float64;

        # at this point the thermal and hydrological state of the soil and snow is calculated
        # energy content soil domain in [J/m^2]
        # distinguished by sensible and latent part
        E_soil_sens_old = BALANCE["energy"]["E_soil_sens"]::Float64;
        E_soil_lat_old = BALANCE["energy"]["E_soil_lat"]::Float64;
        E_soil_old = BALANCE["energy"]["E_soil"]::Float64;

        dE_soil_lat=BALANCE["energy"]["dE_soil_lat"]::Float64;
        dE_soil_sens=BALANCE["energy"]["dE_soil_sens"]::Float64;
        dE_soil=BALANCE["energy"]["dE_soil"]::Float64;

        # energy content snow domain in [J/m^2]
        E_snow_sens_old = BALANCE["energy"]["E_snow_sens"]::Float64;
        E_snow_lat_old = BALANCE["energy"]["E_snow_lat"]::Float64;
        E_snow_old = BALANCE["energy"]["E_snow"]::Float64;
        dE_snow_sens=BALANCE["energy"]["dE_snow_sens"]::Float64;
        dE_snow_lat=BALANCE["energy"]["dE_snow_lat"]::Float64;
        dE_snow=BALANCE["energy"]["dE_snow"]::Float64;

        W_soil_old = BALANCE["water"]["W_soil"]::Float64;
        dW_soil = BALANCE["water"]["dW_soil"]::Float64;

        W_snow_old = BALANCE["water"]["W_snow"]::Float64;
        dW_snow = BALANCE["water"]["dW_snow"]::Float64;

        E_soil_sens = sum(  ( c_w .* lwc_cTgrid[soil_cT_domain] + c_i .* (wc-lwc_cTgrid[soil_cT_domain]) + c_m .* soil_cT_mineral + c_o .* soil_cT_organic ) .* T[soil_cT_domain].* general_K_delta[soil_cT_domain]  );
        BALANCE["energy"]["E_soil_sens"]=E_soil_sens;
        E_soil_lat= sum( rho_w .* L_sl .* lwc_cTgrid[soil_cT_domain] .* general_K_delta[soil_cT_domain]  );
        BALANCE["energy"]["E_soil_lat"] = E_soil_lat;
        E_soil = E_soil_sens + E_soil_lat;
        BALANCE["energy"]["E_soil"] = E_soil
        BALANCE["energy"]["dE_soil_sens"] = dE_soil_sens + E_soil_sens - E_soil_sens_old;
        BALANCE["energy"]["dE_soil_lat"] = dE_soil_lat + E_soil_lat - E_soil_lat_old;
        BALANCE["energy"]["dE_soil"] = dE_soil + E_soil - E_soil_old;

        E_snow_sens = sum( c_cTgrid[snow_cT_domain] .* T[snow_cT_domain] .* general_K_delta[snow_cT_domain] );
        BALANCE["energy"]["E_snow_sens"]=E_snow_sens

        E_snow_lat = sum( rho_w .* L_sl .* lwc_cTgrid[snow_cT_domain].* general_K_delta[snow_cT_domain] );
        BALANCE["energy"]["E_snow_lat"] = E_snow_lat

        E_snow = E_snow_sens + E_snow_lat;
        BALANCE["energy"]["E_snow"] = E_snow
        BALANCE["energy"]["dE_snow_sens"] = dE_snow_sens + E_snow_sens - E_snow_sens_old;
        BALANCE["energy"]["dE_snow_lat"] = dE_snow_lat + E_snow_lat - E_snow_lat_old;
        BALANCE["energy"]["dE_snow"] = dE_snow + E_snow - E_snow_old;

        # water content soil domain in [m]
        W_soil = sum( wc .* general_K_delta[soil_cT_domain] );
        BALANCE["water"]["W_soil"] = W_soil
        BALANCE["water"]["dW_soil"] = dW_soil + (W_soil - W_soil_old)*1000.0; # in [mm]

        # water content snow domain in [m]
        W_snow = sum( snow_Snow_i + snow_Snow_w ) + snow_SWEinitial;
        BALANCE["water"]["W_snow"] = W_snow
        BALANCE["water"]["dW_snow"] = dW_snow + (W_snow - W_snow_old)*1000.0; # in [mm]
        return BALANCE
    end

    function loadSoilTypes(PARA)
        # specify one soil type per row: residualWC [#], fieldCapacity [#], alpha [1/m], n
        PARA["soil"]["soilTypes"] = [ [ 0.00 PARA["soil"]["fieldCapacity"] 4.00 2.0 ];	# sand
                                      [ 0.05 PARA["soil"]["fieldCapacity"] 0.65 1.7 ];	# silt
                                      [ 0.00 0.00                          4.00 2.0 ] ];# pond (freeze curve of sand but fieldCapacity = 0)

        return PARA
    end

    function getAltitude( PARA, GRID )
            #tsvd account for lake surface (snow can accumulate on lake ice)
            if !isempty( GRID["lake"]["water"]["cT_domain_ub"] )# lake water is present
                    altitude = PARA["location"]["initial_altitude"]::Float64 - GRID["general"]["K_grid"][GRID["lake"]["water"]["cT_domain_ub"]::Int64]::Float64;
            elseif !isempty( GRID["lake"]["ice"]["cT_domain_ub"] ) # lake ice is present
                    altitude = PARA["location"]["initial_altitude"]::Float64 - GRID["general"]["K_grid"][GRID["lake"]["ice"]["cT_domain_ub"]::Int64]::Float64;

            else
                    altitude = PARA["location"]["initial_altitude"]::Float64 - GRID["general"]["K_grid"][GRID["soil"]["cT_domain_ub"]::Int64]::Float64; # no lake
            end
            return altitude::Float64
    end

    function getSurfaceAltitude(PARA, GRID)
            surface_altitude = PARA["location"]["initial_altitude"]::Float64 - GRID["general"]["K_grid"][GRID["air"]["cT_domain_lb"]::Int64+1]::Float64;

            #if ~isempty(GRID["snow"]["cT_domain_ub"])
            #    surface_altitude = PARA["location"]["altitude"] + sum( GRID["general"]["K_delta"](GRID["snow"]["cT_domain"]) );
            #else
            #    surface_altitude = PARA["location"]["altitude"];
            #end
            return surface_altitude::Float64
    end

    function getActiveLayerDepthAltitude( PARA, GRID, T::Array{Float64,1})
            # Function that gives the cT index of the bottom cell of the bucket.
            ald_altitude = NaN;
            if isempty(GRID["snow"]["cT_domain_ub"]) && T[GRID["soil"]["cT_domain_ub"]]>0.0; # Condition to work on infiltration
                    T_i=T[GRID["soil"]["cT_domain"]]::Array{Float64,1};
                    i=1;
                    i_max=200;
                    @inbounds @fastmath while  T_i[i]>0.0 && i<=i_max
                            i=i+1;
                    end
                    ald_cT_index=i-1;
                    ald_altitude = PARA["location"]["initial_altitude"]::Float64-GRID["general"]["K_grid"][GRID["soil"]["cT_domain_ub"]::Int64 + ald_cT_index]::Float64;
            end
            return ald_altitude::Float64
    end

    function getWaterTableAltitude(T::Array{Float64,1}, wc::Array{Float64,1}, GRID, PARA)
            T_i=T[GRID["soil"]["cT_domain"]]::Array{Float64,1};
            K_delta=GRID["general"]["K_delta"][GRID["soil"]["cT_domain"]]::Array{Float64,1};  #in m
            if !isempty(GRID["snow"]["cT_domain_ub"]) || T_i[1]<=0.0
                    waterTable=NaN;
            else
                    waterTable = PARA["location"]["altitude"]::Float64;
                    water = wc;
                    porosity=(1.0 - GRID["soil"]["cT_mineral"]::Array{Float64,1} - GRID["soil"]["cT_organic"]::Array{Float64,1});
                    i=1;
                    @inbounds @fastmath while water[i]<porosity[i] && T[i]>0.0
                            waterTable = waterTable - K_delta[i];
                            i=i+1;
                    end
                    if T[i]<=0.0
                            waterTable=NaN;
                            #disp('Dry to permafrost - waterTable=NaN \n')# Loop was stopped by the temperature condition so no water table
                    end
            end
            return waterTable::Float64
    end

    function getSoilAltitude(PARA, GRID)
            #tsvd GRID.lake.cT_domain replaced by GRID["lake"]["water"]["cT_domain"]
            if !isempty( GRID["lake"]["water"]["cT_domain_ub"] )
                    soil_altitude = PARA["location"]["initial_altitude"] - GRID["general"]["K_grid"][GRID["lake"]["water"]["cT_domain_lb"]+1];
            else
                    soil_altitude = PARA["location"]["initial_altitude"] - GRID["general"]["K_grid"][GRID["soil"]["cT_domain_ub"]];
            end
            return soil_altitude
    end



end
