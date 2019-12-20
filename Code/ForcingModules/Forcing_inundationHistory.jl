module Forcing_inundationHistory
include("../Common/matlab.jl")
include("../ForcingModules/CryoGridInputDataSets.jl")
include("../Common/CryoGridTypes.jl")

using MAT

mutable struct forcing
    DATA::CryoGridTypes.forcingdata
    TEMP::CryoGridTypes.forcingtemporary
    PARA::CryoGridTypes.forcingparameter

    initialize::Function
    interpolate_forcing::Function
    initalize_from_file::Function
    load_forcing_from_mat::Function

    function forcing()
        this = new()
        this.DATA = CryoGridTypes.forcingdata([0.0], [0.0], [0.0], [0.0], [0.0], [0.0]);
        this.TEMP = CryoGridTypes.forcingtemporary([0.0], [0.0], [0.0]);
        this.PARA = CryoGridTypes.forcingparameter([0], [0], [0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0]);


        this.initialize = function(this::forcing, SL_no::Int64=666, TF_no::Int64=2, lat::Float64=69.945, lon::Float64=-134.0, startForcing::Float64=-2000.0)
           this.PARA.SL_no .= SL_no; #2
           this.PARA.TF_no .= TF_no; #3;
           this.PARA.T_freeze .= -1.6;
           this.PARA.T_IceSheet .= 0.0;
           this.PARA.IS .= 100.0;
           this.PARA.benthicSalt .= 989.63;
           this.PARA.saltForcingSwitch .= 1; # 0 to switch off salt Forcing
           this.PARA.latitude .= lat;
           this.PARA.longitude .= lon;
           this.PARA.altitude .= 14.634;
           this.PARA.heatFlux_lb .= 0.05;
           this.PARA.startForcing .= startForcing;
           this.PARA.dtForcing .= 100.0;
           this.PARA.endForcing .= 0.0;

           this.PARA.domain_depth .= 10000.0; #?

           return this
        end

        this.load_forcing_from_mat = function(this::forcing)
           #addpath(genpath('forcing'))

           #overwrite given altitude and heat flux by values from EaseGrid
           #and Davies et. al. respectively
           this.PARA.altitude .= CryoGridInputDataSets.getElevation(this.PARA.latitude, this.PARA.longitude);
           this.PARA.heatFlux_lb .= CryoGridInputDataSets.getQ_Davies(this.PARA.latitude, this.PARA.longitude);

           if this.PARA.SL_no[1] == 666
               this = generateForcing_testing(this);
           else
               this = generateForcing_fromData(this);
           end

           #println("generated Forcing data")
           #println([this.DATA.timeForcing[1], this.DATA.timeForcing[20], this.DATA.timeForcing[end - 20], this.DATA.timeForcing[end]])
           #println([this.DATA.TForcing[1], this.DATA.TForcing[20], this.DATA.TForcing[end - 20], this.DATA.TForcing[end]])

           #calculate time in days for the main programm
           this.PARA.start_time .= this.PARA.startForcing .* 365.25;
           this.PARA.end_time .= this.PARA.endForcing .*365.25;

           return this
        end

        this.interpolate_forcing = function(t, this::forcing) #t comes in days!
           t = t/365.25; #convert to t in years
           this.TEMP.TForcing = matlab.interp1(this.DATA.timeForcing, this.DATA.TForcing, t, "linear");
           #forcing.TEMP.surfaceState = matlab.interp1(forcing.DATA.timeForcing, forcing.DATA.surfaceState, t, 'nearest');

           this.TEMP.saltConcForcing = matlab.interp1(this.DATA.timeForcing, this.DATA.saltConcForcing, t, "nearest");

           return this
        end

        return this
    end
end


function generateForcing_testing(this::forcing)
    timeForcing = this.PARA.startForcing[1]:this.PARA.dtForcing[1]:this.PARA.endForcing[1];
    this.DATA.timeForcing = timeForcing;
    forcingData = zeros(size(this.DATA.timeForcing));
    saltConcForcing  = zeros(size(this.DATA.timeForcing));
    coverage = zeros(size(this.DATA.timeForcing));
    time_inundation = 0;
    surfaceState = zeros(size(this.DATA.timeForcing));

    testcase = Int64(this.PARA.TF_no[1]);

    if testcase == 1 #one-third glacier, one third aerial, one third flooded

        i_start = 1;
        #does this exist?
        i_end = Int64(floor(length(forcingData)/3));
        forcingData[i_start:i_end] .= this.PARA.T_IceSheet;
        coverage[i_start:i_end] .= this.PARA.IS[1] + 100.0;
        surfaceState[i_start:i_end] .= -1;

        i_start = i_end + 1;
        i_end = 2*i_end;
        forcingData[i_start:i_end] = LinRange(-17.0, -9.0, i_end - i_start+1);
        coverage[i_start:i_end] .= 0.0;
        surfaceState[i_start:i_end] .= 1;

        i_start = i_end + 1;
        i_end = length(forcingData);
        forcingData[i_start:i_end] = LinRange(this.PARA.T_freeze[1], 0, i_end - i_start+1);
        coverage[i_start:i_end] = LinRange(-30.0, -2.0, i_end - i_start+1);
        surfaceState[i_start:i_end] .= 0;
        saltConcForcing[i_start:i_end] .= this.PARA.benthicSalt;

        time_inundation = i_end - i_start + 1;

    elseif testcase == 2 #subaerial ramp
        forcingData = LinRange(-20.0 ,20.0 , length(this.DATA.timeForcing));
        coverage = zeros(size(this.DATA.timeForcing));
        time_inundation = 0;
        surfaceState = ones(size(this.DATA.timeForcing));

    elseif testcase == 3 #zero temperature, saltConc change with ramp

        i_start = Int64(floor(length(forcingData)/2.0));
        forcingData = 0.0 .*ones(size(this.DATA.timeForcing));
        coverage = zeros(size(this.DATA.timeForcing));
        coverage[i_start:end] .= 20.0;
        time_inundation = length(this.DATA.timeForcing) - i_start;
        surfaceState = ones(size(this.DATA.timeForcing));
        surfaceState[i_start:end] .= 0;

        i_start = Int64(floor(length(forcingData)/3));
        i_end = 2*i_start;
        saltConcForcing[i_start:i_end] = LinRange(0,this.PARA.benthicSalt[1], i_end - i_start + 1);
        saltConcForcing[i_end:end] .= this.PARA.benthicSalt;

    elseif testcase == 4 #subsea ramp

        forcingData = LinRange(5.0,-10.0, length(this.DATA.timeForcing));
        coverage = LinRange(2.0,200.0, length(this.DATA.timeForcing));
        time_inundation = length(this.DATA.timeForcing);
        surfaceState = ones(size(this.DATA.timeForcing));

        saltConcForcing[:] .= this.PARA.benthicSalt;

    elseif testcase == 5 #zero temperature, saltConc change without ramp

        i_start = Int64(floor(length(forcingData)/2));
        forcingData = -10.0 .* ones(size(this.DATA.timeForcing));
        coverage = zeros(size(this.DATA.timeForcing));
        coverage[i_start:end] .= 20.0;
        time_inundation = length(this.DATA.timeForcing) - i_start;
        surfaceState = ones(size(this.DATA.timeForcing));
        surfaceState[i_start:end] .= 0.0;

        saltConcForcing[i_start:end] .= this.PARA.benthicSalt;

    elseif testcase == 6 #rapid temperature and saltConc change

        i_start = Int64(floor(length(forcingData)/6.0));

        saltConcForcing[1:i_start] .= 0.0;
        saltConcForcing[i_start:end] .= this.PARA.benthicSalt;

        i_start = 2*i_start;
        forcingData[1:i_start] .= -10.0;
        forcingData[i_start:end] .= 0.0;

    else
        prinln("no such testcase!")
    end

    this.DATA.TForcing = forcingData;

    #this.DATA.coverage = coverage;
    #this.DATA.time_inundation = time_inundation;
    this.DATA.surfaceState = surfaceState;

    if this.PARA.saltForcingSwitch == 0
        this.DATA.saltConcForcing = zeros(size(saltConcForcing));
    else
        this.DATA.saltConcForcing = saltConcForcing;
    end

    return this
end


function generateForcing_fromData(this::forcing)
    #generateForcing accepts (lat, lon, zsb) for the location and
    #returns the timestamp, the duration
    #of inundation, surface temperature and an index defining the surface
    #(submarine, subaerial, subglacial).

    lat = this.PARA.latitude;
    lon = this.PARA.longitude;
    zsb = this.PARA.altitude;

    timeForcing = [ this.PARA.startForcing[1]:this.PARA.dtForcing[1]:this.PARA.endForcing[1];];
    this.DATA.timeForcing = timeForcing;
    saltConcForcing = timeForcing .* 0.0; #"outer" salt concentration is zero for subaerial and subglacial conditions
    coverage = timeForcing .* 0.0;    #positive values are thickness of ice-coverage
                                      #negative values are depth of covering sea
                                      #(over current elevation)
                                      #0 accounts for no coverage

    #define sea level history
    if this.PARA.SL_no[1] == 1  #Lambeck et al.
        file = matopen("ForcingData/subseaInput.mat")
        elevation = read(file, "elevation");
        nominal_elevation = read(file, "nominal_elevation");
        time = read(file, "time");

        seaLevel = matlab.interp1([time[1:325,1].*1000.0; 50000.0], [nominal_elevation[1:325,1]; nominal_elevation[325,1]], - timeForcing, "linear");

    elseif this.PARA.SL_no[1] == 2 #Walbroeck et al.
        file = matopen("ForcingData/sealevel_Waelbroeck.mat")  #starting at -429.5 kyrs BP
        time_W = read(file, "time_W");
        RSL_W = read(file, "RSL_W");

        time_W = [0 ; time_W];
        RSL_W = [0 ; RSL_W]; #add RSL=0 for t=0

        seaLevel = matlab.interp1(time_W .* 1000.0 ,RSL_W,timeForcing, "linear");

    elseif this.PARA.SL_no[1] == 3 #Grant et al.
        file = matopen("ForcingData/sealevel_Grant_prob.mat")  #starting at -492 kyrs BP
        time_G = read(file, "time_G");
        RSL_G = read(file, "RSL_G");

        seaLevel = matlab.interp1(dropdims(rot180(time_G .* 1000.0), dims=2),dropdims(rot180(RSL_G .- RSL_G[1]), dims=2), timeForcing, "linear"); #set RSL=0 at year t=0 (i.e. substract -RSL_G(1))
    else
        #disp("Invalid switch for sea level history")
    end
    this.DATA.seaLevel = seaLevel;


    #define air temperature forcing
    if this.PARA.TF_no[1] == 1 #use GISP-2 data
        file = matopen("ForcingData/GISP2.mat");
        GISP2_T = read(file, "GISP2_T");

        baseT = -5.0;   #air temperature offset for GISP data
        @inbounds for i = size(GISP2_T,1):-1:2
            if GISP2_T[i,1] == GISP2_T[i-1,1]
                GISP2_T[i,:]=[];
            end
        end
        deltaT = matlab.interp1([0.0; GISP2_T[:,1]*1000.0; 50000.0], [GISP2_T[1,2]; GISP2_T[:,2]; GISP2_T[end,2]], - timeForcing, "linear");
        deltaT = deltaT - deltaT[end,1];
        T = baseT + deltaT;
        forcingData = T;

    elseif this.PARA.TF_no[1] == 2 #use CLIMBER-2 data
        file = matopen("ForcingData/CLIMBER2/SATcor_ANN_C2ip_500k.mat")          #CLIMBER-2 SAT, corrected for thermal offset
        SAT_ANNipt2 = read(file, "SAT_ANNipt2");
        lat_C2ip = read(file, "lat_C2ip");
        lon_C2ip = read(file, "lon_C2ip2");
        time_C2ip = read(file, "time_C2ip");

        min_lat_dist, this_ind_lat = findmin(abs.(lat_C2ip .- lat));
        min_lon_dist, this_ind_lon = findmin(abs.(lon_C2ip .- lon));

        #do we have something like squeeze??
        SAT_C2 = SAT_ANNipt2[this_ind_lon[2],this_ind_lat[2],:]; # time series of CLIMBER-2 SAT for chosen grid cell
        SAT_C2_interp = matlab.interp1(vec(time_C2ip), SAT_C2, timeForcing, "linear");  #interpolate linear to given time steps
        #clear SAT_ANNipt2  % free memory
        forcingData = SAT_C2_interp;
    else
        #disp("Invalid switch TF for specifying temperature forcing")

    end
    this.DATA.airTemp = forcingData;

    #set temperature data for inundated sites to sea water temperature
    ind_inundation = seaLevel .> zsb; #index of years of inundation
    time_inundation = this.PARA.dtForcing * sum(ind_inundation); #duration of inundation

    @inbounds for i = 1:length(forcingData) #adapt forcing temperature for inundated sites
        if zsb[1] < seaLevel[i] #site is inundated
            waterDepth = seaLevel[i] - zsb[1];    #depth water column
            if(waterDepth > 30) #below 30m T sea bottom equals T_freeze)
                T_seaWater = forcing.PARA.T_freeze;
            elseif waterDepth <= 30 && waterDepth > 2 #linear scaling between 30m t0 2m water depth
                T_seaWater = 1.0 ./ 14.0 * (this.PARA.T_freeze / 2.0 * waterDepth - this.PARA.T_freeze);
            else  #between 2m and 0m T sea bottom equals 0�C
                T_seaWater = 0.0;
            end
            forcingData[i] = T_seaWater[1];

            coverage[i] = - waterDepth[1];

            saltConcForcing[i] = this.PARA.benthicSalt[1];
        end
    end

    #set temperature data for (warm-based) ice sheet covered sites to bottom ice sheet temperature
    if this.PARA.IS[1] > 0 #use of CLIMBER-2 ice sheet data
        file = matopen("ForcingData/CLIMBER2/Ice_Thickness_C2ipt_500kyrs.mat") #get CLIMBER-2 ice sheet thickness matrix (same grid resolution as SAT matrix
        ITHipt = read(file, "ITHipt");

        if !@isdefined this_ind_lat
            min_lat_dist, this_ind_lat = findmin(abs.(lat_C2ip .- lat));
            min_lon_dist, this_ind_lon = findmin(abs.(lon_C2ip .- lon));

            file = matopen("ForcingData/CLIMBER2/SATcor_ANN_C2ip_500k.mat")          #CLIMBER-2 SAT, corrected for thermal offset
            time_C2ip = read(file, "time_C2ip");
        end

        ITH = ITHipt[this_ind_lon[2],this_ind_lat[2],:]; #time series of CLIMBER-2 ITH for chosen grid cell
        ITH = matlab.interp1(vec(time_C2ip), ITH, timeForcing, "linear");
        #clear ITHipt  % free memory
        ind_IceCover = ITH .> this.PARA.IS;
        forcingData[ind_IceCover] .= this.PARA.T_IceSheet; #set temperature to ice sheet bottom temperatures for sites with a ice sheet coverage larger than S_IS
        ind_IS = ITH .> this.PARA.IS[1]; #index of years of IceSheet coverage thicker than S_IS

        coverage[ind_IceCover] = ITH[ind_IceCover];
    end


    surfaceState = timeForcing .* 0 .+ 1;      #define surface state as 1 everywhen (subaerial==1)
    surfaceState[ind_inundation] .= 0;         #set inundation periods to 0 (submarine==0)

    if @isdefined ind_IS
        surfaceState[ind_IS] .= -1;            #set subglacial periods to -1 (subglacial ==-1)
    end


    this.DATA.TForcing = forcingData;
    #this.DATA.coverage = coverage;
    #this.DATA.time_inundation = time_inundation;
    this.DATA.surfaceState = surfaceState;
    if this.PARA.saltForcingSwitch == 0
        this.DATA.saltConcForcing = zeros(size(saltConcForcing));
    else
        this.DATA.saltConcForcing = saltConcForcing;
    end

    return this
end

end
