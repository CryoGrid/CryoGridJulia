module SEDIMENT_T_saltConc
#provides sediment module, with state variables T and saltConc
#a coupled system of partial differential equations for temperature and salt concentration is solved.
include("../Common/matlab.jl")
include("../StrataModules/CryoGridConstants.jl")
include("CryoGridTempSaltFunctionalities.jl")
include("../Common/CryoGridTypes.jl")
using MAT
using Statistics


mutable struct stratum
    #mandatory variables
    CONST::CryoGridTypes.constants        #constants
    PARA::CryoGridTypes.parameter         #external service parameters, all other
    STATVAR::CryoGridTypes.statvar      #energy, water content, etc.
    TEMP::CryoGridTypes.temporary         #derivatives in prognostic timestep and optimal timestep
    CHILD::stratum      #pointer to child stratum
    PREVIOUS   #pointer to previous stratum
    NEXT       #pointer to next stratum
    IA_PREVIOUS #pointer to interaction with previous stata
    IA_NEXT #pointer to interaction with next stratum

    #mandatory functionsinitializeFromFile::Function
    initialize::Function
    initialize_statvar::Function
    get_boundary_condition_u::Function
    get_boundary_condition_l::Function
    get_derivatives_prognostic::Function
    get_timestep::Function
    compute_diagnostic_first_cell::Function
    advance_prognostic::Function
    compute_diagnostic::Function
    child::Function
    initializeFromFile::Function

    #-------------------------------------------------------------------------------
    #mandatory functions for each class
    function stratum() #constructor function
        #default values
        this = new()
        this.CONST = CryoGridTypes.constants([0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0])
        this.PARA = CryoGridTypes.parameter([0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0],[0.0], [0.0],[0.0], [0.0], [0.0])
        this.STATVAR = CryoGridTypes.statvar([0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0])
        this.TEMP = CryoGridTypes.temporary([0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0], [1])

        #initialize stratum
        this.initialize = function(this::stratum, sedimenttype, upperPos::Float64=0.0, lowerPos::Float64=-1000.0, layerThick=vec(2 .*ones(1,500)))
            initialize_CONST(this);
            initialize_PARA(this);

            this.STATVAR.upperPos .= upperPos;
            this.STATVAR.lowerPos .= lowerPos;
            if length(layerThick) > 1.0
                this.STATVAR.layerThick = layerThick;
            else
                this.STATVAR.layerThick = dropdims(layerThick .* ones(Int64(floor(abs(lowerPos - upperPos)/layerThick)) - 1, 1), dims=2)
            end

            this.STATVAR.soilType .= sedimenttype.soilType;
            this.STATVAR.mineral .= sedimenttype.mineral;
            this.STATVAR.organic .= sedimenttype.organic;
            this.STATVAR.saltConc .= sedimenttype.salinity;

            return this
        end

        #initialize stratum from file
        #don't use initialize_STATVAR with this, that is already included
        this.initializeFromFile = function(this::stratum, parameters)
            initialize_CONST(this);
            initialize_PARA(this);

            midpoints = parameters[:,1];
            layers = midpoints + 0.5 .* [-(midpoints[1:end-1]-midpoints[2:end]); -(midpoints[end-1] - midpoints[end])];
            layers = [midpoints[1] + 0.5 .*(midpoints[1]-midpoints[2]); layers];

            if layers[end] >= 0.0
                this.STATVAR.upperPos .= -layers[1];
                this.STATVAR.lowerPos .= -layers[end];
            else
                this.STATVAR.upperPos .= layers[1];
                this.STATVAR.lowerPos .= layers[end];
            end

            this.STATVAR.layerThick = abs.(layers[2:end] - layers[1:end-1]);

            this.STATVAR.soilType = parameters[:,6];
            this.STATVAR.mineral = parameters[:,2];
            this.STATVAR.organic = parameters[:,4];
            this.STATVAR.saltConc = parameters[:,5];
            this.STATVAR.porosity = parameters[:,3];
            this.STATVAR.water = 1.0 .- this.STATVAR.organic .- this.STATVAR.mineral;

            this.PARA.alpha   = [4.06] .* parameters[:,10] .+ [6.5e-1] .* -(parameters[:,10] .- 1);
            this.PARA.n       = [2.03] .* parameters[:,10] .+ [1.67] .* -(parameters[:,10] .- 1);

            # make sure Tmelt is in Kelvin
            if mean(parameters[:,9]) <= 0 #then the values are in degree C
                this.STATVAR.Tmelt = parameters[:,9] .+ 273.15;
            else
                this.STATVAR.Tmelt = parameters[:,9]
            end

            this.PARA.a = parameters[:,7];
            this.PARA.b = parameters[:,8];

            #Calculate initial condition T0
            T = CryoGridTempSaltFunctionalities.steadyState(this); #steady state
            #or spinup - put that here
            this.STATVAR.T = T;  #[degree C]

            #conductivity, heat capacity and liquid water content
            saltDiff, thermCond, c_eff = CryoGridTempSaltFunctionalities.getThermalProps_wSalt(this);

            #update struct
            this.STATVAR.saltDiff = saltDiff;
            this.STATVAR.thermCond = thermCond;
            this.STATVAR.c_eff = c_eff;

            return this
        end

        #initialize state variable, such as steady state
        this.initialize_statvar = function(this::stratum)
            initialize_STATVAR(this);

            return this
        end

        #get boundary condition if this is the upper element of the stratigrahy
        this.get_boundary_condition_u = function(this::stratum, forcing) #functions specific for individual class, allow changing from Dirichlet to SEB
            #dirichlet condition with forcing temperature
            T_ub = forcing.TEMP.TForcing;
            thermCond = this.STATVAR.thermCond;
            layerThick = this.STATVAR.layerThick;
            T = this.STATVAR.T;

            #get uppermostGridCell
            #is always 1 if not read from file (version 2017)
            this.TEMP.uppermostGridCell = forcing.TEMP.uppermostGridCell;

            upperIndex = min(this.TEMP.uppermostGridCell[1], length(T)[1]); #we need to useC4D28500 this function before T was initialized as a vector.

            #calculate resulting flux and save temporary variables for conductivity and spatial derivative calculations
            this.TEMP.T_ub = T_ub; #for conductivity and initial steady state
            this.TEMP.heatFlux_ub = thermCond[upperIndex]*(T[upperIndex] .- T_ub) / abs(layerThick[upperIndex] ./ 2.0);#for spatial derivative

            if forcing.TEMP.saltConcForcing[1] == 0.0
                saltFlux_ub = 0.0;
            else
                saltFlux_ub = this.STATVAR.saltDiff[upperIndex] .* (this.STATVAR.saltConc[upperIndex] .- forcing.TEMP.saltConcForcing[1]) ./ abs(this.STATVAR.layerThick[upperIndex] ./2.0);
            end

            this.TEMP.saltFlux_ub = [saltFlux_ub];
            this.TEMP.saltConc_ub = forcing.TEMP.saltConcForcing;

            return this
        end

        #boundary conditions for elements in the middle of the stratigrphy are
        #calculated in the interactions

        #get boundary condition if this is the lowest element of the stratigraphy
        this.get_boundary_condition_l = function(this::stratum, forcing)
            #geothermal heat flux prevails
            this.TEMP.heatFlux_lb = forcing.PARA.heatFlux_lb;
            this.TEMP.saltFlux_lb = [0.0];
            return this
        end

        #calculate spatial derivative
        this.get_derivatives_prognostic = function(this::stratum)
            divT, divsaltConc = get_derivative_temperature_salt(this); #this gives a vector
            this.TEMP.divT = divT;
            this.TEMP.divsaltConc = divsaltConc;
            return this
        end

        #get minimal timestep for this stratum element
        this.get_timestep = function(this::stratum)
            courant_number_temperature = minimum(0.5 * this.STATVAR.c_eff./this.STATVAR.thermCond[1:end-1] .* (this.STATVAR.layerThick).^2);
            courant_number_salt = minimum(0.5 * 1 ./ this.STATVAR.saltDiff[1:end-1] .* (this.STATVAR.layerThick).^2);
            timestep_min_salt = this.PARA.dsaltConc_max[1] ./ maximum(this.TEMP.divsaltConc);
            timestep = min(courant_number_salt, courant_number_temperature, timestep_min_salt)/(3600.0*24.0);
            #println(timestep)
            #timestep = max(timestep, 1) #no time steps smaller than 1 day
            #convert estimate from seconds to days

            return timestep
        end

        #calculate advance in time of state variables
        this.advance_prognostic = function(this::stratum, timestep::Float64) #real timestep derived as minimum of several classes in days
            timestep = timestep*3600.0*24.0; #convert timestep from days to seconds

            this.STATVAR.T = this.STATVAR.T + timestep .* this.TEMP.divT;
            this.STATVAR.saltConc = this.STATVAR.saltConc + timestep .* this.TEMP.divsaltConc;
            return this
        end

        #calculate diagnostic variable that only happen if this stratum element is the
        #upper element in the stratigraphy
        this.compute_diagnostic_first_cell = function(this::stratum, forcing)
            #marine sedimentation could be happening here

            return this
        end

        #calculate diagnostic variables
        this.compute_diagnostic = function(this::stratum, forcing)
            #update conductivity, heat capacity and liquid water content
            this = CryoGridTempSaltFunctionalities.getThermalProps_noSaltDiffusion(this);

            #isostatic movement could happen here
            return this
        end

        return this
    end
end

#-------------------------------------------------------------------------------
#individual functions of this class

function initialize_CONST(this::stratum)
    this = CryoGridConstants.getGlobalConst(this);
    return this
end

function initialize_PARA(this::stratum)
    this = CryoGridConstants.getGlobalPara(this)
    return this
end

function initialize_STATVAR(this::stratum)
    midptDepth = this.STATVAR.upperPos .+ this.STATVAR.layerThick[1] ./ 2.0 .-  cumsum(this.STATVAR.layerThick, dims=1);

    #calculate porosity
    porosityZero = 1.0 .- (this.STATVAR.mineral + this.STATVAR.organic);
    bulkDensityZero = 1.0 ./ ((porosityZero .+ 0.6845) ./ 1.8);
    bulkDensity = bulkDensityZero .+ 0.0037 .* abs.(midptDepth) .^ 0.766; #abs to avoid imaginary numbers!

    porosity = 1.80 .* bulkDensity .^ (-1.0) .- 0.6845;
    porosity[porosity .< 0.03] .= 0.03;

    #update state variables to be vectors
    this.STATVAR.porosity = porosity;
    this.STATVAR.mineral = this.STATVAR.mineral .* (1.0 .- porosity) ./ (this.STATVAR.mineral .+ this.STATVAR.organic);
    this.STATVAR.organic = this.STATVAR.organic .* (1.0 .- porosity) ./ (this.STATVAR.mineral .+ this.STATVAR.organic);
    this.STATVAR.water = 1.0 .- this.STATVAR.organic .- this.STATVAR.mineral;
    this.STATVAR.saltConc = this.STATVAR.saltConc .* ones(size(porosity));

    if this.STATVAR.soilType[1] == 0 #silt
        this.PARA.alpha   = [6.5e-1] .* ones(size(porosity));
        this.PARA.n       = [1.67] .* ones(size(porosity));
    else #sand
        this.PARA.alpha   = [4.06] .* ones(size(porosity));
        this.PARA.n       = [2.03] .* ones(size(porosity));
    end

    #calculate state variables that depend on the more constant ones
    Tmelt, a, b = CryoGridTempSaltFunctionalities.calculateTmelt(this);
    this.STATVAR.Tmelt = Tmelt; #in Kelvin!
    this.PARA.a = a;
    this.PARA.b = b;

    #Calculate initial condition T0
    T = CryoGridTempSaltFunctionalities.steadyState(this); #steady state
    #or spinup - put that here
    this.STATVAR.T = T;  #[degree C]

    #conductivity, heat capacity and liquid water content
    saltDiff, thermCond, c_eff = CryoGridTempSaltFunctionalities.getThermalProps_wSalt(this);

    #update struct
    this.STATVAR.saltDiff = saltDiff;
    this.STATVAR.thermCond = thermCond;
    this.STATVAR.c_eff = c_eff;

    return this
end

function get_derivative_temperature_salt(this)

    #read current state
    T = this.STATVAR.T;
    saltConc = this.STATVAR.saltConc;

    #read boundary conditions
    saltFlux_ub = this.TEMP.saltFlux_ub;
    T_flux_ub = this.TEMP.heatFlux_ub;

    saltFlux_lb = this.TEMP.saltFlux_lb;
    T_flux_lb = this.TEMP.heatFlux_lb;

    #read constants and grid;
    thermCond = this.STATVAR.thermCond;
    L_f = this.CONST.L_f;
    saltDiff = this.STATVAR.saltDiff;

   	layerThick = this.STATVAR.layerThick;
    midptThick = (layerThick[1:end-1] .+ layerThick[2:end])/2;

    #heat flux divergence
    div_F_T = zeros(size(T));

    #above current uppermost Gridcell
    uppermostGridCell = this.TEMP.uppermostGridCell[1];
    @inbounds for i=1:uppermostGridCell-1
        div_F_T[i]=thermCond[i]*(this.TEMP.T_ub[1] .- T[i]) ./ midptThick[i] .- T_flux_ub[1] ./ layerThick[i];
    end

    #upper boundary condition and assuming that k constant near the boundary
    div_F_T[uppermostGridCell]= (thermCond[uppermostGridCell+1] .* (T[uppermostGridCell+1] .- T[uppermostGridCell]) ./ midptThick[uppermostGridCell] .- T_flux_ub[1]) ./ layerThick[uppermostGridCell];

    #use FD for 2nd derivation in space
    div_F_T[uppermostGridCell+1:end-1]=  (thermCond[uppermostGridCell+2:end-1] .* (T[uppermostGridCell+2:end] .- T[uppermostGridCell+1:end-1]) ./ midptThick[uppermostGridCell+1:end] .- thermCond[uppermostGridCell+1:end-2] .* (T[uppermostGridCell+1:end-1] .- T[uppermostGridCell:end-2]) ./ midptThick[uppermostGridCell:end-1]) ./ layerThick[uppermostGridCell+1:end-1];

    #lower BC (dT_dt=geothermal heat flux) and assuming that k constant near the boundary
    div_F_T[end]= (T_flux_lb[1] .- thermCond[end-1] .* (T[end] .- T[end-1]) ./ midptThick[end]) ./layerThick[end];


    #ion flux divergence
    div_F_saltConc = zeros(size(div_F_T));

    @inbounds for i=1:uppermostGridCell-1
        div_F_saltConc[i] = (saltDiff[i+1] .* (saltConc[i+1] .- saltConc[i]) ./ midptThick[i] .- saltFlux_ub[1]) ./layerThick[i];
    end

    #upper boundary condition and assuming that k constant near the boundary
    div_F_saltConc[uppermostGridCell]= (saltDiff[uppermostGridCell+1] .* (saltConc[uppermostGridCell+1] .- saltConc[uppermostGridCell]) ./ midptThick[uppermostGridCell] .- saltFlux_ub[1]) ./layerThick[uppermostGridCell];

    #use FD for 2nd derivation in space
    div_F_saltConc[uppermostGridCell+1:end-1] = (saltDiff[uppermostGridCell+2:end-1] .* (saltConc[uppermostGridCell+2:end] .- saltConc[uppermostGridCell+1:end-1]) ./midptThick[uppermostGridCell+1:end] .- saltDiff[uppermostGridCell+1:end-2] .* (saltConc[uppermostGridCell+1:end-1] .- saltConc[uppermostGridCell:end-2]) ./ midptThick[uppermostGridCell:end-1]) ./layerThick[uppermostGridCell+1:end-1];

    #lower BC zero flux and assuming that D_C constant near the boundary
    div_F_saltConc[end] = (saltFlux_lb[1] .- saltDiff[end-1] .* (saltConc[end] .- saltConc[end-1]) ./ midptThick[end]) ./layerThick[end];

    #Derivative of water content
    dliqWater_dT, dliqWater_saltConc = CryoGridTempSaltFunctionalities.getDerivativeWaterContent(this);
    liqWater = this.STATVAR.liqWater;

    #put everything together

    A = this.STATVAR.c_eff;
    B = L_f .* dliqWater_saltConc;
    D = div_F_T;
    E = (liqWater .+ saltConc .* dliqWater_saltConc);
    F = saltConc .* dliqWater_dT;
    G = div_F_saltConc;

    divT = (.- B .* G .+ D .* E) ./ (A .* E .- B .* F);

    divsaltConc = (.- F .* D + A .* G) ./ (A .* E .- B .* F);

    return divT, divsaltConc
end

end #module
