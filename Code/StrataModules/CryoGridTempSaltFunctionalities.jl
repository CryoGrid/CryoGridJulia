module CryoGridTempSaltFunctionalities
    include("../Common/matlab.jl")
    using MAT

    function calculateTmelt(this)
        #calculate Tmelt in K

        file = matopen("StrataModules/lookup_table_fc.mat");
        fc_parameters = read(file, "fc_parameters");

        a_sand = matlab.interp1(fc_parameters[:,1], fc_parameters[:,2], this.STATVAR.saltConc, "linear").* this.STATVAR.porosity.^(-1.0);
        b_sand = matlab.interp1(fc_parameters[:,1], fc_parameters[:,3], this.STATVAR.saltConc, "linear");
        a_silt = matlab.interp1(fc_parameters[:,1], fc_parameters[:,4], this.STATVAR.saltConc, "linear").* this.STATVAR.porosity.^(-1.0);
        b_silt = matlab.interp1(fc_parameters[:,1], fc_parameters[:,5], this.STATVAR.saltConc, "linear");
        a = a_silt .*(this.STATVAR.soilType .== 1) + a_sand .*(this.STATVAR.soilType .== 0);
        b = b_silt .*(this.STATVAR.soilType .== 1) + b_sand .*(this.STATVAR.soilType .== 0);

        Tmelt = this.CONST.Tmelt ./ this.CONST.L_f .* (-this.CONST.R .* this.STATVAR.saltConc.*this.CONST.Tmelt) .+ [273.15];

        return Tmelt, a, b
    end

    function steadyState(this)

        #get thermal conductivities
        k_water = this.CONST.k_w;
        k_ice = this.CONST.k_i;
        k_organic = this.CONST.k_o;
        k_mineral = this.CONST.k_m;

        if ~hasfield(typeof(this.PREVIOUS), :STATVAR)
            TForcing = this.TEMP.T_ub;
        else
            TForcing = this.PREVIOUS.STATVAR.T[length(this.PREVIOUS.STATVAR.T)];
        end

        Q = this.PARA.heatFlux_lb;

        midpoints = this.STATVAR.upperPos .- this.STATVAR.layerThick[1] / 2.0 .- cumsum(this.STATVAR.layerThick);

        Tmelt = this.STATVAR.Tmelt; #in Kelvin
        Tmelt_inDegreeC = Tmelt .- 273.15; #convert to °C to compare against temperatures.

        #allocate memory
        T_0 = zeros(size(midpoints));
        #apply forcing to surface
        T_0[1] = TForcing[1]; #TForcing is a Vector with 1 element

        @inbounds for i = 2:length(T_0)
            T = T_0[i-1];

            a1 = 1.0 ./ this.STATVAR.porosity[i] - this.PARA.a[i] .* Tmelt_inDegreeC[i] - this.PARA.b[i] .* Tmelt_inDegreeC[i].^2.0;

            if T < Tmelt_inDegreeC[i]
                liqWater = (1.0 ./ (a1 + this.PARA.a[i] .* T + this.PARA.b[i] .* T.^2.0));
            else
                liqWater = this.STATVAR.porosity[i];
            end

            thermCond = (liqWater .* sqrt.(k_water) + (this.STATVAR.porosity[i] - liqWater) .* sqrt.(k_ice) + this.STATVAR.mineral[i]*sqrt.(k_mineral) + this.STATVAR.organic[i]*sqrt.(k_organic)).^2.0;

        	T_0[i] = T_0[i-1] .- 1.0 ./ thermCond[1] .* Q[1] .* (midpoints[i] .- midpoints[i-1]);
        end

        return T_0
    end

    function getThermalProps_noSaltDiffusion(this)
        #read constants
        #W / m K
        k_water = this.CONST.k_w;
        k_ice = this.CONST.k_i;
        k_organic = this.CONST.k_o;
        k_mineral = this.CONST.k_m;

        #J / (K m^3)
        c_water = this.CONST.c_w;
        c_organic = this.CONST.c_o;
        c_mineral = this.CONST.c_m;
        c_ice = this.CONST.c_i;

        L_f = this.CONST.L_f;

        #fraction of volume
        mineral = this.STATVAR.mineral;
        porosity = this.STATVAR.porosity;
        organic = this.STATVAR.organic;

        a = this.PARA.a;
        b = this.PARA.b;
        Tmelt = this.STATVAR.Tmelt; #in Kelvin
        Tmelt_inDegreeC = Tmelt .- 273.15; #convert to °C to compare against temperatures.

        #temperature and upper boundary
        T_ub = this.TEMP.T_ub;
        T = this.STATVAR.T;



        #-------------- assign thermal conductivity -------------------------%
        #conductivity lives on the edges
        #interpolate liqWater to the edges
        layerDepth = [this.STATVAR.upperPos; this.STATVAR.upperPos .- cumsum(this.STATVAR.layerThick)];
        thermCond = zeros(size(layerDepth));

        #calculate for upper boundary
        i = 1;
        a1 = 1.0 ./ porosity[i] - a[i]*abs(Tmelt_inDegreeC[i])^b[i];

        if (T_ub[1] .+ T[i]) / 2.0 <= Tmelt_inDegreeC[i]
            liqWater = 1.0 / (a1 .+ a[i]*abs((T_ub[1] .+ T[i]) / 2.0) .^ b[i]);
        else
            liqWater = porosity[i];
        end

        thermCond[i] = (liqWater[1] .* sqrt.(k_water[1]) .+ (porosity[i] .- liqWater[1]) .* sqrt.(k_ice[1]) .+ mineral[i] .* sqrt.(k_mineral[1]) .+ organic[i] .* sqrt.(k_organic[1])) .^2.0;

        #calculate for inner edges
        @inbounds for i = 2:length(layerDepth)-1
            #take the average of liqWater's of grid cells above and below
            a1 = 1.0 / porosity[i-1] - a[i-1]*abs(Tmelt_inDegreeC[i-1])^b[i-1];
            if T[i-1] <= Tmelt_inDegreeC[i-1]
                liqWater = 1.0 / (a1 + a[i-1]*abs(T[i-1])^b[i-1]);
            else
                liqWater = porosity[i-1];
            end

            a1 = 1.0 / porosity[i] - a[i]*abs(Tmelt_inDegreeC[i])^b[i];
            if T[i] <= Tmelt_inDegreeC[i]
                liqWater = (liqWater + 1.0 / (a1 + a[i]*abs(T[i])^b[i])) / 2.0;
            else
                liqWater = (liqWater+porosity[i]) / 2.0;
            end

            thermCond[i] = (liqWater[1] .* sqrt.(k_water[1]) .+ ((porosity[i] .+ porosity[i-1]) ./ 2.0 .- liqWater[1]) .* sqrt.(k_ice[1]) .+  (mineral[i-1] .+ mineral[i]) ./ 2.0 .* sqrt.(k_mineral[1]) .+ (organic[i-1] .+ organic[i]) ./ 2.0 .* sqrt.(k_organic[1])) .^2.0;
        end

        #extrapolate for lower edge - this is needed for calculation of
        #interaction
        i = length(thermCond);
        thermCond[i] = thermCond[i-1];

        #------------------ determine bulk conductivity and capacity------%
        # also save liqWater here
        #Heat Capacity lives on the midpoints of the cells
        c_eff = zeros(size(this.STATVAR.layerThick));
        liqWater = zeros(size(this.STATVAR.layerThick));
        midptDepth = this.STATVAR.upperPos .- this.STATVAR.layerThick[1] ./ 2.0 .- cumsum(this.STATVAR.layerThick);

        @inbounds for i = 1:length(midptDepth)
            a1 = 1.0 / porosity[i] - a[i]* abs(Tmelt_inDegreeC[i])^b[i];

            #determine water content
            if T[i] <= Tmelt_inDegreeC[i]
                liqWater[i] = 1.0 / (a1 + a[i]*abs(T[i])^b[i]);
            else
                liqWater[i] = porosity[i];
            end

            if T[i] <= Tmelt_inDegreeC[i]
                d_liqWater = a[i] * b[i] * abs(T[i])^(b[i] - 1.0) / (a1 + a[i]*abs(T[i])^b[i]) / (a1 + a[i]*abs(T[i])^b[i]);
            else
                d_liqWater = 0.0;
            end

            c_eff[i] = c_mineral[1] * mineral[i] + c_water[1] * liqWater[i] + c_ice[1] * (porosity[i]-liqWater[i]) + c_organic[1] * organic[i] +  L_f[1] * d_liqWater[1];
        end

        #---------------- update this struct ---------------------------%
        this.STATVAR.thermCond = thermCond;
        this.STATVAR.c_eff = c_eff;
        this.STATVAR.water = liqWater; #argh!!!
        this.STATVAR.liqWater = liqWater;

        return this
    end

    function getThermalProps_wSalt(this)

        #read constants
        c_w = this.CONST.c_w;
        c_o = this.CONST.c_o;
        c_m = this.CONST.c_m;
        c_i = this.CONST.c_i;

        L_f = this.CONST.L_f;
        tau = this.CONST.tau;

        k_a = this.CONST.k_a;
        k_w = this.CONST.k_w;
        k_i = this.CONST.k_i;
        k_o = this.CONST.k_o;
        k_m = this.CONST.k_m;

        layerThick = this.STATVAR.layerThick;
        saltDiff0 = this.PARA.saltDiff0;

        #Calculate liqWater and derivatives
        liqWater, Tmelt = watercontent_temperature_salt(this); #Tmelt in Kelvin!
        this.STATVAR.liqWater = liqWater;
        this.STATVAR.Tmelt = Tmelt;

        porosity = this.STATVAR.porosity;
        ice = porosity .- liqWater;
        mineral = this.STATVAR.mineral;
        organic = this.STATVAR.organic;
        air = 1.0 .- liqWater .- ice .- mineral .- organic;


        dliqWater_dT = getDerivativeWaterContent(this);

        c = c_m[1] .* mineral .+ c_o[1] .* organic + c_w[1] .* liqWater .+ c_i[1] .* ice;
        c_eff = c .+ L_f .* dliqWater_dT;

        saltDiff = saltDiff0 .* liqWater ./ tau;
        saltDiff = [saltDiff[1], (saltDiff[1:end-1] .* saltDiff[2:end]) ./ (saltDiff[1:end-1] .+ saltDiff[2:end]), saltDiff[end]];


        thermCond = (liqWater .* k_w[1] .^0.5 .+ ice .* k_i[1] .^0.5 + mineral .* k_m[1] .^0.5 .+ organic .* k_o[1] .^0.5 + air .* k_a[1] .^0.5).^2;
        thermCond = [thermCond[1], (thermCond[1:end-1] .* thermCond[2:end] .* (layerThick[1:end-1] .+ layerThick[2:end])) ./ (thermCond[1:end-1] .* layerThick[2:end] .+ thermCond[2:end] .* layerThick[1:end-1]), thermCond[end]];

        return saltDiff, thermCond, c_eff
    end

    function watercontent_temperature_salt(this, T=NaN, saltConc=NaN)
        #T is required to be in Kelvin!
        #convert to Kelvin, if you use this.STATVAR.T or
        #give T in Kelvin via optionalArgs

        if ~isnan(T) #T and saltConc given as optional args
            #T is in Kelvin
        else
            #Convert from degreeC to Kelvin
            T = this.STATVAR.T .+ 273.15;
            saltConc = this.STATVAR.saltConc;
        end

        theta_sat = (this.STATVAR.porosity[1:end-1] .+ this.STATVAR.porosity[2:end]) ./2.0;

        #see del'amico
        alpha = this.PARA.alpha;         #in Van genuchten model [m^-1]
        n = this.PARA.n;                 #in van genuchten model [-]


        R = this.CONST.R;
        g = this.CONST.g;
        rho_w = this.CONST.rho_w;
        L_f = this.CONST.L_f;
        Tmelt = this.CONST.Tmelt; # in K

        Tmelt = Tmelt .+ Tmelt ./ L_f .* ( .- R .* saltConc .* Tmelt); #in Kelvin

        #Eq. 17 from Dall Amico, but using the solute potential of the UNFROZEN
        #soil as the last term - I think this is correct now consistent with
        #freezing_point_depression in free water

        water_pot = L_f ./ (rho_w .* g .* Tmelt) .* (T .- Tmelt) .* (T .< Tmelt);
        #the first term is equivalentto Eq. 2 in Dall Amico - this is the effect of the matric potential
        #the third term is the solute potential of the unfrozen soil - somehow this is needed as a reference potential
        #the sethermCond term is the solute potential of the freezing soil - this is where the magig happens
        #as it increases with decreasing liqWater - but not sure why this is PLUS -
        #but it must be plus to offset the efect of the first term which is minus
        #and becoming more negative with decreasing T - otherwise the freezing
        #happens even faster

        water_pot = water_pot .* (water_pot .< 0.0);
        #this is needed to get unchanged water contents above Tmelt

        liqWater = theta_sat ./ (1.0 .+ (alpha[1] .* abs.(water_pot)) .^n[1]) .^(1.0 .- 1.0 ./n[1]);
        #this is the final formula, but waterPot depends also on liqWater, so it's a
        #coupled eq. system

        return liqWater, Tmelt
    end

    function getDerivativeWaterContent(this)
        #Within this function, T is in Kelvin

        T = this.STATVAR.T + 273.15; #convert from C to Kelvin
        saltConc = this.STATVAR.saltConc;
        Tmelt = this.STATVAR.Tmelt; #in Kelvin

        delta = this.PARA.delta; #timestep for derivative
        R = this.CONST.R;
        L_f = this.CONST.L_f;
        Tmelt_free_water = this.CONST.Tmelt;

        saltConc_unfrozen = - 0.5 .* (T .- Tmelt_free_water) .* L_f ./ R ./ Tmelt_free_water .^2.0;

        #derivative wrt salt content
        saltConc_left_bound  = minimum(saltConc .- delta, saltConc_unfrozen .- delta);
        saltConc_right_bound = minimum(saltConc .+ delta, saltConc_unfrozen);
        dliqWater_dsaltConc = (saltConc .< saltConc_unfrozen) .* (watercontent_temperature_salt(this, T, saltConc_left_bound) .- watercontent_temperature_salt(this, T, saltConc_right_bound)) ./(saltConc_left_bound .- saltConc_right_bound);

        #derivative wrt temperature
        T_left_bound = minminimum(T .- delta, Tmelt .- delta);
        T_right_bound = minimum(T .+ delta, Tmelt);
        dliqWater_dT = (T < Tmelt) .* (watercontent_temperature_salt(this, T_left_bound, saltConc) .- watercontent_temperature_salt(this, T_right_bound, saltConc)) ./ (T_left_bound .- T_right_bound);

        return dliqWater_dT, dliqWater_dsaltConc
    end


end
