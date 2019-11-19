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

        Tmelt = this.CONST.Tmelt ./ this.CONST.L_f .* (-this.CONST.R .* this.STATVAR.saltConc.*this.CONST.Tmelt);

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

            a1 = 1.0 ./ this.STATVAR.porosity[i] - this.PARA.a[i] .* Tmelt[i] - this.PARA.b[i] .* Tmelt[i].^2.0;

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
        thermCond = zeros(length(this.STATVAR.layerThick)+1,1);
        thermCond = dropdims(thermCond, dims=2);

        #calculate for upper boundary
        i = 1;
        a1 = 1.0 ./ porosity[i] - a[i]*abs(Tmelt[i])^b[i];

        if (T_ub[1] .+ T[i]) / 2.0 <= Tmelt_inDegreeC[i]
            liqWater = 1.0 / (a1 .+ a[i]*abs((T_ub .+ T[i]) / 2.0) .^ b[i]);
        else
            liqWater = porosity[i];
        end

        thermCond[i] = (liqWater[1] .* sqrt.(k_water[1]) .+ (porosity[i] .- liqWater[1]) .* sqrt.(k_ice[1]) .+ mineral[i] .* sqrt.(k_mineral[1]) .+ organic[i] .* sqrt.(k_organic[1])) .^2.0;

        #calculate for inner edges
        @inbounds for i = 2:length(layerDepth) - 1
            #take the average of liqWater's of grid cells above and below
            a1 = 1.0 / porosity[i-1] - a[i-1]*abs(Tmelt[i-1])^b[i-1];
            if T[i-1] <= Tmelt_inDegreeC[i-1]
                liqWater = 1.0 / (a1 + a[i-1]*abs(T[i-1])^b[i-1]);
            else
                liqWater = porosity[i-1];
            end

            a1 = 1.0 / porosity[i] - a[i]*abs(Tmelt[i])^b[i];
            if T[i] <= Tmelt_inDegreeC[i]
                liqWater = (liqWater + 1.0 / (a1 + a[i]*abs(T[i])^b[i])) / 2.0;
            else
                liqWater = (liqWater+porosity[i]) / 2.0;
            end

            thermCond[i] = (liqWater[1] .* sqrt.(k_water[1]) .+ ((porosity[i] .+ porosity[i-1]) ./ 2.0 .- liqWater[1]) .* sqrt.(k_ice[1]) .+  (mineral[i-1] .+ mineral[i]) ./ 2.0 .* sqrt.(k_mineral[1]) .+ (organic[i-1] .+ organic[i]) ./ 2.0 .* sqrt.(k_organic[1])) .^2.0;
        end

        #extrapolate for lower edge - this is needed for calculation of
        #interaction
        thermCond[i+1] = thermCond[i];

        #------------------ determine bulk conductivity and capacity------%
        # also save liqWater here
        #Heat Capacity lives on the midpoints of the cells
        c_eff = zeros(size(this.STATVAR.layerThick));
        liqWater = zeros(size(this.STATVAR.layerThick));
        midptDepth = this.STATVAR.upperPos .- this.STATVAR.layerThick[1] ./ 2.0 .- cumsum(this.STATVAR.layerThick);

        @inbounds for i = 1:length(midptDepth)
            a1 = 1.0 / porosity[i] - a[i]* abs(Tmelt[i])^b[i];

            #determine water content
            if T[i] <= Tmelt_inDegreeC[i]
                liqWater[i] = 1.0 / (a1 + a[i]*abs(T[i])^b[i]);
            else
                liqWater[i] = porosity[i];
            end

            if T[i] <= Tmelt_inDegreeC[i]
                println("T[i]")
                println(T[i])
                println("a[i]")
                println(a[i])
                println("b[i]")
                println(b[i])
                d_liqWater = a[i] * b[i] * abs(T[i]^(b[i] - 1.0) / (a1 + a[i]*abs(T[i])^b[i])) / (a1 + a[i]*abs(T[i])^b[i]);
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

end
