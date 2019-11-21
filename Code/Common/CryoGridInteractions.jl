module CryoGridInteractions

    mutable struct IA_HEAT

        PREVIOUS
        NEXT

        get_boundary_condition_m::Function

        function IA_HEAT()
            this = new();

            this.PREVIOUS = [];
            this.NEXT = [];

            this.get_boundary_condition_m = function(this::IA_HEAT)
                #get the two classes involved
                upperStratum = this.PREVIOUS;
                lowerStratum = this.NEXT;

                #calculate the heat flux based on the temperature difference at the border
                flux = (upperStratum.STATVAR.T[end] .- lowerStratum.STATVAR.T[1]) .* upperStratum.STATVAR.thermCond[end] .* lowerStratum.STATVAR.thermCond[1] ./(upperStratum.STATVAR.thermCond[end] .* lowerStratum.STATVAR.layerThick[1] ./ 2.0 .+ lowerStratum.STATVAR.thermCond[1] .* upperStratum.STATVAR.layerThick[end] ./ 2.0 );

                #assign the flux to the lower/upper boundary of the upper/lower class, respectively
                upperStratum.TEMP.heatFlux_lb = [-flux];
                lowerStratum.TEMP.heatFlux_ub = [-flux];

                #add temperature at upper boundary of the lower class for conductivity calculations
                lowerStratum.TEMP.T_ub = [upperStratum.STATVAR.T[end]];

                return this
            end
            return this
        end
    end
end
