 # CG_main.jl
#
# Main file of CryoGrid, julia version
#
# Based on the matlab version by S.Westermann
#
#
#
# structure of this file
# 1. load modules and stratum
# 2. set run variable
# 3. initialize forcing
# 4. initialize stratigraphy
# 5. run time-loop
# 6. save output


# load used modules and stratum
include("StrataModules/SEDIMENT_T.jl")
include("ForcingModules/Forcing_inundationHistory.jl")
include("Common/CryoGridTypes.jl")
include("Common/CryoGridInitialization.jl")
include("Common/CryoGridInteractions.jl")


# set run variables

# forcing variables are currently set in the focing class



# initialize forcing
forcing = Forcing_inundationHistory.forcing();
forcing = forcing.initialize(forcing);
forcing = forcing.load_forcing_from_mat(forcing);

# initialize stratigraphy
terrestrial = SEDIMENT_T.stratum();
marine = SEDIMENT_T.stratum();

# sediment has fields soilType, mineral, organic, salinity
saltysilt = CryoGridTypes.sediment([1], [0.6], [0.0], [890.0]);
sand = CryoGridTypes.sediment([0], [0.6], [0.0], [0.0]);

marine.initialize(marine, saltysilt);
terrestrial.initialize(terrestrial, sand);

interaction_heat = CryoGridInteractions.IA_HEAT();

#TODO!
#assemble the model stratigraphy and define interactions between classes
#TOP_CLASS, BOTTOM_CLASS, TOP, BOTTOM = assemble_stratigraphy(class_list, stratigraphy_list, grid, forcing);

TOP = CryoGridInitialization.Top();
TOP = TOP.init_top(TOP, marine);

BOTTOM = CryoGridInitialization.Bottom();
BOTTOM = BOTTOM.init_bottom(BOTTOM, terrestrial);

TOP.NEXT.NEXT= BOTTOM.PREVIOUS;
BOTTOM.PREVIOUS.PREVIOUS = TOP.NEXT;
TOP.NEXT.PREVIOUS = [];
BOTTOM.PREVIOUS.NEXT = [];

TOP.NEXT.IA_PREVIOUS = [];
TOP.NEXT.IA_NEXT = interaction_heat;
BOTTOM.PREVIOUS.IA_PREVIOUS = TOP.NEXT.IA_NEXT;
BOTTOM.PREVIOUS.IA_NEXT = [];


#TOP_CLASS = add_CHILD_snow(TOP_CLASS, class_list, stratigraphy_list);
CURRENT = TOP.NEXT;
while ~isequal(CURRENT.NEXT, BOTTOM)
    global CURRENT
    CURRENT.initialize_statvar(CURRENT);
    CURRENT = CURRENT.NEXT;
end

# run time loop
t = forcing.PARA.start_time;

CURRENT = TOP.NEXT;

#t is in days, timestep should also be in days
while t <= forcing.PARA.end_time

    forcing = forcing.interpolate_forcing(t, forcing);
    #---------boundary conditions
    #proprietary function for each class, i.e. the "real upper boundary"
    #only evaluated for the first cell/block

    TOP.NEXT.get_boundary_condition_u(TOP.NEXT, forcing);
    CURRENT = TOP.NEXT;

    #function independent of classes, each class must comply with this function
    #evaluated for every interface between two cells/blocks
    while ~isequal(CURRENT.NEXT, BOTTOM)
        CURRENT.get_boundary_condition_m(CURRENT.IA_NEXT);
        CURRENT = CURRENT.NEXT;
    end

    #proprietary function for each class, i.e. the "real lower boundary"
    #only evaluated for the last cell/block
    CURRENT.get_boundary_condition_l(CURRENT,  forcing);  #At this point, CURRENT is equal to BOTTOM_CLASS
    #--------------------------

    #calculate spatial derivatives for every cell in the stratigraphy
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT, BOTTOM)
        CURRENT.get_derivatives_prognostic(CURRENT);
        CURRENT = CURRENT.NEXT;
    end

    #calculate minimum timestep required for all cells in days
    CURRENT = TOP.NEXT;
    timestep = 10.0; #in days!
    while ~isequal(CURRENT, BOTTOM)
        timestep = min(timestep, CURRENT.get_timestep(CURRENT));
        CURRENT = CURRENT.NEXT;
    end
    #timestep = min(timestep, (out.OUTPUT_TIME-t).*day_sec);
    timestep = min(timestep, (out.OUTPUT_TIME-t)); #in days!!
    #make sure to hit the output times!

    #calculate prognostic variables
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT, BOTTOM)
        CURRENT.advance_prognostic(CURRENT, timestep);
        CURRENT = CURRENT.NEXT;
    end


    #calculate diagnostic variables
    #some effects only happen in the first cell
    #calculate bottom to top
    TOP.NEXT.compute_diagnostic_first_cell(TOP.NEXT, forcing);
    CURRENT = BOTTOM.PREVIOUS;
    while ~isequal(CURRENT, TOP)
        CURRENT.compute_diagnostic(CURRENT, forcing);
        CURRENT = CURRENT.PREVIOUS;
    end

    TOP_CLASS = TOP.NEXT; #TOP_CLASS and BOTTOM_CLASS for convenient access
    BOTTOM_CLASS = BOTTOM.PREVIOUS;

    #store the output according to the defined OUT clas
    out = store_OUT(out, t, TOP_CLASS, BOTTOM, forcing, run_number);

    #calculate new time
    t = t + timestep; #./day_sec;

    if out.BREAK == 1
        break
    end

end
