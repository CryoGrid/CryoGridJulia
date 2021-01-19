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
include("Results/OUT_subseaPF.jl")
include("CG_main.jl")


# set run variables
savename = "testlocation1_50k"

# forcing variables are currently set in the focing class
# initialize forcing
forcing = Forcing_inundationHistory.forcing();
forcing = forcing.initialize(forcing, 666, 2, 69.945, -134.0, -50000.0)
forcing = forcing.load_forcing_from_mat(forcing);

#initialize out
out = OUT_subseaPF.OUT();
out = out.initalize_from_file(out);
out = out.complete_init_out(out, forcing);


# initialize stratigraphy
terrestrial = SEDIMENT_T.stratum();
marine = SEDIMENT_T.stratum();

# sediment has fields soilType, mineral, organic, salinity
saltysilt = CryoGridTypes.sediment([1], [0.6], [0.0], [890.0]);
sand = CryoGridTypes.sediment([0], [0.6], [0.0], [0.0]);

marine.initialize(marine, saltysilt, 0.0, -250.0, 2.0);
terrestrial.initialize(terrestrial, sand, -250.0, -2000.0, 2.0);

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
TOP.NEXT.PREVIOUS = TOP;
BOTTOM.PREVIOUS.NEXT = BOTTOM;

TOP.NEXT.IA_PREVIOUS = [];
TOP.NEXT.IA_NEXT = interaction_heat;
TOP.NEXT.IA_NEXT.PREVIOUS = TOP.NEXT;
TOP.NEXT.IA_NEXT.NEXT = TOP.NEXT.NEXT;
BOTTOM.PREVIOUS.IA_PREVIOUS = TOP.NEXT.IA_NEXT;
BOTTOM.PREVIOUS.IA_NEXT = [];


#TOP_CLASS = add_CHILD_snow(TOP_CLASS, class_list, stratigraphy_list);
CURRENT = TOP.NEXT;
while ~isequal(CURRENT, BOTTOM)
    global CURRENT
    CURRENT.initialize_statvar(CURRENT); # this uses this.TEMP.T_ub but it's still zero
    CURRENT = CURRENT.NEXT;
end



@time CG_main.main(TOP, BOTTOM, forcing, out, savename)
#@code_warntype CG_main.main(TOP, BOTTOM, forcing, out, savename)
