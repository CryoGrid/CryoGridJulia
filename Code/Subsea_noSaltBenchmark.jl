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

lat_list = [69.945, 73.635, 73.68, 74.985];
lon_list = [-134.0, 117.0, 117.0, 117.0];

@inbounds for i=1:4

    # set run variables
    savename = "submarine_Benchmark_Testlocation" * string(i);

    # forcing variables are currently set in the focing class
    # initialize forcing
    forcing = Forcing_inundationHistory.forcing();
    forcing = forcing.initialize(forcing, 3, 2, lat_list[i], lon_list[i], -450000.0)
    forcing = forcing.load_forcing_from_mat(forcing);

    #initialize out
    out = OUT_subseaPF.OUT();
    out = out.initalize_from_file(out);
    out = out.complete_init_out(out, forcing);


    # initialize stratigraphy
    marine = SEDIMENT_T.stratum();

    # sediment has fields soilType, mineral, organic, salinity
    saltysilt = CryoGridTypes.sediment([1], [0.6], [0.0], [890.0]);

    grid = cat([0.0:0.5:50.0;], [51.0:1.0:200.0;], [202.0:2.0:500.0;], [505.0:5.0:1300;], [1350.0:50.0:3000.0;], [3500.0:500.0:6000.0;], dims=1)
    layerThick = grid[2:end] - grid[1:end-1]

    marine.initialize(marine, saltysilt, 0.0, -6000.0, layerThick);

    TOP = CryoGridInitialization.Top();
    TOP = TOP.init_top(TOP, marine); # after this TOP.NEXT = marine

    BOTTOM = CryoGridInitialization.Bottom();
    BOTTOM = BOTTOM.init_bottom(BOTTOM, marine); # after this BOTTOM.PREVIOUS = marine

    TOP.NEXT.NEXT= BOTTOM;
    BOTTOM.PREVIOUS.PREVIOUS = TOP;


    forcing = forcing.interpolate_forcing(forcing.PARA.start_time, forcing);
    TOP.NEXT.get_boundary_condition_u(TOP.NEXT, forcing);


    #TOP_CLASS = add_CHILD_snow(TOP_CLASS, class_list, stratigraphy_list);
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT, BOTTOM)
        #global CURRENT
        CURRENT.initialize_statvar(CURRENT);
        CURRENT = CURRENT.NEXT;
    end


    @time CG_main.main(TOP, BOTTOM, forcing, out, savename)
    #@code_warntype CG_main.main(TOP, BOTTOM, forcing, out, savename)
end
