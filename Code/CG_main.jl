 # CG_main
# Main file of CryoGrid, julia version
#
# Based on the matlab version by S.Westermann
#
#
#

module CG_main
# load used modules and stratum
#include("StrataModules/SEDIMENT_T.jl")
#include("ForcingModules/Forcing_inundationHistory.jl")
#include("Common/CryoGridTypes.jl")
#include("Common/CryoGridInitialization.jl")
#include("Common/CryoGridInteractions.jl")
#include("Results/OUT_subseaPF.jl")

function main(TOP, BOTTOM, forcing, out, savename)
# run time loop
t = forcing.PARA.start_time;

CURRENT = TOP.NEXT;

#t is in days, timestep should also be in days
while t <= forcing.PARA.end_time
    #global forcing
    #global CURRENT
    #global TOP
    #global BOTTOM
    #global t
    #global out
    #global savename


    forcing = forcing.interpolate_forcing(t, forcing);
    #---------boundary conditions
    #proprietary function for each class, i.e. the "real upper boundary"
    #only evaluated for the first cell/block

    TOP.NEXT.get_boundary_condition_u(TOP.NEXT, forcing);
    CURRENT = TOP.NEXT;

    #function independent of classes, each class must comply with this function
    #evaluated for every interface between two cells/blocks
    while ~isequal(CURRENT.NEXT, BOTTOM)
        CURRENT.IA_NEXT.get_boundary_condition_m(CURRENT.IA_NEXT);
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
    timestep = min(timestep, (out.TEMP.out_time[1] - t[1])); #in days!!
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
    out = out.store_out(out, t, TOP, BOTTOM, forcing, savename);

    #calculate new time
    t = t .+ timestep; #./day_sec;

    #println("current ground temperature")
    #println(TOP.NEXT.STATVAR.T[1])
    if out.BREAK == true
        t = forcing.PARA.end_time + [1.0] ;
    end

end
end
end #module
