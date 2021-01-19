module OUT_inversion
#this_SUBSEAPF Manage and save output data for subsea permafrost
#manages a struct with forcing data, info on the performance of the
#current run and stores the results at certain times
#saves the accumulated results at the final time or if the run breaks
include("../Common/CryoGridTypes.jl")
include("../Common/matlab.jl")
using Statistics
using MAT
#using Libc
#using Filesystem

mutable struct OUT

    RES::CryoGridTypes.outresults
    TEMP::CryoGridTypes.outtemporary
    PARA::CryoGridTypes.outparameter
    FORCING#::CryoGridTypes.forcingdata
    BREAK::Bool
    RUNINFO::CryoGridTypes.runinfo #save: anzahl echte zeitschritte, rechenzeit (gesamt)
    #provide_variables::Function
    initalize_from_file::Function
    complete_init_out::Function
    store_out::Function

    function OUT()
        this = new()

        this.initalize_from_file = function(this::OUT)
            this.TEMP = CryoGridTypes.outtemporary([0.0], [0.0], [1.0]);

            this.PARA = CryoGridTypes.outparameter([0.0], [0.0], [0.0], [0.0]);

            this.BREAK = false;
            this.RUNINFO = CryoGridTypes.runinfo([time()], [0.0], [1000.0], [0.0], [0.0])
            return this
        end

        this.complete_init_out = function(this::OUT, forcing, depthInterp=vec(0:0.1:10), output_timestep=7)

            start_saving = max(forcing.PARA.startForcing[1], -365.25*5.0);
            this.TEMP.out_time .= start_saving; #start saving on the first time step
            # this value is in days in the forcing!
            #save a maximum of 5 years
            this.PARA.save_time .= forcing.PARA.endForcing[1];
            this.PARA.output_timestep .= output_timestep; #days

            #save forcing data
            this.FORCING = deepcopy(forcing.PARA);
            timestamp = cat([this.TEMP.out_time[1]:this.PARA.output_timestep[1]:this.PARA.save_time[1];], this.PARA.save_time[1], dims=1); #years for saving
            emptyRes = NaN64 .* ones(length(depthInterp), length(timestamp));

            this.RES = CryoGridTypes.outresults(timestamp,  depthInterp, deepcopy(emptyRes), deepcopy(emptyRes), deepcopy(emptyRes), deepcopy(emptyRes), deepcopy(emptyRes));

            return this
        end

        this.store_out = function(this::OUT, t, TOP, BOTTOM, forcing, savename::String)

            #update runinfo every step
            this.RUNINFO.timesteps[1] += 1.0;
            #this.RUNINFO.dt_min = min(this.RUNINFO.dt_min, run_info.current_timestep) ;
            #this.RUNINFO.dt_max = max(this.RUNINFO.dt_max, run_info.current_timestep);

            if t==this.TEMP.out_time
                out_index = Int64(this.TEMP.out_index[1]);

                T = Array{Float64,1}()
                saltConc = Array{Float64,1}()
                thermCond = Array{Float64,1}()
                c_eff = Array{Float64,1}()
                liqWater = Array{Float64,1}()
                midpointDepth = Array{Float64,1}()
                layerDepth = Array{Float64,1}()

                CURRENT = TOP.NEXT;
                while ~isequal(CURRENT, BOTTOM)
                    #global CURRENT

                    T = cat(T, CURRENT.STATVAR.T, dims=1);
                    if forcing.PARA.saltForcingSwitch[1] == 0 #no salt diffusion
                        saltConc = cat(saltConc, CURRENT.STATVAR.saltConc./CURRENT.STATVAR.porosity .*CURRENT.STATVAR.liqWater, dims=1);
                    else
                        saltConc = cat(saltConc, CURRENT.STATVAR.saltConc, dims=1);
                    end
                    thermCond = cat(thermCond, CURRENT.STATVAR.thermCond, dims=1);
                    c_eff = cat(c_eff, CURRENT.STATVAR.c_eff, dims=1);
                    liqWater = cat(liqWater, CURRENT.STATVAR.liqWater, dims=1);

                    layers = [CURRENT.STATVAR.upperPos; CURRENT.STATVAR.upperPos .- cumsum(CURRENT.STATVAR.layerThick)];
                    layerDepth = cat(layerDepth, layers, dims=1);

                    midpoints = CURRENT.STATVAR.upperPos .+ CURRENT.STATVAR.layerThick[1] / 2.0 .- cumsum(CURRENT.STATVAR.layerThick);
                    midpointDepth = cat(midpointDepth, midpoints, dims=1);

                    CURRENT = CURRENT.NEXT;
                end

                this.RES.T[:,out_index] = reverse(matlab.interp1(reverse(midpointDepth), reverse(T), reverse(this.RES.depthInterp), "linear", NaN))
                this.RES.saltConc[:,out_index] = reverse(matlab.interp1(reverse(midpointDepth), reverse(saltConc), reverse(this.RES.depthInterp), "linear", NaN))
                this.RES.thermCond[:,out_index] = reverse(matlab.interp1(reverse(layerDepth), reverse(thermCond), reverse(this.RES.depthInterp), "linear", NaN))
                this.RES.c_eff[:,out_index] = reverse(matlab.interp1(reverse(midpointDepth), reverse(c_eff), reverse(this.RES.depthInterp), "linear", NaN))
                this.RES.liqWater[:,out_index] = reverse(matlab.interp1(reverse(midpointDepth), reverse(liqWater), reverse(this.RES.depthInterp), "linear", NaN))

                if t==this.PARA.save_time || this.BREAK == true
                    #if savename is empty, don't save the results

                    this.RUNINFO.endtime .= time() .- this.RUNINFO.starttime;

                    if ~(savename=="")
                        folder = "Results/" * savename;

                        mkpath(folder);

                        println("Save in folder " * folder )
                        file = savename * ".mat";
                        saveloc = folder * "/" * file;
                        f = matopen(saveloc, "w")
                        write(f, "RESULTS", this.RES)
                        write(f, "RUNINFO", this.RUNINFO)
                        write(f, "PARA", this.PARA)
                        write(f, "FORCING", this.FORCING)
                        close(f)
                    end
                    #don't save anything if there is no savename, just give back the results
                    this.BREAK = true;
                else
                    this.TEMP.out_index[1] += 1.0;
                    this.TEMP.out_time .= this.RES.time[Int64(this.TEMP.out_index[1])];
                end

            end
            return this
        end
        return this
    end
end
end
