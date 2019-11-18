module CryoGridInputDataSets

    using MAT

    function getElevation(Lat::Array{Float64,1}, Lon::Array{Float64,1})
        file = matopen("../ForcingData/subsea6p25kmCellSize_Shelf.mat")
        data = read(file, "LL6p25km");

        zsb = zeros(size(Lat));

        @inbounds for i = 1:length(Lat)
            distances = sqrt.((data[:,1] .- Lat[i]).^2.0 + (data[:,2] .- Lon[i]).^2.0);
            min_distance, this_ind = findmin(distances);
            zsb[i] = data[this_ind,3];
        end

        return zsb
    end

    function getQ_Davies(Lat::Array{Float64,1}, Lon::Array{Float64,1})
        #check that lat and lon are equal in length
        if size(Lat)!=size(Lon)
            return NaN
        end
        Qout = zeros(size(Lat));

        #load geothermal data
    	file = matopen("../ForcingData/Data_HFip.mat") #geothermal heat flux data from Davis et al. (2013, Q3) - data are defined on the same grid as CLIMBER2 SAT and ice sheet data (lon x lat) (240x41)
        lat_HFip = read(file, "lat_HFip");
        lon_HFip = read(file, "lon_HFip");
        HFip = read(file, "HFip");
        #find indices
        @inbounds for i = 1:length(Lat)
            #check this in julia!
            minlatdist, indLat = findmin((abs.(lat_HFip .- Lat[i])));
            minlondist, indLon = findmin((abs.(lon_HFip .- Lon[i])));
            Qout[i] = HFip[indLon[2],indLat[2]] ./ 1000.0;
        end
        return Qout
    end

end
