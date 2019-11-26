module CryoGridConstants

    #reads in global constants
    function getGlobalConst(stratum)

        stratum.CONST.day_sec .= [86400.0];	#[sec]	number of seconds in a day
        stratum.CONST.year_sec .= [31536000.0]; #[sec]	number of seconds in a year

        stratum.CONST.rho_w .= [1000.0];      #[kg/m3]	density water
        stratum.CONST.rho_i .= [920.0];   	#[kg/m3]	density ice

        stratum.CONST.L_f .= [3.34e8];      #[J/m3]	volumetric latent heat of fusion, freezing
        stratum.CONST.L_v .= [0.0];    		#[J/m3] volumetric latent heat of vaporization

        stratum.CONST.c_i .= [1.9e6];       #[J/m3K]	volumetric heat capacity ice
        stratum.CONST.c_w .= [4.2e6];       #[J/m3K]	volumetric heat capacity water
        stratum.CONST.c_o .= [2.5e6];       #[J/m3K]	volumetric heat capacity organic content
        stratum.CONST.c_m .= [2.0e6];       #[J/m3K]	volumetric heat capacity mineral content
        stratum.CONST.c_a .= [1005.0];      #[J/m3K]	volumetric heat capacity air

        stratum.CONST.k_a .= [0.025];       #[W/mK]	thermal conductivity air [Hillel(1982)]
        stratum.CONST.k_i .= [2.2];         #[W/mK]	thermal conductivity ice [Hillel(1982)]
        stratum.CONST.k_w .= [0.57];        #[W/mK]	thermal conductivity water [Hillel(1982)]
        stratum.CONST.k_o .= [0.25];        #[W/mK]	thermal conductivity organic [Hillel(1982)]
        stratum.CONST.k_m .= [3.0];         #[W/mK]	thermal conductivity mineral - default value, should be changed if necessary

        stratum.CONST.sigma .= [5.67e-8];   #[W/m2K4] Stefan-Boltzmann const.

        stratum.CONST.R	.= [8.314459];	     #[J/K mol]	universal gas constant
        stratum.CONST.R_a .= [287.058];     #gas constant
        stratum.CONST.kappa .= [0.4];       #karman constant

        stratum.CONST.Tmelt .= [273.15];	#[K]	freezing temperature of free water

        stratum.CONST.g .= [9.81];          #[m/sec2]	gravitational acceleration Earth surface

        stratum.CONST.benthicSalt .= [989.63]; #[mol/m3]	moles of ions, Millero et al. (2008)
        stratum.CONST.tau .= [1.5];	        #[-] tortuosity %1.5 standard

        return stratum
    end

    #reads in global parameter
    function getGlobalPara(stratum)
        stratum.PARA.heatFlux_lb    .= [0.05];
        stratum.PARA.albedo  .= [0.15];
        stratum.PARA.epsilon .= [0.99];
        stratum.PARA.rs      .= [100.0];
        stratum.PARA.z       .= [2.0]; #measurement height [m]
        stratum.PARA.z0      .= [1.0e-3]; #roughness length [m]

        stratum.PARA.max_timestep .= [3600.0]; #[sec]
        stratum.PARA.dE_max  .= [0.5e5]; #[J/m3]
        stratum.PARA.dt_max  .= [100.0];
        stratum.PARA.dsaltConc_max .= [5.0e-3];

        stratum.PARA.alpha   .= [6.5e-1]; #6.50E-01 for silt, 4.06 for sand
        stratum.PARA.n       .= [1.67]; # 1.67 for silt, 2.03 for sand
        stratum.PARA.delta   .= [1.0e-7];
        stratum.PARA.saltDiff0 .= [8.0e-10];


        return stratum
    end

end
