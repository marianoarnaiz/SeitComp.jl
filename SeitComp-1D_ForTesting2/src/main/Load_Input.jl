"This is a function Helps load the inputs into memory"

function Load_Input(FILE::String)
        try
                include(FILE);
                #SWF
                global Obs_T=Obs_SWD[:,1]; #Observed periods (s)
                global Nmodes=Int((size(Obs_SWD,2)-1)/4); #Number of modes (Currently 1 to 3)
                # For P wave Receiver functions
                global depth_P= Obs_PRF[2,2:end];
                global dist_P=Obs_PRF[1,2:end];
                # For S wave Receiver functions
                global depth_S= Obs_SRF[2,2:end];
                global dist_S=Obs_SRF[1,2:end];
                ## Inputs for SKS wave Receiver functions
                global depth_SKS= Obs_SKSRF[2,2:end];
                global dist_SKS=Obs_SKSRF[1,2:end];
                # ## Compute from the initial conditions:
                # # zTPC = Depth [km], Temperature [K], Preassure [kbar], Composition  [%wt]
                # # Model = Depth [km], Vp [km/s], Vs [km/s], Rho[g/cm^3], etc. A Model of the Earth to be used in the Seismic part
                # # qv = Surface Heatflow
                # zTPC,Model,qv=Make_Model(Temperature_Nodes,Temp_Depth_Nodes,Depths,LAB_index,Mantle_Composition);

                printstyled(" Input for modelling correcly loaded!",color=:green)
        catch e2
                printstyled(" An error ocurred. Please check Input_Template.jl file",color=:red)
        end
return
end
