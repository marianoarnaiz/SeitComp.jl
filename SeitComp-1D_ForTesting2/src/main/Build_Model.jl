"This a function to only compute the properties model from the input.
it will not run, but it will give you an idea of the model that will
be used for the forward computation. I hope this saves you some time"

function Build_Model()
        try
                ## Compute from the initial conditions:
                # zTPC = Depth [km], Temperature [K], Preassure [kbar], Composition  [%wt]
                # Model = Depth [km], Vp [km/s], Vs [km/s], Rho[g/cm^3], etc. A Model of the Earth to be used in the Seismic part
                # qv = Surface Heatflow
                global zTPC,Model,qv=Make_Model(Temperature_Nodes,Temp_Depth_Nodes,Depths,LAB_index,Mantle_Composition);
                printstyled(" Geophysical Model correctly created!",color=:green)
                println(" ")
                printstyled("   🔥   Surface Heat Flow:  $qv mW/m^2 ",color=:light_red)
                println(" ")
                ## Compute the Qs for the T and P
                global Qsv,J1v,J2v=Qs_Burgers(zTPC,Grain_Size_Model);
                Model[:,5]=Qsv;
                ## Compute Melt Fraction
                global Melt,Latent_Heat,T_Solidus,T_Liquidus=Melt_fraction(zTPC,Melt_Crust,Melt_Mantle)
##              ## Correct Properties/Model for Qs and Melt Fraction
                ## We need to Correct the Geophysical Model (Variable Model) for
                # the effect of Melt and Attenuation. We will do it here so this
                # function does is not just a calling function!
                ## Correction for Melt
                if sum(Melt) == 0.0
                        printstyled("   🪨   No melt was found in the model.",color=:light_blue)
                        println(" ")
                else
                        # # Correct the Model for melt: Vp,Vs,Rho and Qs are corrected here
                        # global Model=Melt_Model(zTPC,Model,Melt)
                        # # Tell me where there is melt
                        # melt_around = mean(zTPC[findall(Melt .> 0),1]);
                        maxmelt=maximum(Melt);

                        println(" ")
                        printstyled("   🌋   Melt was found in the model!",color=:light_red)
                        println(" ")
                        printstyled("       Max melt= $maxmelt %",color=:light_red)
                        println(" ")
                end
##              ## Apply effect of attenuation to all the model!
                Model=Attenuate_Model(zTPC,Model)
        catch e2
                printstyled(" An error ocurred. Please check your Input_Template.jl",color=:red)
        end
return
end
