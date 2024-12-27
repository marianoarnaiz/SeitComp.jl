"This is a function Helps run all the forward calculations"

function Run0()
        try
                ## Compute SWF Dispersion
                #SWF
                T, SWV,EG_FUNC_R_Z, EG_FUNC_R_R, EG_FUNC_L_T = SW_Disp(Model,FEM_MESH,Obs_T,Nmodes,Earth,Elastic,Anisotropy);

                ## Compute P wave RF.
                #P RF
                P_matrix=zeros(2048,size(dist_P,1)+1);
                for i=1:size(dist_P,1)
                time_P, P_rf, Used_model = P_RF_forward(Model,depth_P[i],dist_P[i],Gaussian_factor,Clean_Model,rotate_P);
                P_matrix[:,1]=time_P;
                P_matrix[:,i+1]=P_rf;
                end

                ## Compute S wave RF.
                S_matrix=zeros(2048,size(dist_S,1)+1);
                try
                for i=1:size(dist_S,1)
                time_S, S_rf, Used_model = S_RF_forward(Model,depth_S[i],dist_S[i],Gaussian_factor,Clean_Model);
                S_matrix[:,1]=time_S;
                S_matrix[:,i+1]=S_rf;
                end
                catch e
                    println(" ")
                    printstyled(" S to P RF are not compatible with the Model.",color=:yellow)
                    println(" ")
                    printstyled(" You might need to change lines 47 and 52",color=:yellow)
                    println(" ")
                    printstyled(" in src/seis/S_RF_forward.jl to make it compatible.",color=:yellow)
                    println(" ")
                    printstyled(" Here a cutoff depth is imposed (250 km)",color=:yellow)
                    println(" ")
                end


                ## Compute SKS wave RF.
                #SKS RF
                SKS_matrix=zeros(4096,size(dist_SKS,1)+1);
                for i=1:size(dist_SKS,1)
                time_SKS, SKS_rf, Used_model = SKS_RF_forward(Model,depth_SKS[i],dist_SKS[i],Gaussian_factor,Clean_Model);
                SKS_matrix[:,1]=time_SKS;
                SKS_matrix[:,i+1]=SKS_rf;
                end
        catch e2
                printstyled(" An error ocurred. Please check Input_Template.jl file",color=:red)
        end
return
end
