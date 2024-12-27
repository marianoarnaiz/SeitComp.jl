"This is a function Helps run all the forward calculations"

function Write_Outs()
        try
                writedlm("results/Model.txt", round.(Model,digits=5));
                writedlm("results/Cal_PRF.txt", round.([[0 dist_P'] ; [0 depth_P'] ; P_matrix],digits=5));
                writedlm("results/Cal_SRF.txt", round.([[0 dist_S'] ; [0 depth_S'] ; S_matrix],digits=5));
                writedlm("results/Cal_SKSRF.txt", round.([[0 dist_SKS'] ; [0 depth_SKS'] ; SKS_matrix],digits=5));
                writedlm("results/Cal_Stacked_PRF.txt", round.(Stacked_PRF,digits=5));
                writedlm("results/Cal_Stacked_SRF.txt", round.(Stacked_SRF,digits=5));
                writedlm("results/Cal_Stacked_SKSRF.txt", round.(Stacked_SKSRF,digits=5));


                # fis SWV for writting
                if Nmodes==1
                        writedlm("results/Cal_SWD.txt", round.([T SWV[:,:,1]],digits=5));
                end
                if Nmodes==2
                        writedlm("results/Cal_SWD.txt", round.([T SWV[:,:,1] SWV[:,:,2]],digits=5));
                end
                if Nmodes==3
                        writedlm("results/Cal_SWD.txt", round.([T SWV[:,:,1] SWV[:,:,2] SWV[:,:,3]],digits=5));
                end
                if Nmodes==4
                        writedlm("results/Cal_SWD.txt", round.([T SWV[:,:,1] SWV[:,:,2] SWV[:,:,3] SWV[:,:,4]],digits=5));
                end
                if Nmodes==5
                        writedlm("results/Cal_SWD.txt", round.([T SWV[:,:,1] SWV[:,:,2] SWV[:,:,3] SWV[:,:,4] SWV[:,:,5]],digits=5));
                end



                #let me fix zTPC for you
                for i in eachindex(zTPC)
                    if !isassigned(zTPC, i)
                        zTPC[i] = -99999
                    end
                end
                writedlm("results/zTPC.txt", zTPC);
                printstyled(" Files are in results/ !",color=:green)
        catch e20
                printstyled(" An error ocurred. Please re-run all the computations file",color=:red)
        end
return
end
