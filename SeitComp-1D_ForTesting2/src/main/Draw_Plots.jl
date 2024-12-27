s"This is a function Helps set up all the enviroment. JUST RUN ONCE!"

using DelimitedFiles, GMT, Plots
default(titlefont = ("times"), legend_font_family= ("times"), linewidth=2, yminorgrid = true)

function Draw_Plots()

    println(" ")
    printstyled(" Plotting Figures. Go to SeitComp-1D/figs",color=:green)
    println(" ")

    ## Fix some inputs
    #reload models just in case
    Model2=readdlm("data/PREM.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB
    ## Plot the model and all the properties
    LAB=argmin(abs.(zTPC[:,2].-Temperature_Nodes[2]));
    Cut=argmin(abs.(zTPC[:,1].-300));
    ## Get the Melt solidus and liquidus
    Z_SL=zTPC[:,1]
    T_Sol=T_Solidus;
    T_Liq=T_Liquidus;

    ## Compute Geophysical Model Anomalies
    iRefVp=Interpolations.interpolate((Model2[:,1],), Model2[:,2], Gridded(Linear()))
    iRefVs=Interpolations.interpolate((Model2[:,1],), Model2[:,3], Gridded(Linear()))
    iRefRho=Interpolations.interpolate((Model2[:,1],), Model2[:,4], Gridded(Linear()))

    RefVp=iRefVp(Model[:,1]);
    RefVs=iRefVs(Model[:,1]);
    RefRho=iRefRho(Model[:,1]);

    ## Plot the Geophysical Model
    p1=Plots.plot(Model[:,2], Model[:,1],yflip=true,xlim=(5,15),ylim=(-8,2890),label="Vp",lc=RGB(68/255, 1/255, 84/255), title="Geophysical Model")
    Plots.plot!(Model2[:,2], Model2[:,1],yflip=true,xlim=(5,15),ylim=(-8,2890),label="PREM Vp",lc=RGB(68/255, 1/255, 84/255),ls=:dot,xlabel="Vp (km/s)", ylabel="Depth (km)")
    p2=Plots.plot(Model[:,3], Model[:,1],yflip=true,xlim=(2.5,8),ylim=(-8,2890),label="Vs",lc=RGB(42/255, 120/255, 142/255), ylabel="Depth (km)")
    Plots.plot!(Model2[:,3], Model2[:,1],label="PREM Vs",lc=RGB(42/255, 120/255, 142/255),ls=:dot,xlabel="Vs (km/s)")
    p3=Plots.plot(Model[:,4], Model[:,1],yflip=true,xlim=(2,6),ylim=(-8,2890),label="Rho",lc=RGB(122/255, 209/255, 81/255))
    Plots.plot!(Model2[:,4], Model2[:,1],label="PREM Rho",lc=RGB(122/255, 209/255, 81/255),ls=:dot,minorgrid=true,xlabel="Rho (g/cm^3)", ylabel="Depth (km)")
    Plots.plot(p1,p2,p3, layout=(3,1),size=(800,1200),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
    savefig("figs/0a.-Geophysical_Model.pdf")

    ## Plot  the Geophysical Model Anomalies
    p4=Plots.plot(100*(Model[:,2]-RefVp)./RefVp, Model[:,1],yflip=true,xlim=(-5,5),ylim=(-8,2890),label="Vp",lc=RGB(68/255, 1/255, 84/255), title="Geophysical Anomalies")
    Plots.plot!(RefVp-RefVp, Model[:,1],yflip=true,xlim=(-5,5),ylim=(-8,2890),label="PREM Vp",lc=RGB(68/255, 1/255, 84/255),ls=:dot,xlabel="δVp (km/s)", ylabel="Depth (km)")
    p5=Plots.plot(100*(Model[:,3]-RefVs)./RefVs, Model[:,1],yflip=true,xlim=(-5,5),ylim=(-8,2890),label="Vs",lc=RGB(42/255, 120/255, 142/255), ylabel="Depth (km)")
    Plots.plot!(RefVs-RefVs, Model[:,1],label="PREM Vs",lc=RGB(42/255, 120/255, 142/255),ls=:dot,xlabel="δVs (km/s)")
    p6=Plots.plot(100*(Model[:,4]-RefRho)./RefRho, Model[:,1],yflip=true,xlim=(-5,5),ylim=(-8,2890),label="Rho",lc=RGB(122/255, 209/255, 81/255))
    Plots.plot!(RefRho-RefRho, Model[:,1],label="PREM Rho",lc=RGB(122/255, 209/255, 81/255),ls=:dot,minorgrid=true,xlabel="δRho (g/cm^3)", ylabel="Depth (km)")
    Plots.plot(p4,p5,p6, layout=(1,3),size=(800,1200),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
    savefig("figs/0d.-Geophysical_Anomalies.pdf")
    Plots.plot(p4,p5,p6,layout=(3,1),size=(800,1200),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
    savefig("figs/0c.-Geophysical_Anomalies.pdf")
    Plots.plot(p1,p4,p2,p5,p3,p6, layout=(3,2),size=(800,1200),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
    savefig("figs/0b.-Geophysical_Model.pdf")


    ## Plot a Model for all the mantle
    p1=Plots.plot(Model[:,2], Model[:,1],yflip=true,xlim=(2,16),ylim=(-8,2890),label="Vp",lc=RGB(68/255, 1/255, 84/255), title="Mantle Model")
    Plots.plot!(Model2[:,2], Model2[:,1],yflip=true,xlim=(2,16),ylim=(-8,2890),label="PREM Vp",lc=RGB(68/255, 1/255, 84/255),ls=:dot)
    Plots.plot!(Model[:,3], Model[:,1],label="Vs",lc=RGB(42/255, 120/255, 142/255))
    Plots.plot!(Model2[:,3], Model2[:,1],label="PREM Vs",lc=RGB(42/255, 120/255, 142/255),ls=:dot)
    Plots.plot!(Model[:,4], Model[:,1],label="Rho",lc=RGB(122/255, 209/255, 81/255))
    Plots.plot!(Model2[:,4], Model2[:,1],label="PREM Rho",lc=RGB(122/255, 209/255, 81/255),ls=:dot,minorgrid=true,xlabel="Vx (km/s) & Rho (g/cm^3)")

    p2=Plots.plot(zTPC[:,2], zTPC[:,1],yflip=true,xlim=(200,4500),ylim=(-8,2890),label="T (K)",lc=RGB(254/255, 176/255, 120/255),xlabel="T (K)",legend = :bottomleft)
      Plots.plot!(T_Sol, Z_SL,label=false,lc=RGB( 85/255, 107/255,  47/255),annotations = ([3160], [1500], "S."),annotationfontfamily="times",annotationfontsize=10)
      Plots.plot!(T_Liq, Z_SL,label=false,lc=RGB(70/255, 130/255, 180/255),annotations = ([3775], [1500], "L."),annotationfontfamily="times",annotationfontsize=10)
      pp = twiny()
      Plots.plot!(pp,zTPC[:,3]*1e-4, zTPC[:,1],yflip=true,xlim=(0,150),ylim=(-8,2890),label="P (GPa)",lc=RGB(0/255, 0/255, 4/255),legend = :topleft,xlabel="P (GPa)",grid=false)
      t1, t2 = 200, 4500;
      y = [Temp_Depth_Nodes[2]/1000]
      Plots.plot!([t1; t2], [y; y],label=false,lc=RGB(44/255, 17/255, 95/255),ls=:dash,annotations = ([4200], [y.-50], "LAB"),annotationfontfamily="times",annotationfontsize=10)
      y = [zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-14)),1]]
      Plots.plot!([t1; t2], [y; y],label=false,lc=RGB(114/255, 31/255, 129/255),ls=:dash,annotations = ([4200], [y.-50], "410 D."),annotationfontfamily="times",annotationfontsize=10)
      y = [zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-24)),1]]
      Plots.plot!([t1; t2], [y; y],label=false,lc=RGB(183/255, 55/255, 121/255),ls=:dash,annotations = ([4200], [y.-50], "660 D."),annotationfontfamily="times",annotationfontsize=10)
      y = mean([zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-120)),1],zTPC[argmin(abs.(zTPC[:,2] .-2500)),1]])
      Plots.plot!([t1; t2], [y; y],label=false,lc=RGB(221/255, 81/255, 58/255),ls=:dash,annotations = ([4200], [y.-50], "D''"),annotationfontfamily="times",annotationfontsize=10)

    p3=Plots.plot(Model[:,5], Model[:,1],yflip=true,xlim=(0,600),ylim=(-8,2890),label="Qs",lc=RGB(13/255, 8/255, 135/255),xlabel="Qs",legend = :bottomleft)
    pp2 = twiny()
    Plots.plot!(pp2,Melt, zTPC[:,1],yflip=true,xlim=(-1,3),ylim=(-8,2890),label="Melt (%)",lc=RGB(204/255, 71/255, 120/255),legend = :topright,xlabel="Melt (%)",grid=false)

    Plots.plot(p1,p2,p3, layout=(1,3),size=(1200,800),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
    savefig("figs/4.-Mantle_Model.pdf")

    ## Model of the Lithosphere
    p1=Plots.plot(Model[:,2], Model[:,1],yflip=true,xlim=(2,10),ylim=(-8,350),label="Vp",lc=RGB(68/255, 1/255, 84/255), title="Lithosphere Model",legend=:bottomleft)
    Plots.plot!(Model2[:,2], Model2[:,1],yflip=true,ylim=(-8,350),label="PREM Vp",lc=RGB(68/255, 1/255, 84/255),ls=:dot)
    Plots.plot!(Model[:,3], Model[:,1],label="Vs",lc=RGB(42/255, 120/255, 142/255))
    Plots.plot!(Model2[:,3], Model2[:,1],label="PREM Vs",lc=RGB(42/255, 120/255, 142/255),ls=:dot)
    Plots.plot!(Model[:,4], Model[:,1],label="Rho",lc=RGB(122/255, 209/255, 81/255))
    Plots.plot!(Model2[:,4], Model2[:,1],label="PREM Rho",lc=RGB(122/255, 209/255, 81/255),ls=:dot,minorgrid=true,xlabel="Vx (km/s) & Rho (g/cm^3)")

    p2=Plots.plot(zTPC[:,2], zTPC[:,1],yflip=true,xlim=(200,2500),ylim=(-8,350),label="T (K)",lc=RGB(254/255, 176/255, 120/255),xlabel="T (K)",legend = :bottomleft)
    Plots.plot!(T_Sol, Z_SL,label="Solidus",lc=RGB( 85/255, 107/255,  47/255))
    Plots.plot!(T_Liq, Z_SL,label="Liquidus",lc=RGB(70/255, 130/255, 180/255))
    pp = twiny()
    Plots.plot!(pp,zTPC[:,3]*1e-4, zTPC[:,1],yflip=true,xlim=(0,20),ylim=(-8,350),label="P (GPa)",lc=RGB(0/255, 0/255, 4/255),legend = :topright,xlabel="P (GPa)",grid=false)
    t1, t2 = 200, 4000;
    y = [Temp_Depth_Nodes[2]/1000]
    Plots.plot!([t1; t2], [y; y],label="LAB",lc=RGB(44/255, 17/255, 95/255),ls=:dash)
    y = [zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-14)),1]]
    Plots.plot!([t1; t2], [y; y],label="410 Disc.",lc=RGB(114/255, 31/255, 129/255),ls=:dash)
    y = [zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-24)),1]]
    Plots.plot!([t1; t2], [y; y],label="660 Disc.",lc=RGB(183/255, 55/255, 121/255),ls=:dash)
    y = mean([zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-120)),1],zTPC[argmin(abs.(zTPC[:,2] .-2500)),1]])
    Plots.plot!([t1; t2], [y; y],label="D''",lc=RGB(221/255, 81/255, 58/255),ls=:dash)

    p3=Plots.plot(Model[:,5], Model[:,1],yflip=true,xlim=(0,600),ylim=(-8,350),label="Qs",lc=RGB(13/255, 8/255, 135/255),xlabel="Qs",legend = :bottomleft)
    pp2 = twiny()
    Plots.plot!(pp2,Melt, zTPC[:,1],yflip=true,xlim=(-1,3),ylim=(-8,350),label="Melt (%)",lc=RGB(204/255, 71/255, 120/255),legend = :topright,xlabel="Melt (%)",grid=false)

    Plots.plot(p1,p2,p3, layout=(1,3),size=(1200,800),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
    savefig("figs/2.-Lithos_Model.pdf")

    ## Model of the Transition Mantle

    p1=Plots.plot(Model[:,2], Model[:,1],yflip=true,xlim=(2,14),ylim=(350,800),label="Vp",lc=RGB(68/255, 1/255, 84/255), title="T. Zone Model",legend=:bottomleft)
    Plots.plot!(Model2[:,2], Model2[:,1],yflip=true,ylim=(350,800),label="PREM Vp",lc=RGB(68/255, 1/255, 84/255),ls=:dot)
    Plots.plot!(Model[:,3], Model[:,1],label="Vs",lc=RGB(42/255, 120/255, 142/255))
    Plots.plot!(Model2[:,3], Model2[:,1],label="PREM Vs",lc=RGB(42/255, 120/255, 142/255),ls=:dot)
    Plots.plot!(Model[:,4], Model[:,1],label="Rho",lc=RGB(122/255, 209/255, 81/255))
    Plots.plot!(Model2[:,4], Model2[:,1],label="PREM Rho",lc=RGB(122/255, 209/255, 81/255),ls=:dot,minorgrid=true,xlabel="Vx (km/s) & Rho (g/cm^3)")

    p2=Plots.plot(zTPC[:,2], zTPC[:,1],yflip=true,xlim=(200,4000),ylim=(350,800),label="T (K)",lc=RGB(254/255, 176/255, 120/255),xlabel="T (K)",legend = :bottomleft)
    Plots.plot!(T_Sol, Z_SL,label="Solidus",lc=RGB( 85/255, 107/255,  47/255))
    Plots.plot!(T_Liq, Z_SL,label="Liquidus",lc=RGB(70/255, 130/255, 180/255))
    pp = twiny()
    Plots.plot!(pp,zTPC[:,3]*1e-4, zTPC[:,1],yflip=true,xlim=(0,35),ylim=(350,800),label="P (GPa)",lc=RGB(0/255, 0/255, 4/255),legend = :topright,xlabel="P (GPa)",grid=false)
    t1, t2 = 200, 4000;
    y = [Temp_Depth_Nodes[2]/1000]
    Plots.plot!([t1; t2], [y; y],label="LAB",lc=RGB(44/255, 17/255, 95/255),ls=:dash)
    y = [zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-14)),1]]
    Plots.plot!([t1; t2], [y; y],label="410 Disc.",lc=RGB(114/255, 31/255, 129/255),ls=:dash)
    y = [zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-24)),1]]
    Plots.plot!([t1; t2], [y; y],label="660 Disc.",lc=RGB(183/255, 55/255, 121/255),ls=:dash)
    y = mean([zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-120)),1],zTPC[argmin(abs.(zTPC[:,2] .-2500)),1]])
    Plots.plot!([t1; t2], [y; y],label="D''",lc=RGB(221/255, 81/255, 58/255),ls=:dash)

    p3=Plots.plot(Model[:,5], Model[:,1],yflip=true,xlim=(0,600),ylim=(350,800),label="Qs",lc=RGB(13/255, 8/255, 135/255),xlabel="Qs",legend = :bottomleft)
    pp2 = twiny()
    Plots.plot!(pp2,Melt, zTPC[:,1],yflip=true,xlim=(-1,3),ylim=(350,800),label="Melt (%)",lc=RGB(204/255, 71/255, 120/255),legend = :topright,xlabel="Melt (%)",grid=false)

    Plots.plot(p1,p2,p3, layout=(1,3),size=(1200,800),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
    savefig("figs/3.-T_Zone_Model.pdf")



    ## Model of the Crust
    p1=Plots.plot(Model[:,2], Model[:,1],yflip=true,xlim=(2,10),ylim=(-8,100),label="Vp",lc=RGB(68/255, 1/255, 84/255), title="Crust Model",legend=:bottomleft)
    Plots.plot!(Model2[:,2], Model2[:,1],yflip=true,ylim=(-8,100),label="PREM Vp",lc=RGB(68/255, 1/255, 84/255),ls=:dot)
    Plots.plot!(Model[:,3], Model[:,1],label="Vs",lc=RGB(42/255, 120/255, 142/255))
    Plots.plot!(Model2[:,3], Model2[:,1],label="PREM Vs",lc=RGB(42/255, 120/255, 142/255),ls=:dot)
    Plots.plot!(Model[:,4], Model[:,1],label="Rho",lc=RGB(122/255, 209/255, 81/255))
    Plots.plot!(Model2[:,4], Model2[:,1],label="PREM Rho",lc=RGB(122/255, 209/255, 81/255),ls=:dot,minorgrid=true,xlabel="Vx (km/s) & Rho (g/cm^3)")

    p2=Plots.plot(zTPC[:,2], zTPC[:,1],yflip=true,xlim=(200,1500),ylim=(-8,100),label="T (K)",lc=RGB(254/255, 176/255, 120/255),xlabel="T (K)",legend = :bottomleft)
    Plots.plot!(T_Sol, Z_SL,label="Solidus",lc=RGB( 85/255, 107/255,  47/255))
    Plots.plot!(T_Liq, Z_SL,label="Liquidus",lc=RGB(70/255, 130/255, 180/255))
    pp = twiny()
    Plots.plot!(pp,zTPC[:,3]*1e-4, zTPC[:,1],yflip=true,xlim=(0,5),ylim=(-8,100),label="P (GPa)",lc=RGB(0/255, 0/255, 4/255),legend = :topright,xlabel="P (GPa)",grid=false)
    t1, t2 = 200, 4000;
    y = [Temp_Depth_Nodes[2]/1000]
    Plots.plot!([t1; t2], [y; y],label="LAB",lc=RGB(44/255, 17/255, 95/255),ls=:dash)
    y = [zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-14)),1]]
    Plots.plot!([t1; t2], [y; y],label="410 Disc.",lc=RGB(114/255, 31/255, 129/255),ls=:dash)
    y = [zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-24)),1]]
    Plots.plot!([t1; t2], [y; y],label="660 Disc.",lc=RGB(183/255, 55/255, 121/255),ls=:dash)
    y = mean([zTPC[argmin(abs.(zTPC[:,3]*1e-4 .-120)),1],zTPC[argmin(abs.(zTPC[:,2] .-2500)),1]])
    Plots.plot!([t1; t2], [y; y],label="D''",lc=RGB(221/255, 81/255, 58/255),ls=:dash)

    p3=Plots.plot(Model[:,5], Model[:,1],yflip=true,xlim=(0,600),ylim=(-8,100),label="Qs",lc=RGB(13/255, 8/255, 135/255),xlabel="Qs",legend = :bottomleft)
    pp2 = twiny()
    Plots.plot!(pp2,Melt, zTPC[:,1],yflip=true,xlim=(-1,3),ylim=(-8,100),label="Melt (%)",lc=RGB(204/255, 71/255, 120/255),legend = :topright,xlabel="Melt (%)",grid=false)
    Plots.plot(p1,p2,p3, layout=(1,3),size=(1200,800),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
    savefig("figs/1.-Crust_Model.pdf")

    ## Only Properties



     # RF observed
     RF1=readdlm("data/Obs_Stacked_PRF.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB
     RF2=readdlm("data/Obs_Stacked_SRF.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB
     RF3=readdlm("data/Obs_Stacked_SKSRF.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB

     #RF Stacks
     ps3=Plots.plot(Stacked_PRF[:,1],Stacked_PRF[:,2]./maximum(Stacked_PRF[:,2]),lc=:black,label="P-to-S Synthetic",ylabel="Norm. A",xlabel="Time After P (s)",minorgrid=true,xlim=(0,80),ylim=(-1,1))
     Plots.plot!(RF1[:,1],RF1[:,2]./maximum(RF1[:,2]),lc=:blue,lw=1,label="Observed")
     ps2=Plots.plot(Stacked_SRF[:,1],Stacked_SRF[:,2]./maximum(Stacked_SRF[:,2]),lc=:black,label="S-to-P Synthetic",ylabel="Norm. A",xlabel="Time Before S (s)",minorgrid=true,xlim=(0,80),ylim=(-1,1))
     Plots.plot!(RF2[:,1],RF2[:,2]./maximum(RF2[:,2]),lc=:red,lw=1,label="Observed")
     ps1=Plots.plot(Stacked_SKSRF[:,1],Stacked_SKSRF[:,2]./maximum(Stacked_SKSRF[:,2]),lc=:black,label="SKS-to-P Synthetic",ylabel="Norm. A",xlabel="Time Before S (s)",minorgrid=true,xlim=(0,80),ylim=(-1,1))
     Plots.plot!(RF3[:,1],RF3[:,2]./maximum(RF3[:,2]),lc=:red,lw=1,label="Observed")

     #RF Matrix
     po1=heatmap(Obs_SKSRF[3:end,1],dist_SKS,Obs_SKSRF[3:end,2:end]',c=:balance,clims=(-maximum([maximum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]) abs(minimum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]))]),maximum([maximum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]) abs(minimum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]))])),xlim=(0,100) ,xlabel = "Time Before SKS (s)", ylabel = "Distance (°)",title="Observed RF")
     po2=heatmap(Obs_SRF[3:end,1],dist_S,Obs_SRF[3:end,2:end]',c=:balance,clims=(-maximum([maximum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]) abs(minimum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]))]),maximum([maximum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]) abs(minimum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]))])),xlim=(0,50),xlabel = "Time Before S (s)", ylabel = "Distance (°)")
     po3=heatmap(Obs_PRF[3:end,1],dist_P,Obs_PRF[3:end,2:end]',c=:balance,clims=(-maximum([maximum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]) abs(minimum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]))]),maximum([maximum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]) abs(minimum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]))])),xlim=(0,90),xlabel = "Time After P (s)", ylabel = "Distance (°)")
     pc1=heatmap(SKS_matrix[:,1],dist_SKS,SKS_matrix[:,2:end]',c=:balance,clims=(-maximum([maximum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]) abs(minimum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]))]),maximum([maximum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]) abs(minimum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]))])),xlim=(0,100),xlabel = "Time Before SKS (s)",title="Synthetic RF")
     pc2=heatmap(S_matrix[:,1],dist_S,S_matrix[:,2:end]',c=:balance,clims=(-maximum([maximum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]) abs(minimum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]))]),maximum([maximum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]) abs(minimum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]))])),xlim=(0,50),xlabel = "Time Before S (s)")
     pc3=heatmap(P_matrix[:,1],dist_P,P_matrix[:,2:end]',c=:balance,clims=(-maximum([maximum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]) abs(minimum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]))]),maximum([maximum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]) abs(minimum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]))])),xlim=(0,90),xlabel = "Time After P (s)")

     Plots.plot(ps1,po1, pc1,ps2, po2, pc2,ps3, po3, pc3, layout = (3,3),size=(2400,900),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
     savefig("figs/6.-Receiver_Functions.pdf")

    # SW misfit
    ## Dispersion Curves
    # fix the input to match the output format
      Obs_SWD_M=zeros(size(SWV));
      if Nmodes>1
          for u=1:Nmodes
              A=2+Nmodes*(u-1)+(u-1);
              B=A+Nmodes;
              Obs_SWD_M[:,:,u]=Obs_SWD[:,A:B];
          end
      else
          Obs_SWD_M=Obs_SWD[:,2:end];
      end


    ## Plot the dispersion Curves

    ## Plot the dispersion Curves
    # Rc
    p1=Plots.plot(T,SWV[:,1,:].*reverse(Obs_SWD_M[:,1,:],dims=1)./reverse(Obs_SWD_M[:,1,:],dims=1),xlim=(0,300),ylim=(2,10),legend=false,ylabel="Rayleigh c (km/s)",markerstrokewidth=0.1,minorgrid=true,title="Dispersion Curves")
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,1,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0)
    # RU
    p2=Plots.plot(T,SWV[:,2,:].*reverse(Obs_SWD_M[:,2,:],dims=1)./reverse(Obs_SWD_M[:,2,:],dims=1),xlim=(0,300),ylim=(2,7),legend=false,ylabel="Rayleigh U (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,2,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0)
    # Lc
    p3=Plots.plot(T,SWV[:,3,:].*reverse(Obs_SWD_M[:,3,:],dims=1)./reverse(Obs_SWD_M[:,3,:],dims=1),xlim=(0,300),ylim=(2,10),legend=false,ylabel="Love c (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,3,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0)
    # LU
    p4=Plots.plot(T,SWV[:,4,:].*reverse(Obs_SWD_M[:,4,:],dims=1)./reverse(Obs_SWD_M[:,4,:],dims=1),xlim=(0,300),ylim=(2,7),legend=false,xlabel="Period (s)",ylabel="Love U (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,4,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0)

    ## Plot the anomalies
    p5=Plots.scatter(T, (SWV[:,1,:].-reverse(Obs_SWD_M[:,1,:],dims=1))*100 ./reverse(Obs_SWD_M[:,1,:],dims=1) ,xlim=(0,300),legend=false,ylabel="Rayleigh dc (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,ylim=[-ceil(maximum(abs.(SWV[:,1,:].-reverse(Obs_SWD_M[:,1,:],dims=1))*100 ./reverse(Obs_SWD_M[:,1,:],dims=1))),ceil(maximum(abs.(SWV[:,1,:].-reverse(Obs_SWD_M[:,1,:],dims=1))*100 ./reverse(Obs_SWD_M[:,1,:],dims=1)))],title="Misfit" )
    p6=Plots.scatter(T, (SWV[:,2,:].-reverse(Obs_SWD_M[:,2,:],dims=1))*100 ./reverse(Obs_SWD_M[:,2,:],dims=1) ,xlim=(0,300),legend=false,ylabel="Rayleigh dU (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,ylim=[-ceil(maximum(abs.(SWV[:,2,:].-reverse(Obs_SWD_M[:,2,:],dims=1))*100 ./reverse(Obs_SWD_M[:,2,:],dims=1))),ceil(maximum(abs.(SWV[:,2,:].-reverse(Obs_SWD_M[:,2,:],dims=1))*100 ./reverse(Obs_SWD_M[:,2,:],dims=1)))] )
    p7=Plots.scatter(T, (SWV[:,3,:].-reverse(Obs_SWD_M[:,3,:],dims=1))*100 ./reverse(Obs_SWD_M[:,3,:],dims=1) ,xlim=(0,300),legend=false,ylabel="Love dc (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,ylim=[-ceil(maximum(abs.(SWV[:,3,:].-reverse(Obs_SWD_M[:,3,:],dims=1))*100 ./reverse(Obs_SWD_M[:,3,:],dims=1))),ceil(maximum(abs.(SWV[:,3,:].-reverse(Obs_SWD_M[:,3,:],dims=1))*100 ./reverse(Obs_SWD_M[:,3,:],dims=1)))] )
    p8=Plots.scatter(T, (SWV[:,4,:].-reverse(Obs_SWD_M[:,4,:],dims=1))*100 ./reverse(Obs_SWD_M[:,4,:],dims=1) ,xlim=(0,300),legend=false,ylabel="Love dU (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,xlabel="Period (s)",ylim=[-ceil(maximum(abs.(SWV[:,4,:].-reverse(Obs_SWD_M[:,4,:],dims=1))*100 ./reverse(Obs_SWD_M[:,4,:],dims=1))),ceil(maximum(abs.(SWV[:,4,:].-reverse(Obs_SWD_M[:,4,:],dims=1))*100 ./reverse(Obs_SWD_M[:,4,:],dims=1)))] )
    # Put it all together
    Plots.plot(p1,p5,p2,p6,p3,p7,p4,p8, layout=(4,2),size=(900,1200),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
    savefig("figs/5a.-Surface_Wave_Dispersion.pdf")


    ## Plot the dispersion Curves LOG
    # Rc
    p1=Plots.plot(T,SWV[:,1,:].*reverse(Obs_SWD_M[:,1,:],dims=1)./reverse(Obs_SWD_M[:,1,:],dims=1),xlim=(20,300),ylim=(2,10),legend=false,ylabel="Rayleigh c (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,1,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0,xaxis=(:log, [10, :auto]))
    # RU
    p2=Plots.plot(T,SWV[:,2,:].*reverse(Obs_SWD_M[:,2,:],dims=1)./reverse(Obs_SWD_M[:,2,:],dims=1),xlim=(20,300),ylim=(2,7),legend=false,ylabel="Rayleigh U (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,2,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0,xaxis=(:log, [10, :auto]))
    # Lc
    p3=Plots.plot(T,SWV[:,3,:].*reverse(Obs_SWD_M[:,3,:],dims=1)./reverse(Obs_SWD_M[:,3,:],dims=1),xlim=(20,300),ylim=(2,10),legend=false,ylabel="Love c (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,3,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0,xaxis=(:log, [10, :auto]))
    # LU
    p4=Plots.plot(T,SWV[:,4,:].*reverse(Obs_SWD_M[:,4,:],dims=1)./reverse(Obs_SWD_M[:,4,:],dims=1),xlim=(20,300),ylim=(2,7),legend=false,xlabel="Period (s)",ylabel="Love U (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,4,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0,xaxis=(:log, [10, :auto]))

    ## Plot the anomalies
    p5=Plots.scatter(T, (SWV[:,1,:].-reverse(Obs_SWD_M[:,1,:],dims=1))*100 ./reverse(Obs_SWD_M[:,1,:],dims=1) ,xlim=(20,300),legend=false,ylabel="Rayleigh dc (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,xaxis=(:log, [10, :auto]),ylim=[-ceil(maximum(abs.(SWV[:,1,:].-reverse(Obs_SWD_M[:,1,:],dims=1))*100 ./reverse(Obs_SWD_M[:,1,:],dims=1))),ceil(maximum(abs.(SWV[:,1,:].-reverse(Obs_SWD_M[:,1,:],dims=1))*100 ./reverse(Obs_SWD_M[:,1,:],dims=1)))] )
    p6=Plots.scatter(T, (SWV[:,2,:].-reverse(Obs_SWD_M[:,2,:],dims=1))*100 ./reverse(Obs_SWD_M[:,2,:],dims=1) ,xlim=(20,300),legend=false,ylabel="Rayleigh dU (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,xaxis=(:log, [10, :auto]),ylim=[-ceil(maximum(abs.(SWV[:,2,:].-reverse(Obs_SWD_M[:,2,:],dims=1))*100 ./reverse(Obs_SWD_M[:,2,:],dims=1))),ceil(maximum(abs.(SWV[:,2,:].-reverse(Obs_SWD_M[:,2,:],dims=1))*100 ./reverse(Obs_SWD_M[:,2,:],dims=1)))] )
    p7=Plots.scatter(T, (SWV[:,3,:].-reverse(Obs_SWD_M[:,3,:],dims=1))*100 ./reverse(Obs_SWD_M[:,3,:],dims=1) ,xlim=(20,300),legend=false,ylabel="Love dc (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,xaxis=(:log, [10, :auto]),ylim=[-ceil(maximum(abs.(SWV[:,3,:].-reverse(Obs_SWD_M[:,3,:],dims=1))*100 ./reverse(Obs_SWD_M[:,3,:],dims=1))),ceil(maximum(abs.(SWV[:,3,:].-reverse(Obs_SWD_M[:,3,:],dims=1))*100 ./reverse(Obs_SWD_M[:,3,:],dims=1)))] )
    p8=Plots.scatter(T, (SWV[:,4,:].-reverse(Obs_SWD_M[:,4,:],dims=1))*100 ./reverse(Obs_SWD_M[:,4,:],dims=1) ,xlim=(20,300),legend=false,ylabel="Love dU (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,xaxis=(:log, [10, :auto]),xlabel="Period (s)",ylim=[-ceil(maximum(abs.(SWV[:,4,:].-reverse(Obs_SWD_M[:,4,:],dims=1))*100 ./reverse(Obs_SWD_M[:,4,:],dims=1))),ceil(maximum(abs.(SWV[:,4,:].-reverse(Obs_SWD_M[:,4,:],dims=1))*100 ./reverse(Obs_SWD_M[:,4,:],dims=1)))] )
    # Put it all together
    Plots.plot(p1,p5,p2,p6,p3,p7,p4,p8, layout=(4,2),size=(900,1200),bottom_margin = 1Plots.cm, left_margin=1Plots.cm,right_margin=1Plots.cm,top_margin=1Plots.cm)
    savefig("figs/5b.-Surface_Wave_Dispersion_log.pdf")

    # THE END
    println(" ")
    printstyled(" C'est Fini",color=:green)
    println(" ")
     #
return
end


## A funtion to put to X axis
function twiny(sp::Plots.Subplot)
    sp[:top_margin] = max(sp[:top_margin], 30Plots.px)
    Plots.plot!(sp.plt, inset = (sp[:subplot_index], bbox(0,0,1,1)))
    twinsp = sp.plt.subplots[end]
    twinsp[:xaxis][:mirror] = true
    twinsp[:background_color_inside] = RGBA{Float64}(0,0,0,0)
    Plots.link_axes!(sp[:yaxis], twinsp[:yaxis])
    twinsp
end
twiny(plt::Plots.Plot = current()) = twiny(plt[1])

## A function to get the Solidus and Liquidus curves to plot!
function Make_Sol_Liq(zTPC)

indx_m=findall(x->x=="M", zTPC[:,8])
# P=zTPC[indx_m[1]:argmin(abs.(zTPC[:,3]*1e-4 .- 14)),3].*1e-4; # pressure in the nodes in GPa
# T=zTPC[indx_m[1]:argmin(abs.(zTPC[:,3]*1e-4 .- 14)),2]; # T in K
# z=zTPC[indx_m[1]:argmin(abs.(zTPC[:,3]*1e-4 .- 14)),1]; # Z in km
# #
P=zTPC[indx_m,3]*1e-4 # pressure in the nodes in GPa
T=zTPC[indx_m,2]; # T in K
z=zTPC[indx_m,1]; # Z in km
##
T_Solidus_Matrix=zeros(size(P,1),1); #because 6 cases
T_Liquidus_Matrix=zeros(size(P,1),); #because 6 cases

## Model Made From Katz, Gerya Latisov, and Fu
for i=1:size(P,1)
	# 4: Lithospheric and Asthenospheric Mantle. With latent heat 400 kJ/kg

	# In the Upper Mantle: Katz et al 2003
    if P[i]< 10
			 # Solidus Temperature
        	 T_Solidus_Matrix[i,1]=(1085.7+132.9*P[i]-5.1*P[i]^2)+273.15;
			 # Liquidus temperature
	 		 T_Liquidus_Matrix[i,1]=1968 + 67.1*P[i] + -4.35*P[i]^2 + 0.17*P[i]^3 + -2.25e-03*P[i]^4
	# Correct for higher presures P > 10 as Gerya does, to avoid comeback of the function
	elseif  P[i] >= 10.0 && P[i] <=25
			 # Solidus Temperature
        	 T_Solidus_Matrix[i,1]=(1086-5.7*P[i]+390*log(P[i]))+273.15#2178+0.030819*((P[i]*1000)-10000);
			 # Liquidus temperature
			 #T_Liquidus_Matrix[i,1]=(1780+45*P[i]-2*P[i]^2)+273.15
	 		 T_Liquidus_Matrix[i,1]=1968 + 67.1*P[i] + -4.35*P[i]^2 + 0.17*P[i]^3 + -2.25e-03*P[i]^4
			 # Correct for higher presures P > 10 as Gerya does, to avoid comeback of the function
    # In the Lower Mantle: Fu et al 2018
	elseif P[i] >= 25.0
		# Solidus Temperature
		T_Solidus_Matrix[i,1]=2328+14.6*P[i]-0.0412*P[i]^2;
		# Liquidus temperature
		T_Liquidus_Matrix[i,1]=2876+15.7*P[i]-0.0421*P[i]^2;

	end
end

return z,T_Solidus_Matrix, T_Liquidus_Matrix
end
