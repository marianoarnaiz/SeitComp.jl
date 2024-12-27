"This is a function Helps set up all the enviroment. JUST RUN ONCE!"

using DelimitedFiles, GMT, Plots

function Draw_Plots()
    println(" ")
    printstyled(" Plotting Figures. Go to Mantle_Mod_1D_v0.0/figs",color=:green)
    println(" ")


    #reload models just in case
    Model2=readdlm("data/AK135.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB
    ## Plot the model and all the properties

    LAB=argmin(abs.(zTPC[:,2].-Temperature_Nodes[2]));
    Cut=argmin(abs.(zTPC[:,1].-300));

    ## Model of the Mantle
    gmt("psxy -R-8/2890/2/16 -W1,blue -JX16/5.25 -Ax  -Bxa500f100 -Bya2f1+lVp(b)_Vs(r)_Rho(g) -K -P  > figs/Model_Mantle.ps", [Model[:,1] Model[:,2]])
    gmt("psxy -R-8/2890/2/16 -W1,blue,..  -JX16/5.25 -Ax  -K -O  >> figs/Model_Mantle.ps", [Model2[:,1] Model2[:,2]])
    gmt("psxy -R-8/2890/2/16 -W1,red -JX16/5.25 -Ax -K -O  >> figs/Model_Mantle.ps", [Model[:,1] Model[:,3]])
    gmt("psxy -R-8/2890/2/16 -W1,red,..  -JX16/5.25 -Ax -K -O  >> figs/Model_Mantle.ps", [Model2[:,1] Model2[:,3]])
    gmt("psxy -R-8/2890/2/16 -W1,darkgreen -JX16/5.25 -Ax -K -O >> figs/Model_Mantle.ps", [Model[:,1] Model[:,4]])
    gmt("psxy -R-8/2890/2/16 -W1,darkgreen,..  -JX16/5.25 -Ax -K -O >> figs/Model_Mantle.ps", [Model2[:,1] Model2[:,4]])
    gmt("psxy -R-8/2890/200/4000 -W1,red -JX16/5.25 -Bxa500f100 -Bya500f100+lT_[k] -Ax -K -O -Y7 >> figs/Model_Mantle.ps", [zTPC[:,1] zTPC[:,2]])
    gmt("psxy -R-8/2890/0/150 -W1,springgreen -JX16/5.25 -Bxa500f100 -Bya30f10+lP_[GPa] -Ax -K -O -Y7 >> figs/Model_Mantle.ps", Float64.([zTPC[:,1] zTPC[:,3]*0.0001]))
    gmt("psxy -R-8/2890/0/10 -W1,maroon4 -JX16/5.25 -Bxa500f100 -Bya1f0.5+lComposition -Ax -K -O -Y7 >> figs/Model_Mantle.ps", [zTPC[LAB:end,1] zTPC[LAB:end,4]])
    gmt("psxy -R-8/2890/0/10 -W1,black -JX16/5.25 -Ax -O >> figs/Model_Mantle.ps", [zTPC[LAB:end,1] zTPC[LAB:end,5]])
    gmt( "psconvert figs/Model_Mantle.ps -A -E300 -Tg")

    ## Model of the Lithosphere
    gmt("psxy -R-8/300/2/10 -W1,blue -JX16/5.25 -Ax  -Bxa50f25 -Bya2f1+lVp(b)_Vs(r)_Rho(g) -K -P  > figs/Model_Lithos.ps", [Model[:,1] Model[:,2]])
    gmt("psxy -R-8/300/2/10 -W1,blue,..  -JX16/5.25 -Ax  -K -O  >> figs/Model_Lithos.ps", [Model2[:,1] Model2[:,2]])
    gmt("psxy -R-8/300/2/10 -W1,red -JX16/5.25 -Ax -K -O  >> figs/Model_Lithos.ps", [Model[:,1] Model[:,3]])
    gmt("psxy -R-8/300/2/10 -W1,red,..  -JX16/5.25 -Ax -K -O  >> figs/Model_Lithos.ps", [Model2[:,1] Model2[:,3]])
    gmt("psxy -R-8/300/2/10 -W1,darkgreen -JX16/5.25 -Ax -K -O >> figs/Model_Lithos.ps", [Model[:,1] Model[:,4]])
    gmt("psxy -R-8/300/2/10 -W1,darkgreen,..  -JX16/5.25 -Ax -K -O >> figs/Model_Lithos.ps", [Model2[:,1] Model2[:,4]])
    gmt("psxy -R-8/300/200/2000 -W1,red -JX16/5.25 -Bxa50f25 -Bya500f100+lT_[k] -K -O -Y7 >> figs/Model_Lithos.ps", [zTPC[:,1] zTPC[:,2]])
    gmt("psxy -R-8/300/0/15 -W1,springgreen -JX16/5.25 -Bxa50f25 -Bya5f1+lP_[GPa] -O -Y7 >> figs/Model_Lithos.ps", Float64.([zTPC[:,1] zTPC[:,3]*0.0001]))
    #gmt("psxy -R-8/300/0/100 -W1,maroon4 -JX16/5.25 -Ba50f25/a20f5:Minerals(%): -Ax -K -O -Y7 >> figs/Model_Lithos.ps", Float64.([zTPC[1:LAB,1] zTPC[1:LAB,4]]))
    #gmt("psxy -R-8/300/0/100 -W1,black -JX16/5.25 -Ax -K -O >> figs/Model_Lithos.ps", Float64.([zTPC[1:LAB,1] zTPC[1:LAB,5]]))
    #gmt("psxy -R-8/300/0/100 -W1,darkgreen -JX16/5.25 -Ax -K -O >> figs/Model_Lithos.ps", Float64.([zTPC[1:LAB,1] zTPC[1:LAB,6]]))
    #gmt("psxy -R-8/300/0/100 -W1,darkblue -JX16/5.25 -Ax -K -O >> figs/Model_Lithos.ps", Float64.([zTPC[1:LAB,1] zTPC[1:LAB,7]]))
    gmt( "psconvert figs/Model_Lithos.ps -A -E300 -Tg")

    ## Model of the Crust

    gmt("psxy -R-8/80/2/10 -W1,blue -JX16/5.25 -Ax  -Bxa10f5 -Bya2f1+lVp(b)_Vs(r)_Rho(g) -K -P  > figs/Model_Crust.ps", [Model[:,1] Model[:,2]])
    gmt("psxy -R-8/80/2/10 -W1,blue,..  -JX16/5.25 -Ax  -K -O  >> figs/Model_Crust.ps", [Model2[:,1] Model2[:,2]])
    gmt("psxy -R-8/80/2/10 -W1,red -JX16/5.25 -Ax -K -O  >> figs/Model_Crust.ps", [Model[:,1] Model[:,3]])
    gmt("psxy -R-8/80/2/10 -W1,red,..  -JX16/5.25 -Ax -K -O  >> figs/Model_Crust.ps", [Model2[:,1] Model2[:,3]])
    gmt("psxy -R-8/80/2/10 -W1,darkgreen -JX16/5.25 -Ax -K -O >> figs/Model_Crust.ps", [Model[:,1] Model[:,4]])
    gmt("psxy -R-8/80/2/10 -W1,darkgreen,..  -JX16/5.25 -Ax -K -O >> figs/Model_Crust.ps", [Model2[:,1] Model2[:,4]])
    gmt("psxy -R-8/80/200/1500 -W1,red -JX16/5.25 -Bxa10f5 -Bya500f100+lT_[k]  -K -O -Y7 >> figs/Model_Crust.ps", [zTPC[:,1] zTPC[:,2]])
    gmt("psxy -R-8/80/0/5 -W1,springgreen -JX16/5.25  -Bxa5f1 -Bya5f1+lP_[GPa] -O -Y7 >> figs/Model_Crust.ps", Float64.([zTPC[:,1] zTPC[:,3]*0.0001]))
    #gmt("psxy -R-8/80/0/100 -W1,maroon4 -JX16/5.25 -Ba10f5/a20f5:Minerals(%): -Ax -K -O -Y7 >> figs/Model_Crust.ps", Float64.([zTPC[1:LAB,1] zTPC[1:LAB,4]]))
    #gmt("psxy -R-8/80/0/100 -W1,black -JX16/5.25 -Ax -K -O >> figs/Model_Crust.ps", Float64.([zTPC[1:LAB,1] zTPC[1:LAB,5]]))
    #gmt("psxy -R-8/80/0/100 -W1,darkgreen -JX16/5.25 -Ax -K -O >> figs/Model_Crust.ps", Float64.([zTPC[1:LAB,1] zTPC[1:LAB,6]]))
    #gmt("psxy -R-8/300/0/100 -W1,darkblue -JX16/5.25 -Ax -K -O >> figs/Model_Lithos.ps", Float64.([zTPC[1:LAB,1] zTPC[1:LAB,7]]))
    gmt( "psconvert figs/Model_Crust.ps -A -E300 -Tg")

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
      # R wave c
      gmt("psxy -R-0/300/2/10 -W0.5,red -JX14/5  -Bxa50f10 -Bya1f1+lR_wave_c_(km/s)  -K -P -X3.5  > figs/SWF_Curves.ps", [T SWV[:,1,1]])
      gmt("psxy -R-0/300/2/10 -Sc0.1 -W0.5 -JX14/5 -K -O  >> figs/SWF_Curves.ps", [Obs_SWD[:,1] (Obs_SWD_M[:,1,1])])
      for i=2:Nmodes
          gmt("psxy -R-0/300/2/10 -W0.5,red -JX14/5 -K -O  >> figs/SWF_Curves.ps", [T SWV[:,1,i]])
          gmt("psxy -R-0/300/2/10 -Sc0.1 -W0.5 -JX14/5 -K -O  >> figs/SWF_Curves.ps", [Obs_SWD[:,1] (Obs_SWD_M[:,1,i])])
      end

      # R wave U
      gmt("psxy -R-0/300/2/7 -W0.5,red -JX14/5 -Bxa50f10 -Bya1f1+lR_wave_U_(km/s)  -K -O -Y7  >> figs/SWF_Curves.ps", [T SWV[:,2,1]])
      gmt("psxy -R-0/300/2/7 -Sc0.1 -W0.5 -JX14/5 -K -O  >> figs/SWF_Curves.ps", [Obs_SWD[:,1] (Obs_SWD_M[:,2,1])])
      for i=2:Nmodes
          gmt("psxy -R-0/300/2/7 -W0.5,red -JX14/5 -K -O  >> figs/SWF_Curves.ps", [T SWV[:,2,i]])
          gmt("psxy -R-0/300/2/7 -Sc0.1  -W0.5 -JX14/5 -K -O  >> figs/SWF_Curves.ps", [Obs_SWD[:,1] (Obs_SWD_M[:,2,i])])
      end
      #gmt("psxy -R-0/200/2/7 -W0.5,black -JX14/5 -K -O  >> figs/SWF_Curves.ps", [T SWV[:,2,3]])

      # L wave c
      gmt("psxy -R-0/300/2/10 -W0.5,red -JX14/5 -Bxa50f10 -Bya1f1+lL_wave_c_(km/s)  -K -O -Y7  >> figs/SWF_Curves.ps", [T SWV[:,3,1]])
      gmt("psxy -R-0/300/2/10 -Sc0.1 -W0.5 -JX14/5 -K -O  >> figs/SWF_Curves.ps", [Obs_SWD[:,1] (Obs_SWD_M[:,3,1])])
      for i=2:Nmodes
          gmt("psxy -R-0/300/2/10 -W0.5,red -JX14/5 -K -O  >> figs/SWF_Curves.ps", [T SWV[:,3,i]])
          gmt("psxy -R-0/300/2/10 -Sc0.1 -W0.5 -JX14/5 -K -O  >> figs/SWF_Curves.ps", [Obs_SWD[:,1] (Obs_SWD_M[:,3,i])])
      end
      #gmt("psxy -R-0/200/2/10 -W0.5,black -JX14/5 -K -O  >> figs/SWF_Curves.ps", [T SWV[:,3,3]])

      # L wave U
      gmt("psxy -R-0/300/2/7 -W0.5,red -JX14/5 -Bxa50f10 -Bya1f1+lL_wave_U_(km/s) -K -O -Y7  >> figs/SWF_Curves.ps", [T SWV[:,4,1]])
      gmt("psxy -R-0/300/2/7 -Sc0.1 -W0.5 -JX14/5 -K -O  >> figs/SWF_Curves.ps", [Obs_SWD[:,1] (Obs_SWD_M[:,4,1])])
      for i=2:Nmodes-1
          gmt("psxy -R-0/300/2/7 -W0.5,red -JX14/5 -K -O  >> figs/SWF_Curves.ps", [T SWV[:,4,i]])
          gmt("psxy -R-0/300/2/7 -Sc0.1 -W0.5  -JX14/5 -K -O  >> figs/SWF_Curves.ps", [Obs_SWD[:,1] (Obs_SWD_M[:,4,i])])
      end
      gmt("psxy -R-0/300/2/7 -W0.5,red -JX14/5 -K -O  >> figs/SWF_Curves.ps", [T SWV[:,4,Nmodes]])
      gmt("psxy -R-0/300/2/7 -Sc0.1 -W0.5  -JX14/5 -O  >> figs/SWF_Curves.ps", [Obs_SWD[:,1] (Obs_SWD_M[:,4,Nmodes])])
      gmt( "psconvert figs/SWF_Curves.ps -A -E300 -Tg")

      #gmt("psxy -R-0/200/2/7 -W0.5,black -JX14/5 -K -O  >> figs/SWF_Curves.ps", [T SWV[:,4,3]])

      ## RF sections
      # # P,S and SKS RF Sections (all 3 in the same Figure!)
       gmt("psxy -R-0/90/25/95 -W1,red  -JX14/5.5 -Bxa10f5+lt_after_P_(s) -Bya10f5+lDistance -K -P -X3.5  > figs/RF_Sections.ps", [P_matrix[:,1] 40*P_matrix[:,2].+dist_P[1]])
       gmt("psxy -R-0/90/25/95 -W0.5,black  -JX14/5.5 -K -O  >> figs/RF_Sections.ps", [Obs_PRF[3:end,1] 40*Obs_PRF[3:end,2].+dist_P[1]])
       for i=2:size(dist_P,1)
           gmt("psxy -R-0/90/25/95 -W1,red -JX14/5.5 -K -O  >> figs/RF_Sections.ps", [P_matrix[:,1] 40*P_matrix[:,i+1].+dist_P[i]])
           gmt("psxy -R-0/90/25/95 -W0.5,black -JX14/5.5 -K -O  >> figs/RF_Sections.ps", [Obs_PRF[3:end,1] 40*Obs_PRF[3:end,i+1].+dist_P[i]])
       end

       gmt("psxy -R-0/90/50/95 -W1,red -JX14/5.5 -Bxa10f5+lt_before_S_(s) -Bya10f5+lDistance -K -O -Y9.5  >> figs/RF_Sections.ps", [S_matrix[:,1] 40*S_matrix[:,2].+dist_S[1]])
       gmt("psxy -R-0/90/50/95 -W0.5,black -JX14/5.5  -K -O >> figs/RF_Sections.ps", [Obs_SRF[3:end,1] 40*Obs_SRF[3:end,2].+dist_S[1]])
       for i=2:size(dist_S,1)
           gmt("psxy -R-0/90/50/95 -W1,red -JX14/5.5  -K -O  >> figs/RF_Sections.ps", [S_matrix[:,1] 40*S_matrix[:,i+1].+dist_S[i]])
           gmt("psxy -R-0/90/50/95 -W0.5,black -JX14/5.5  -K -O  >> figs/RF_Sections.ps", [Obs_SRF[3:end,1] 40*Obs_SRF[3:end,i+1].+dist_S[i]])

       end

       gmt("psxy -R-0/90/80/125 -W1,red -JX14/5.5 -Bxa10f5+lt_before_SKS_(s) -Bya10f5+lDistance -K -O -Y9.5  >> figs/RF_Sections.ps", [SKS_matrix[:,1] 40*SKS_matrix[:,2].+dist_SKS[1]])
       gmt("psxy -R-0/90/80/125 -W0.5,black -JX14/5.5  -K -O  >> figs/RF_Sections.ps", [Obs_SKSRF[3:end,1] 40*Obs_SKSRF[3:end,2].+dist_SKS[1]])

       for i=2:size(dist_SKS,1)-1
           gmt("psxy -R-0/90/80/125 -W1,red -JX14/5.5  -K -O  >> figs/RF_Sections.ps", [SKS_matrix[:,1] 40*SKS_matrix[:,i+1].+dist_SKS[i]])
           gmt("psxy -R-0/90/80/125 -W0.5,black -JX14/5.5  -K -O  >> figs/RF_Sections.ps", [Obs_SKSRF[3:end,1] 40*Obs_SKSRF[3:end,i+1].+dist_SKS[i]])
       end
       gmt("psxy -R-0/90/80/125 -W0.5,black -JX14/5.5  -O  >> figs/RF_Sections.ps", [Obs_SKSRF[3:end,1] 40*Obs_SKSRF[3:end,size(dist_SKS,1)+1].+dist_SKS[size(dist_SKS,1)]])
       gmt( "psconvert figs/RF_Sections.ps -A -E300 -Tg")
     #RF matrix
     po1=heatmap(Obs_SKSRF[3:end,1],dist_SKS,Obs_SKSRF[3:end,2:end]',c=:balance,clims=(-maximum([maximum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]) abs(minimum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]))]),maximum([maximum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]) abs(minimum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]))])),xlim=(0,100) ,xlabel = "Time Before SKS (s)", ylabel = "Distance (°)",fontsize=5)
     po2=heatmap(Obs_SRF[3:end,1],dist_S,Obs_SRF[3:end,2:end]',c=:balance,clims=(-maximum([maximum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]) abs(minimum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]))]),maximum([maximum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]) abs(minimum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]))])),xlim=(0,50),xlabel = "Time Before S (s)", ylabel = "Distance (°)",fontsize=5)
     po3=heatmap(Obs_PRF[3:end,1],dist_P,Obs_PRF[3:end,2:end]',c=:balance,clims=(-maximum([maximum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]) abs(minimum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]))]),maximum([maximum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]) abs(minimum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]))])),xlim=(0,90),xlabel = "Time After P (s)", ylabel = "Distance (°)",fontsize=5)
     pc1=heatmap(SKS_matrix[:,1],dist_SKS,SKS_matrix[:,2:end]',c=:balance,clims=(-maximum([maximum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]) abs(minimum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]))]),maximum([maximum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]) abs(minimum([Obs_SKSRF[3:end,2:end] SKS_matrix[:,2:end]]))])),xlim=(0,100),xlabel = "Time Before SKS (s)",fontsize=5)
     pc2=heatmap(S_matrix[:,1],dist_S,S_matrix[:,2:end]',c=:balance,clims=(-maximum([maximum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]) abs(minimum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]))]),maximum([maximum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]) abs(minimum([Obs_SRF[3:end,2:end] S_matrix[:,2:end]]))])),xlim=(0,50),xlabel = "Time Before S (s)",fontsize=5)
     pc3=heatmap(P_matrix[:,1],dist_P,P_matrix[:,2:end]',c=:balance,clims=(-maximum([maximum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]) abs(minimum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]))]),maximum([maximum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]) abs(minimum([Obs_PRF[3:end,2:end] P_matrix[:,2:end]]))])),xlim=(0,90),xlabel = "Time After P (s)",fontsize=5)

     Plots.plot(po1, pc1, po2, pc2, po3, pc3, layout = (3,2),size=(1200,600),bottom_margin = 0.5Plots.cm, left_margin=0.5Plots.cm,right_margin=0.25Plots.cm,top_margin=0.25Plots.cm)
     savefig("figs/RF_Matrix.pdf")

     ## Stacked RF

     RF1=readdlm("data/Staked_PRF_ULN.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB
     RF2=readdlm("data/Staked_PRF_TLY.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB

     gmt("psxy -R0/80/-2/2 -W1,black  -JX14/4.0 -Bxa10f5+lt_after_P_(s) -Bya1f1 -K -P -X3.5  > figs/Stacked_RFs.ps", [Stacked_PRF[:,1] Stacked_PRF[:,2]./maximum(Stacked_PRF[:,2])])
     gmt("psxy -R0/80/-2/2 -W0.5,red,--  -JX14/4.0 -Bxa10f5+lt_after_P_(s) -Bya1f1 -O -K >> figs/Stacked_RFs.ps", [RF1[:,1] RF1[:,2]./maximum(RF1[:,2])])
     gmt("psxy -R0/80/-2/2 -W0.5,blue,--  -JX14/4.0 -Bxa10f5+lt_after_P_(s) -Bya1f1 -O -K  >> figs/Stacked_RFs.ps", [RF2[:,1] RF2[:,2]./maximum(RF2[:,2])])
     gmt("psxy -R0/80/-2/2 -W1,black  -JX14/4.0 -Bxa10f5+lt_before_S_(s) -Bya1f1 -K -P -O -Y7.5  >> figs/Stacked_RFs.ps", [Stacked_SRF[:,1] Stacked_SRF[:,2]./maximum(Stacked_SRF[:,2])])
     gmt("psxy -R0/80/-2/2 -W1,black  -JX14/4.0 -Bxa10f5+lt_before_SKS_(s) -Bya1f1 -P -O -Y7.5  >> figs/Stacked_RFs.ps", [Stacked_SKSRF[:,1] Stacked_SKSRF[:,2]./maximum(Stacked_SKSRF[:,2])])
     gmt( "psconvert figs/Stacked_RFs.ps -A -E300 -Tg")
     ## Composition plots

    try
         ternary(Sediments_Composition[:,1:3], frame=(annot=20, ticks=10, grid=:a, alabel="Quartz", blabel="Carbonates", clabel="Clays", suffix=" %"), marker=:circle, fill=:gold ,markeredgecolor=0,fmt=:pdf)
    catch e1
        ternary([-1 -1 -1], frame=(annot=20, ticks=10, grid=:a, alabel="Quartz", blabel="Carbonates", clabel="Clays", suffix=" %"), marker=:circle, fill=:white,fmt=:pdf, )
    end

    try
         ternary!(Crust_Composition[findall( x -> x == "Felsic" , Crust_Composition[:,4]),1:3],frame=(annot=20, ticks=10, grid=:a, alabel="Quartz", blabel="K_Feldespar", clabel="Plagioclase", suffix=" %"), marker=:circle,fill=:pink,markeredgecolor=0, Y=17.5,fmt=:pdf)
    catch e1
         ternary!([-1 -1 -1], frame=(annot=20, ticks=10, grid=:a, alabel="Quartz", blabel="K_Feldespar", clabel="Plagioclase", suffix=" %"), marker=:circle, fmt=:pdf,fill=:white )
    end

    try
        ternary!(Crust_Composition[findall( x -> x == "Mafic" , Crust_Composition[:,4]),1:3], frame=(annot=20, ticks=10, grid=:a, alabel="Anorthosite", blabel="Clinopyroxene", clabel="Ortopyroxene", suffix=" %"), marker=:circle,fill=:grey,markeredgecolor=0, X=17.5,fmt=:pdf)
    catch e1
        ternary!([-1 -1 -1], frame=(annot=20, ticks=10, grid=:a, alabel="Anorthosite", blabel="Clinopyroxene", clabel="Ortopyroxene", suffix=" %"), marker=:circle,fmt=:pdf, fill=:white )
    end

    try
        ternary!(Crust_Composition[findall( x -> x == "Ultramafic" , Crust_Composition[:,4]),1:3], frame=(annot=20, ticks=10, grid=:a, alabel="Olivine", blabel="Clinopyroxene", clabel="Ortopyroxene", suffix=" %"), marker=:circle,fill=:darkgreen,markeredgecolor=0, Y=-17.5,fmt=:pdf,savefig="figs/Composition.pdf" )
    catch e1
        ternary!([-1 -1 -1], frame=(annot=20, ticks=10, grid=:a, alabel="Olivine", blabel="Clinopyroxene", clabel="Ortopyroxene", suffix=" %"), marker=:circle, fill=:white, Y=-17.5, fmt=:pdf,savefig="figs/Composition.pdf")
    end

    # SW misfit


    ## Plot the dispersion Curves

    ## Plot the dispersion Curves
    # Rc
    p1=Plots.scatter(T,SWV[:,1,:].*reverse(Obs_SWD_M[:,1,:],dims=1)./reverse(Obs_SWD_M[:,1,:],dims=1),xlim=(0,300),ylim=(2,10),legend=false,ylabel="Rayleigh c (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,1,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0,xaxis=(:log, [1, :auto]))
    # RU
    p2=Plots.scatter(T,SWV[:,2,:].*reverse(Obs_SWD_M[:,2,:],dims=1)./reverse(Obs_SWD_M[:,2,:],dims=1),xlim=(0,300),ylim=(2,7),legend=false,ylabel="Rayleigh U (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,2,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0,xaxis=(:log, [1, :auto]))
    # Lc
    p3=Plots.scatter(T,SWV[:,3,:].*reverse(Obs_SWD_M[:,3,:],dims=1)./reverse(Obs_SWD_M[:,3,:],dims=1),xlim=(0,300),ylim=(2,10),legend=false,ylabel="Love c (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,3,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0,xaxis=(:log, [1, :auto]))
    # LU
    p4=Plots.scatter(T,SWV[:,4,:].*reverse(Obs_SWD_M[:,4,:],dims=1)./reverse(Obs_SWD_M[:,4,:],dims=1),xlim=(0,300),ylim=(2,7),legend=false,xlabel="Period (s)",ylabel="Love U (km/s)",markerstrokewidth=0.1,minorgrid=true)
    Plots.scatter!(Obs_SWD[:,1],Obs_SWD_M[:,4,:],xlim=(0,300),ylim=(2,10),legend=false,grid=true,minorgrid=true,markerstrokewidth=0.25,markersize=5,color=:white,strokecolor=:black,markeralpha=0.25,markerstrokealpha=1.0,xaxis=(:log, [1, :auto]))

    ## Plot the anomalies
    p5=Plots.scatter(T, (SWV[:,1,:].-reverse(Obs_SWD_M[:,1,:],dims=1))*100 ./reverse(Obs_SWD_M[:,1,:],dims=1) ,xlim=(0,300),legend=false,ylabel="Rayleigh dc (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,xaxis=(:log, [1, :auto]),ylim=[-ceil(maximum(abs.(SWV[:,1,:].-reverse(Obs_SWD_M[:,1,:],dims=1))*100 ./reverse(Obs_SWD_M[:,1,:],dims=1))),ceil(maximum(abs.(SWV[:,1,:].-reverse(Obs_SWD_M[:,1,:],dims=1))*100 ./reverse(Obs_SWD_M[:,1,:],dims=1)))] )
    p6=Plots.scatter(T, (SWV[:,2,:].-reverse(Obs_SWD_M[:,2,:],dims=1))*100 ./reverse(Obs_SWD_M[:,2,:],dims=1) ,xlim=(0,300),legend=false,ylabel="Rayleigh dU (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,xaxis=(:log, [1, :auto]),ylim=[-ceil(maximum(abs.(SWV[:,2,:].-reverse(Obs_SWD_M[:,2,:],dims=1))*100 ./reverse(Obs_SWD_M[:,2,:],dims=1))),ceil(maximum(abs.(SWV[:,2,:].-reverse(Obs_SWD_M[:,2,:],dims=1))*100 ./reverse(Obs_SWD_M[:,2,:],dims=1)))] )
    p7=Plots.scatter(T, (SWV[:,3,:].-reverse(Obs_SWD_M[:,3,:],dims=1))*100 ./reverse(Obs_SWD_M[:,3,:],dims=1) ,xlim=(0,300),legend=false,ylabel="Love dc (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,xaxis=(:log, [1, :auto]),ylim=[-ceil(maximum(abs.(SWV[:,3,:].-reverse(Obs_SWD_M[:,3,:],dims=1))*100 ./reverse(Obs_SWD_M[:,3,:],dims=1))),ceil(maximum(abs.(SWV[:,3,:].-reverse(Obs_SWD_M[:,3,:],dims=1))*100 ./reverse(Obs_SWD_M[:,3,:],dims=1)))] )
    p8=Plots.scatter(T, (SWV[:,4,:].-reverse(Obs_SWD_M[:,4,:],dims=1))*100 ./reverse(Obs_SWD_M[:,4,:],dims=1) ,xlim=(0,300),legend=false,ylabel="Love dU (%)",grid=true,minorgrid=true,markerstrokewidth=0.1,xaxis=(:log, [1, :auto]),xlabel="Period (s)",ylim=[-ceil(maximum(abs.(SWV[:,4,:].-reverse(Obs_SWD_M[:,4,:],dims=1))*100 ./reverse(Obs_SWD_M[:,4,:],dims=1))),ceil(maximum(abs.(SWV[:,4,:].-reverse(Obs_SWD_M[:,4,:],dims=1))*100 ./reverse(Obs_SWD_M[:,4,:],dims=1)))] )
    # Put it all together
    Plots.plot(p1,p5,p2,p6,p3,p7,p4,p8, layout=(4,2),size=(900,1200),bottom_margin = 0.5Plots.cm, left_margin=0.5Plots.cm,right_margin=0.25Plots.cm,top_margin=0.25Plots.cm)

    savefig("figs/SW_Misfit.pdf")
    # THE END
    println(" ")
    printstyled(" C'est Fini",color=:green)
    println(" ")
     #
return
end
