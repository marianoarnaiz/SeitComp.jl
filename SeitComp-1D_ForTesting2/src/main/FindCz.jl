"This function will plot for you 3 ternary diagrams
for Vp, Vs and Rho for a specific composition. I hope it helps in your modelling"

using ScatteredInterpolation,Interpolations, JLD2,GMT

## Load the Igneous Functions
@load "src/composition2vpvsrho/crust/Felsic_Vp.jld"
@load "src/composition2vpvsrho/crust/Felsic_Vs.jld"
@load "src/composition2vpvsrho/crust/Felsic_Rho.jld"
@load "src/composition2vpvsrho/crust/Mafic_Vp.jld"
@load "src/composition2vpvsrho/crust/Mafic_Vs.jld"
@load "src/composition2vpvsrho/crust/Mafic_Rho.jld"
@load "src/composition2vpvsrho/crust/Ultramafic_Vp.jld"
@load "src/composition2vpvsrho/crust/Ultramafic_Vs.jld"
@load "src/composition2vpvsrho/crust/Ultramafic_Rho.jld"
## Load the Sediments Functions
@load "src/composition2vpvsrho/sediments/Sediments_Vp.jld"
@load "src/composition2vpvsrho/sediments/Sediments_Vs.jld"
@load "src/composition2vpvsrho/sediments/Sediments_Rho.jld"
## Load the Mantle Functions
@load "src/composition2vpvsrho/mantle/Mantle_Vp.jld"
@load "src/composition2vpvsrho/mantle/Mantle_Vs.jld"
@load "src/composition2vpvsrho/mantle/Mantle_Rho.jld"

function FindCz(z,Vp,Vs,Rho,RockType,SecondInput=nothing)
# z is depth in km, in the computed model that is saved in zTPC in K
# RockType is "S" or "I" for sedimentary of Igneous rocks
# SecondInput is Porosity for "S" and "Felsic"/"Mafic"/"Ultramafic" of "I"

# Test for this depth
#z=10;
zindx=findall(zTPC[:,1].==z);
zindx=zindx[1];
T=zTPC[zindx,2]; # Temperature in K
P=zTPC[zindx,3]*1e-4; # Pressure in GPa
#Vp=4.6; Vs=2.8; Rho=2.35;

if RockType=="S" && typeof(SecondInput) == Float64 || typeof(SecondInput) == Int64
   DATA_MATRIX=Sediments_TPC_2_VpVsRho([ones(66)*T ones(66)*P Compmix ones(66)*SecondInput]);

   ternary([Compmix DATA_MATRIX[:,1]], frame=(annot=20, ticks=10, grid=:a, alabel="Quartz", blabel="Carbonates", clabel="Clays", suffix=" %"),contour=(annot=[Vp],contour=[Vp],labels=(distance=3,),pen=(1.0,:blue)))
   ternary!([Compmix DATA_MATRIX[:,2]],contour=(annot=[Vs],contour=[Vs],labels=(distance=3,),pen=(1.0,:red)))
   ternary!([Compmix DATA_MATRIX[:,3]], fmt=:pdf,contour=(annot=[Rho],contour=[Rho],labels=(distance=3,),pen=(1.0,:green)),savefig="figs/FindC_Result.pdf")
   printstyled("Your Triangle with Contours is in /figs",color=:light_cyan)

elseif RockType=="I" && typeof(SecondInput) == String
   w = Vector{String}()
   for i=1:66
      push!(w,SecondInput)
   end
   DATA_MATRIX=Igneous_TPC_2_VpVsRho([ones(66)*T ones(66)*P Compmix w]);
   if SecondInput == "Felsic"

      ternary([Compmix DATA_MATRIX[:,1]], frame=(annot=20, ticks=10, grid=:a, alabel="Quartz", blabel="K_Feldespar", clabel="Plagioclase", suffix=" %"),contour=(annot=[Vp],contour=[Vp],labels=(distance=3,),pen=(1.0,:blue)))
      ternary!([Compmix DATA_MATRIX[:,2]],contour=(annot=[Vs],contour=[Vs],labels=(distance=3,),pen=(1.0,:red)))
      ternary!([Compmix DATA_MATRIX[:,3]], fmt=:pdf,contour=(annot=[Rho],contour=[Rho],labels=(distance=3,),pen=(1.0,:green)),savefig="figs/FindC_Result.pdf")
      printstyled("Your Triangle with Contours is in /figs",color=:light_cyan)

   elseif SecondInput == "Mafic"

      ternary([Compmix DATA_MATRIX[:,1]], frame=(annot=20, ticks=10, grid=:a, alabel="Anorthosite", blabel="Clinopyroxene", clabel="Ortopyroxene", suffix=" %"),contour=(annot=[Vp],contour=[Vp],labels=(distance=3,),pen=(1.0,:blue)))
      ternary!([Compmix DATA_MATRIX[:,2]],contour=(annot=[Vs],contour=[Vs],labels=(distance=3,),pen=(1.0,:red)))
      ternary!([Compmix DATA_MATRIX[:,3]], fmt=:pdf,contour=(annot=[Rho],contour=[Rho],labels=(distance=3,),pen=(1.0,:green)),savefig="figs/FindC_Result.pdf")
      printstyled("Your Triangle with Contours is in /figs",color=:light_cyan)

   elseif SecondInput == "Ultramafic"

      ternary([Compmix DATA_MATRIX[:,1]], frame=(annot=20, ticks=10, grid=:a, alabel="Olivine", blabel="Clinopyroxene", clabel="Ortopyroxene", suffix=" %"),contour=(annot=[Vp],contour=[Vp],labels=(distance=3,),pen=(1.0,:blue)))
      ternary!([Compmix DATA_MATRIX[:,2]],contour=(annot=[Vs],contour=[Vs],labels=(distance=3,),pen=(1.0,:red)))
      ternary!([Compmix DATA_MATRIX[:,3]], fmt=:pdf,contour=(annot=[Rho],contour=[Rho],labels=(distance=3,),pen=(1.0,:green)),savefig="figs/FindC_Result.pdf")
      printstyled("Your Triangle with Contours is in /figs",color=:light_cyan)

   end

   elseif RockType=="M" #P *1e4 because the value is expected in Mbar
      DATA_MATRIX=Mantle_TPC_2_VpVsRho([T*ones(220) P*ones(220)*1e4 Compmi_M]);
      grid1=gmt("xyz2grd -R1.0/11.0/1.0/20.0 -I1", [Compmi_M[:,1] Compmi_M[:,2] DATA_MATRIX[:,1]])
      grid2=gmt("xyz2grd -R1.0/11.0/1.0/20.0 -I1", [Compmi_M[:,1] Compmi_M[:,2] DATA_MATRIX[:,2]])
      grid3=gmt("xyz2grd -R1.0/11.0/1.0/20.0 -I1", [Compmi_M[:,1] Compmi_M[:,2] DATA_MATRIX[:,3]])
      fig1=Plots.contour(grid1,aspect_ratio = 0.5,xlabel = "wt% Al2O3", ylabel = "wt% FeO",title = "Mantle Properties",levels=[Vp],xlim=(1.0,11.0),ylim=(1.0,20.0),xticks=(collect(1:1:11)),yticks=(collect(1:1:20)),linecolor=:blue,contour_labels=true,legend = :none)
      Plots.contour!(grid2,aspect_ratio = 0.5,xlim=(1.0,11.0),ylim=(1.0,20.0),xticks=(collect(1:1:11)),levels=[Vs],yticks=(collect(1:1:20)),linecolor=:red,contour_labels=true,legend = :none)
      Plots.contour!(grid3,aspect_ratio = 0.5,xlim=(1.0,11.0),ylim=(1.0,20.0),xticks=(collect(1:1:11)),levels=[Rho],yticks=(collect(1:1:20)),linecolor=:green,contour_labels=true,show=false,legend = :none)
      savefig(fig1,"figs/FindC_Result.pdf")
      printstyled("Your Map with Contours is in /figs",color=:light_cyan)


else
   printstyled("There is an incompatibility among the inputs",color=:red)
end


return DATA_MATRIX
end


## Composition for ternary diagrams
Compmix=[    0     0   100
     0    10    90
     0    20    80
     0    30    70
     0    40    60
     0    50    50
     0    60    40
     0    70    30
     0    80    20
     0    90    10
     0   100     0
    10     0    90
    10    10    80
    10    20    70
    10    30    60
    10    40    50
    10    50    40
    10    60    30
    10    70    20
    10    80    10
    10    90     0
    20     0    80
    20    10    70
    20    20    60
    20    30    50
    20    40    40
    20    50    30
    20    60    20
    20    70    10
    20    80     0
    30     0    70
    30    10    60
    30    20    50
    30    30    40
    30    40    30
    30    50    20
    30    60    10
    30    70     0
    40     0    60
    40    10    50
    40    20    40
    40    30    30
    40    40    20
    40    50    10
    40    60     0
    50     0    50
    50    10    40
    50    20    30
    50    30    20
    50    40    10
    50    50     0
    60     0    40
    60    10    30
    60    20    20
    60    30    10
    60    40     0
    70     0    30
    70    10    20
    70    20    10
    70    30     0
    80     0    20
    80    10    10
    80    20     0
    90     0    10
    90    10     0
   100     0     0];

## Composition for Mapping Mantle properties
Compmi_M=[ 1.0	1.0
1.0	2.0
1.0	3.0
1.0	4.0
1.0	5.0
1.0	6.0
1.0	7.0
1.0	8.0
1.0	9.0
1.0	10.0
1.0	11.0
1.0	12.0
1.0	13.0
1.0	14.0
1.0	15.0
1.0	16.0
1.0	17.0
1.0	18.0
1.0	19.0
1.0	20.0
2.0	1.0
2.0	2.0
2.0	3.0
2.0	4.0
2.0	5.0
2.0	6.0
2.0	7.0
2.0	8.0
2.0	9.0
2.0	10.0
2.0	11.0
2.0	12.0
2.0	13.0
2.0	14.0
2.0	15.0
2.0	16.0
2.0	17.0
2.0	18.0
2.0	19.0
2.0	20.0
3.0	1.0
3.0	2.0
3.0	3.0
3.0	4.0
3.0	5.0
3.0	6.0
3.0	7.0
3.0	8.0
3.0	9.0
3.0	10.0
3.0	11.0
3.0	12.0
3.0	13.0
3.0	14.0
3.0	15.0
3.0	16.0
3.0	17.0
3.0	18.0
3.0	19.0
3.0	20.0
4.0	1.0
4.0	2.0
4.0	3.0
4.0	4.0
4.0	5.0
4.0	6.0
4.0	7.0
4.0	8.0
4.0	9.0
4.0	10.0
4.0	11.0
4.0	12.0
4.0	13.0
4.0	14.0
4.0	15.0
4.0	16.0
4.0	17.0
4.0	18.0
4.0	19.0
4.0	20.0
5.0	1.0
5.0	2.0
5.0	3.0
5.0	4.0
5.0	5.0
5.0	6.0
5.0	7.0
5.0	8.0
5.0	9.0
5.0	10.0
5.0	11.0
5.0	12.0
5.0	13.0
5.0	14.0
5.0	15.0
5.0	16.0
5.0	17.0
5.0	18.0
5.0	19.0
5.0	20.0
6.0	1.0
6.0	2.0
6.0	3.0
6.0	4.0
6.0	5.0
6.0	6.0
6.0	7.0
6.0	8.0
6.0	9.0
6.0	10.0
6.0	11.0
6.0	12.0
6.0	13.0
6.0	14.0
6.0	15.0
6.0	16.0
6.0	17.0
6.0	18.0
6.0	19.0
6.0	20.0
7.0	1.0
7.0	2.0
7.0	3.0
7.0	4.0
7.0	5.0
7.0	6.0
7.0	7.0
7.0	8.0
7.0	9.0
7.0	10.0
7.0	11.0
7.0	12.0
7.0	13.0
7.0	14.0
7.0	15.0
7.0	16.0
7.0	17.0
7.0	18.0
7.0	19.0
7.0	20.0
8.0	1.0
8.0	2.0
8.0	3.0
8.0	4.0
8.0	5.0
8.0	6.0
8.0	7.0
8.0	8.0
8.0	9.0
8.0	10.0
8.0	11.0
8.0	12.0
8.0	13.0
8.0	14.0
8.0	15.0
8.0	16.0
8.0	17.0
8.0	18.0
8.0	19.0
8.0	20.0
9.0	1.0
9.0	2.0
9.0	3.0
9.0	4.0
9.0	5.0
9.0	6.0
9.0	7.0
9.0	8.0
9.0	9.0
9.0	10.0
9.0	11.0
9.0	12.0
9.0	13.0
9.0	14.0
9.0	15.0
9.0	16.0
9.0	17.0
9.0	18.0
9.0	19.0
9.0	20.0
10.0	1.0
10.0	2.0
10.0	3.0
10.0	4.0
10.0	5.0
10.0	6.0
10.0	7.0
10.0	8.0
10.0	9.0
10.0	10.0
10.0	11.0
10.0	12.0
10.0	13.0
10.0	14.0
10.0	15.0
10.0	16.0
10.0	17.0
10.0	18.0
10.0	19.0
10.0	20.0
11.0	1.0
11.0	2.0
11.0	3.0
11.0	4.0
11.0	5.0
11.0	6.0
11.0	7.0
11.0	8.0
11.0	9.0
11.0	10.0
11.0	11.0
11.0	12.0
11.0	13.0
11.0	14.0
11.0	15.0
11.0	16.0
11.0	17.0
11.0	18.0
11.0	19.0
11.0	20.0];
# Compmi_M=zeros(90,2)
# c1=1;
# c2=1;
# for AL=1.1:1:9.1
#    for FE= 0.1:1:9.1
#          Compmi_M[c1,:]=[AL FE];
#          c1=c1+1;
#    end
# end
