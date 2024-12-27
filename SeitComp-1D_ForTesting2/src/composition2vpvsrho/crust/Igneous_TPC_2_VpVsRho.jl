"This is a function that helps compute the values of Vp,Vs and Rho
for the Igneous (crustal) layers. Based on an interpolation of the values given
by MinVel-master"

using Interpolations, JLD2

## Load the Sediments Functions
@load "src/composition2vpvsrho/crust/Felsic_Vp.jld"
@load "src/composition2vpvsrho/crust/Felsic_Vs.jld"
@load "src/composition2vpvsrho/crust/Felsic_Rho.jld"
@load "src/composition2vpvsrho/crust/Mafic_Vp.jld"
@load "src/composition2vpvsrho/crust/Mafic_Vs.jld"
@load "src/composition2vpvsrho/crust/Mafic_Rho.jld"
@load "src/composition2vpvsrho/crust/Ultramafic_Vp.jld"
@load "src/composition2vpvsrho/crust/Ultramafic_Vs.jld"
@load "src/composition2vpvsrho/crust/Ultramafic_Rho.jld"

#Input is a matrix: [T(K) P(GPa) Mineral1(%) Mineral2(%) Mineral3(%) "Type_of Rock"]
#Felsic rocks: [T(K) P(GPa) Qz(%) Kfeldespar(%) Plagioclase(%) "Felsic"]
#Mafic rocks: [T(K) P(GPa) An(%) Cpx(%) Opc(%) "Mafic"]
#Ultramafic rocks: [T(K) P(GPa) Ol(%) Cpx(%) Opc(%) "Ultramafic"]
function Igneous_TPC_2_VpVsRho(data::Array{Any, 2})

properties=(zeros(Float64,size(data,1),3));
for i=1:size(data,1)

   if sum(data[i,3:5],dims=1)[1] == 100.0
      composition_index=findall( (COMPOSITION_MATRIX[:,1].==data[i,3]) .& (COMPOSITION_MATRIX[:,2].==data[i,4]) .& (COMPOSITION_MATRIX[:,3].==data[i,5]))[1]
      if data[i,6] == "Felsic"
         properties[i,1]=Felsic_Vp(data[i,1],data[i,2],composition_index);
         properties[i,2]=Felsic_Vs(data[i,1],data[i,2],composition_index);
         properties[i,3]=Felsic_Rho(data[i,1],data[i,2],composition_index);
      elseif data[i,6] == "Mafic"
         properties[i,1]=Mafic_Vp(data[i,1],data[i,2],composition_index);
         properties[i,2]=Mafic_Vs(data[i,1],data[i,2],composition_index);
         properties[i,3]=Mafic_Rho(data[i,1],data[i,2],composition_index);
      elseif data[i,6] == "Ultramafic"
         properties[i,1]=Ultramafic_Vp(data[i,1],data[i,2],composition_index);
         properties[i,2]=Ultramafic_Vs(data[i,1],data[i,2],composition_index);
         properties[i,3]=Ultramafic_Rho(data[i,1],data[i,2],composition_index);
      end
   else
      println("There is a problem with the composition at line $i")
   end
end
return properties
end

## Composition Matrix for index search
COMPOSITION_MATRIX=[    0     0   100
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
