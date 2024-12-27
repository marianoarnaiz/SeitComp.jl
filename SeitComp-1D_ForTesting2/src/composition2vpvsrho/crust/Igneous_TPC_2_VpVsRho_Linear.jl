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

      if data[i,6] == "Felsic"
         properties[i,1]=Felsic_Vp(data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]);
         properties[i,2]=Felsic_Vs(data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]);
         properties[i,3]=Felsic_Rho(data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]);
      elseif data[i,6] == "Mafic"
         properties[i,1]=Mafic_Vp(data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]);
         properties[i,2]=Mafic_Vs(data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]);
         properties[i,3]=Mafic_Rho(data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]);
      elseif data[i,6] == "Ultramafic"
         properties[i,1]=Ultramafic_Vp(data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]);
         properties[i,2]=Ultramafic_Vs(data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]);
         properties[i,3]=Ultramafic_Rho(data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]);
      end
   else
      println("There is a problem with the composition at line $i")
   end
end
return properties
end
