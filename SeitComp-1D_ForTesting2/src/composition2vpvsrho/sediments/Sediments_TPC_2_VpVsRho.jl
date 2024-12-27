"This is a function that helps compute the values of Vp,Vs and Rho
for the Sedimentary layers. Based on an interpolation of the values given
by MinVel-master"

using ScatteredInterpolation, JLD2

## Load the Sediments Functions
@load "src/composition2vpvsrho/sediments/Sediments_Vp.jld"
@load "src/composition2vpvsrho/sediments/Sediments_Vs.jld"
@load "src/composition2vpvsrho/sediments/Sediments_Rho.jld"

#data is a matrix: [T(K) P(GPa) QZ(%) Calcite(%) Shales(%)]
function Sediments_TPC_2_VpVsRho(data)

#Initiate every thin
properties=(zeros(size(data,1),3));
Vp=0.0;
Vs=0.0;
Rho=0.0;

#for each line
for i=1:size(data,1)
   #Check the compositons
   if sum(data[i,3:5]) == 100
      # Compute Vp
      V=evaluate(Sediments_Vp, [data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]]);
      Vp=V[1]
      if Vp[1] > 8
         Vp=NaN; Vs=NaN; Rho=NaN;
      else
         #Compute Vs and Rho
         V2=evaluate(Sediments_Vs, [data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]]);
         Vs=V2[1];
         R=evaluate(Sediments_Rho, [data[i,1];data[i,2];data[i,3];data[i,4];data[i,5]]);
         Rho=R[1];
## Correct Massive Vp,Vs and Rho by porotisy percentage
#Consider a general gas like air that fills the rock
         Vpfluid=20*sqrt(data[i,1])*0.001
         Vsfluid=0.001
         # Raymer et al. (1980)
         Vp=((1-data[i,6]*0.01)^2)*Vp+(data[i,6]*0.01)*Vpfluid;
         Vs=((1-data[i,6]*0.01)^2)*Vs+(data[i,6]*0.01)*Vsfluid
         #Athy 1930
         Rho=((data[i,6]*0.01)-1)*-Rho;

         # Safe to output
         properties[i,1]=Vp;
         properties[i,2]=Vs;
         properties[i,3]=Rho;

      end
      else
         println("There is a problem with the composition at linr $i")
      end
end
return properties
end
