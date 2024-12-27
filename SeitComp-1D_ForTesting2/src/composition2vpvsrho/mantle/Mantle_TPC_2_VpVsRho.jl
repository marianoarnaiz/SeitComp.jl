"This is a function that helps compute the values of Vp,Vs and Rho
for the mantle Olivine. Based on an interpolation of the values given
by Perple_X and Perple_X_Litmod"

using Interpolations, JLD2

## Load the Mantle Functions
@load "src/composition2vpvsrho/mantle/Mantle_Vp.jld"
@load "src/composition2vpvsrho/mantle/Mantle_Vs.jld"
@load "src/composition2vpvsrho/mantle/Mantle_Rho.jld"

function Mantle_TPC_2_VpVsRho(input,val=1:3)

properties=(zeros(Float64,size(input,1),3));
for i=1:size(input,1)

   if input[i,2] < 10000
      input[i,2] = 10000;
   end

   Vp=Mantle_Vp(input[i,1],input[i,2],input[i,3],input[i,4]);

   if Vp > 14.0 || Vp < 7.0
      properties[i,1]=NaN;
      properties[i,2]=NaN;
      properties[i,3]=NaN;
   else
      properties[i,1]=Vp
      properties[i,2]=Mantle_Vs(input[i,1],input[i,2],input[i,3],input[i,4]);
      #properties[i,2]=Vs
      properties[i,3]=Mantle_Rho(input[i,1],input[i,2],input[i,3],input[i,4]);
      #properties[i,3]=Rho
   end
end
return properties[:,val]
end
