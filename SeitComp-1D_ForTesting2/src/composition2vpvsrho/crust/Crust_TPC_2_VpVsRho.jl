"This is a function that helps compute the values of Vp,Vs and Rho
for the Igneous layers. Based on an interpolation of the values given
by MinVel-master"

#using ScatteredInterpolation, JLD2

## Load the Sediments Functions
#@load "func/Sediments_Vp.jld"
#@load "func/Sediments_Vs.jld"
#@load "func/Sediments_Rho.jld"

#data is a matrix: [T(K) P(GPa) QZ(%) Calcite(%) Shales(%)]
function Crust_TPC_2_VpVsRho(data,val=1:3)

#Initiate every thin
properties=(zeros(size(data,1),3));
Vp=0.0;
Vs=0.0;
Rho=0.0;

# Coeficients of the equation
Coef=[2.92e+01	2.40e+02	2.90e+02
1.54e-02	2.00e-02	1.67e-02
2.28e-02	9.83e-03	1.35e-02
-2.66e+00	-2.87e+01	-3.43e+01
1.31e-01	1.39e+00	1.65e+00
-3.28e-03	-3.45e-02	-4.05e-02
4.40e-05	4.62e-04	5.42e-04
-3.02e-07	-3.21e-06	-3.76e-06
8.40e-10	9.05e-09	1.06e-08
3.66e-01	2.48e+00	2.71e+00
-4.96e-02	-1.71e-01	-3.50e-01
7.43e-03	6.20e-03	1.18e-02
-3.07e-04	-5.35e-04	-5.52e-04
-8.73e-03	-7.86e-02	-8.33e-02
3.92e-05	5.83e-04	6.72e-04
-1.34e-03	4.39e-03	1.01e-02
2.37e-05	-2.82e-05	-9.34e-05];

#for each line
for i=1:size(data,1)

   # Compute Vp from the composition and account for p and Na2O if necessary
   if data[i,1] < 12000
       Vp = 6.90 - 0.011*data[i,2] + 0.037*data[i,3]+ 0.045*data[i,4];
   elseif data[i,1] >= 12000
      # Na2O(z)
       Na2O= 5.522e-5*((data[i,5]/1000)^3) -0.0061*((data[i,5]/1000)^2) + 0.2007*(data[i,5]/1000) + 1.8111;
       Vp = Coef[1,1] +Coef[2,1]*data[i,3]+Coef[3,1]*data[i,4]+Coef[4,1]*data[i,2]+Coef[5,1]*data[i,2]^2 +Coef[6,1]*data[i,2]^3 +Coef[7,1]*data[i,2]^4 + Coef[8,1]*data[i,2]^5 + Coef[9,1]*data[i,2]^6 + Coef[10,1]*Na2O + Coef[11,1]*Na2O^2 + Coef[12,1]*Na2O^3 + Coef[13,1]*Na2O^4 + Coef[14,1]*data[i,2]*Na2O + Coef[15,1]*(data[i,2]^2)*Na2O + Coef[16,1]*data[i,2]*(Na2O^2) + Coef[17,1]*(data[i,2]^2)*(Na2O^2);
   end

         Vs=Vp2Vs(Vp);
         Rho=Vp2Rho(Vp);
         # Safe to output
         properties[i,1]=Vp;
         properties[i,2]=Vs;
         properties[i,3]=Rho;
  
end
return properties[:,val]
end

## Rho to Vp: From Brotcher 2005
# Density expected in g/cm^3
function Rho2Vp(rho)
Vp=  39.128*rho - 63.0648*rho^2 + 37.083*rho^3 - 9.1819*rho^4 + 0.8228*rho^5
    return Vp
end
##
function Vp2Vs(Vp)
Vs= 0.7858 - 1.2344*Vp + 0.7949*Vp^2 - 0.1238*Vp^3 + 0.0064*Vp^4
    return Vs
end

## Vp to Rho
function Vp2Rho(Vp)
    Rho = 1.6612*Vp - 0.4721*Vp^2 + 0.0671*Vp^3 - 0.0043*Vp^4 + 0.000106*Vp^5;
    return Rho
end
