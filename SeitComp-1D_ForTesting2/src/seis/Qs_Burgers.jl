#By Mariano Arnaiz
"Qs_Burgers: Compute attenuation of S waves (Qs) from the Burgers formulations
Use with care.
This package was developed by Mariano Arnaiz (marianoarnaiz@gmail.com)
at the Universidad Complutense de Madrid
in spring 2022 as part of the WINTER_C project with Dr. Javier Fullea."
## So. Lets begin
using QuadGK, Interpolations, JLD2
#include("Qs_Burguers_Creep.jl")

## Load Grain size functions
@load "src/seis/eGrain_Size_S.jld" # Grain size according to Schierjott 2020
@load "src/seis/eGrain_Size_D.jld" # Grain size according to Dannberg 2017

function Qs_Burgers(zTPC::Array{Any, 2},Grain_Model::String="C")
## Reference Frequency
Period=1; # Central Period in Seconds
Period_ref=1;
## Find T,P and grain size the Mantle
## Find T,P and grain size the Mantle
# get the index of the Mantle layes
indx_m=findall(x->x=="M", zTPC[:,8])
# Get T and P
Tv=zTPC[indx_m,2]; # Temperature in K
Pv=zTPC[indx_m,3]*1e-4 # Pressure in GPa

## Choose Grain size
# evaluate for the grain size
if Grain_Model == "S"
    d=eGrain_Size_S.(zTPC[indx_m,1])
elseif Grain_Model == "D"
    d=eGrain_Size_D.(zTPC[indx_m,1])
elseif Grain_Model == "C"
    d=ones(size(Tv))*8e-3; # Test Grain size
end

# For each mantle layer, computes Qs
Qsv=zeros(size(Tv)); # Qs empty vector
J1v=zeros(size(Tv)); # J1 empty vector
J2v=zeros(size(Tv)); # J1 empty vector

for i=1:size(Tv,1)
	#for S
	if Grain_Model == "S" && Pv[i] <= 24
		Period=0.1
	else
		Period=Period_ref;
	end

	#For D
	if Grain_Model == "D" && Pv[i] <= 24
		Period=0.25
	elseif Grain_Model == "D" && Pv[i] > 24
		Period=0.5;
	end

    Qs,J1,J2 = Qs_Burguers_Creep(Tv[i],Pv[i],d[i],Period,Grain_Model)
    Qsv[i]=Qs;
    J1v[i]=J1;
    J2v[i]=J2;
    #println(" Qs at index $i")
end

# Pad the vectors to include the none mantle nodes
Qsv=[ones(indx_m[1]-1)*500; Qsv];
J1v=[zeros(indx_m[1]-1); J1v];
J1v=[zeros(indx_m[1]-1); J2v];

return Qsv, J1v, J2v
end

## FUNCTIONS TO COMPUTE QS FROM BURGERS
## CREEP MAIN FUNCTION
## Compute Qs, J1 and J2 terms
## FUNCTIONS TO COMPUTE QS FROM BURGERS
function Qs_Burguers_Creep(T,P,d,Period,Grain_Model)

# Compute Angular frequency
ω=2*pi/Period# angular frequency

# Constants and references
E✴=[3.75e5 2.31e5 2.70e5 2.86e5]#3.6e5 #Activation energy E ∗ (J/mol)
mv=3.0  #Viscous grain size exponent
ma=1.31  #Anelastic grain size exponent
τHR=1.0e7  #Reference upper HTB period τHR
τLR=1e-3  #Reference lower HTB period τLR
τMR=3.02e7 #Reference Maxwell period τMR
τPR=3.98e-4 #Reference peak period τPR
α=0.274 #Anelastic frequency exponent
if Grain_Model == "D"
    ΔB=[1.04 1.04 1.04 0.04] #Burgers element strength ∆B
    V✴=[6e-6 6e-6 6e-6 1.0e-6] # # Activation volume V ∗ (m3/mol)
elseif Grain_Model == "S"
    E✴=[3.25e5 2.0e5 2.0e5 2.86e5]#3.6e5 #Activation energy E ∗ (J/mol)
    ΔB=[1.04 1.04 1.04 0.15]
    V✴=[16e-6 16e-6 16e-6 1.0e-6] # # Activation volume V ∗ (m3/mol)
elseif Grain_Model == "C"
    ΔB=[1.04 1.04 1.04 0.2]
    V✴=[6e-6 6e-6 6e-6 1.0e-6] # # Activation volume V ∗ (m3/mol)
end
ΔP=0.057 #Peak height ∆P
𝛔=4.0 #Peak width
R=8.31446261815324 #Gas universal cte
TR=1173.15  #Reference temperature TR (K)
PR=0.2  #Reference pressure PR
dR=13.4e-6 #rReference grainsize dR

# Choose a mineral according to depth
if P <= 14
    Mantle_Mineral=1; # Olivine in the Lithosphere and Asthenosphere
elseif P > 14 && P <= 17.5
    Mantle_Mineral=2; # Wadsleyite in the Upper tranzitional mantle
elseif P > 17.5 && P <= 24
    Mantle_Mineral=3; # Ringwoodite in the lower tranzitional mantle
elseif P > 24
    Mantle_Mineral=4; # perovskite in the Lower mantle
end
# Fist compute the 4 τ: H,L,P,M
# τH
τH= τHR * ((d/dR)^ma) * exp( (E✴[Mantle_Mineral]/R)*((1/T) - (1/TR)) + (V✴[Mantle_Mineral]/R)*((P/T) - (PR/TR)) );
# τL
τL= τLR * ((d/dR)^ma) * exp( (E✴[Mantle_Mineral]/R)*((1/T) - (1/TR)) + (V✴[Mantle_Mineral]/R)*((P/T) - (PR/TR)) );
# τP
τP= τPR * ((d/dR)^ma) * exp( (E✴[Mantle_Mineral]/R)*((1/T) - (1/TR)) + (V✴[Mantle_Mineral]/R)*((P/T) - (PR/TR)) );
# τP
τM= τMR * ((d/dR)^mv) * exp( (E✴[Mantle_Mineral]/R)*((1/T) - (1/TR)) + (V✴[Mantle_Mineral]/R)*((P/T) - (PR/TR)) );

# Define all the functions inside the integrals of J1 and J2
# Inside J1
# Inside the first integral
J1_anel(τi,ωi)=(τi^(α-1))/(1+(ωi^2)*(τi^2));
# Inside the Second integral
J1_p(τi,ωi,τPi)=(1/(τi*(1+(ωi*τi)^2))) * exp(-0.5*(log(τi/τPi)/𝛔)^2);
# Inside J2
# Inside the first integral
J2_anel(τi,ωi)=(τi^(α))/(1+(ωi^2)*(τi^2));
# Inside the Second integral
J2_p(τi,ωi,τPi)=(1/((1+(ωi*τi)^2))) * exp(-0.5*(log(τi/τPi)/𝛔)^2);

# Solve the 4 integrals in one step
#For J1
ij1, err1 = quadgk(𝜏 -> J1_anel(𝜏,ω), τL, τH)
ip1, err3=quadgk(𝜏 ->  J1_p(𝜏,ω,τP),0.0,Inf)
# Jor J2
ij2, err2 = quadgk(𝜏 -> J2_anel(𝜏,ω), τL, τH)
ip2, err4=quadgk(𝜏 ->  J2_p(𝜏,ω,τP),0.0,Inf)

# Now we compute J1 and J2
# J1
J1= 1 + (((α*ΔB[Mantle_Mineral])/((τH^α) - (τL^α))) * ij1) +  (ΔP/(𝛔*sqrt(2*pi)))*ip1;
# J2
J2= (((ω*α*ΔB[Mantle_Mineral])/((τH^α) - (τL^α))) * ij2) +  ((ω*ΔP/(𝛔*sqrt(2*pi)))*ip2) + (1/(ω*τM));

# Compute Qs!
Qs=J1/J2

# Make sure values are not  silly
if Qs> 500
    Qs = 500;
end

if Qs< 50
    Qs = 50;
end
    return Qs,J1,J2
end

## A function to attenuate velocities as a function of the
# extended Burgers model of linear viscoelasticity with creep function:
# dynamic compliance J*(ω). See Jackson and Faul 2010 for details.
function Attenuate_V(Vp,Vs,Rho,J1,J2)
      # Shear modulus with attenuation
      Mu_mod=((J2^2 + J1^2)^(-0.5))*((Vs^2)*Rho )
      #Anharmonic bulk modulus
      K_mod= (((Vp^2)*Rho) -(4/3)*((Vs^2)*Rho ))
      #Attenuated Vs
      Vs_a = sqrt(Mu_mod/Rho)
      #Attenuated Vp
      Vp_a = sqrt((K_mod+(4/3)*Mu_mod)/Rho)
      # in GPa
      Mu_mod=Mu_mod*1e-9;
      K_mod=K_mod*1e-9;
       return Vp_a,Vs_a,K_mod,Mu_mod
 end

## A function to apply the corrections by Karoto to all the Model
function Attenuate_Model(zTPC,Model)
	# Get all the Mantle nodes
	indx_m=findall(x->x=="M", zTPC[:,8]);
	# For each layer
	for i = 1:size(indx_m,1)
		## Correct the model accordingly
		Model[indx_m[i],2]=Model[indx_m[i],2]*(1-(2/9)*cot(pi*0.274)*(1/Qsv[indx_m[i]])) #Vp
		Model[indx_m[i],3]=Model[indx_m[i],3]*(1-(1/2)*cot(pi*0.274)*(1/Qsv[indx_m[i]])) #Vs
	end
	return Model
end

## Compute JU
# This gets the unrelaxed modulus at T,P conditions based on constants from study

 function GetJu(T,P,Vs=nothing,Rho=nothing)
     if Vs == nothing || Rho == nothing
         GUR=	66.5e+9; #Shear modulus at TR, PR 						Pa
         dGdT=  -0.0136e+9 #T-derivative of Shear modulus Pa K^-1
         dGdP=	1.8 # P-derivative of Shear modulus
         TR=1173.15  #Reference temperature TR (K)
         PR=0.2  #Reference pressure PR
         Ju=1/(GUR + dGdT*(T-TR) + dGdP*(P*1e9-PR*1e9))
     else
         Ju=1/(Rho*(Vs^2))
     end
     return Ju
 end

## QL6
# Durek & Ekstrom 1996, A radial model of anelasticity consistent with long-period surface-wave attenuation
QL6=[       0          300
         24.4          300
         24.4          191
           80          191
           80           70
          220           70
          220          165
          670          165
          670          355
         2891          355];

## QML9: Lawrence & Wysession 2006, qlm9: A new radial quality factor (Qu) model for the lower mantle
QLM9=[     0         600
          80         600
          80          80
         220          80
         220         143
         400         143
         400         276
         670         276
         670         362
        1000         362
        1000         325
        1350         325
        1350         287
        1700         287
        1700         307
        2050         307
        2050         383
        2400         383
        2400         459
        2700         459
        2700         452
        2800         452
        2800         278
        2891         278];
## SL8; ; SL8 (Anderson and Hart, 1978)

SL8=[  11  500
   11  100
  200   90
  421   95
  421  330
  671  353
 2200  366
 2400   92
 2843   92
 2891 92];

## OS08:

 OS08=[0 101
712 101
712 621
2900 621];

## QM1

QM1=[0.0  120.0
 88.0  100.0
190.4  116.0
318.5  121.0
421.0  152.0
546.0  210.0
671.0  286.0
878.8  362.0
1017.4  383.0
1155.9  379.0
1329.0  368.0
1502.2  360.0
1675.4  343.0
1848.5  311.0
2021.7  285.0
2194.8  274.0
2368.0  272.0
2541.1  271.0
2714.3  268.0
2886.7  266.0];
