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
Period=50; # Central Period in Seconds
ω=2*pi/Period# angular frequency

## Find T,P and grain size the Mantle

# get the index of the Mantle layes
indx_m=findall(x->x=="M", zTPC[:,8])
# Get T and P
Tv=zTPC[indx_m,2]; # Temperature in K
Pv=zTPC[indx_m,3]*1e-4 # Pressure in GPa
# evaluate for the grain size
if Grain_Model == "S"
    gs=eGrain_Size_S.(zTPC[indx_m,1])
elseif Grain_Model == "D"
    gs=eGrain_Size_D.(zTPC[indx_m,1])
end

## For each mantle layer, computes Qs
Qsv=zeros(size(Tv)); # Qs empty vector
J1v=zeros(size(Tv)); # J1 empty vector
J2v=zeros(size(Tv)); # J1 empty vector

for i=1:size(Tv,1)
    Qs,J1,J2 = Qs_Burguers_Creep(Tv[i],Pv[i],gs[i],ω)
    Qsv[i]=Qs;
    J1v[i]=J1;
    J2v[i]=J2;
    #println(" Qs at index $i")
end

return Qsv, J1v, J2v
end

#Plots.plot(Qsv,zTPC[indx_m,1],yflip=true)
#Vp_a,Vs_a,K_mod,Mu_mod= Attenuate_V(Vp,Vs,Rho,J1,J2)

## FUNCTIONS TO COMPUTE QS FROM BURGERS
## CREEP MAIN FUNCTION
## Compute J1 and J2 terms
function Qs_Burguers_Creep(Ti0,presi0,gsi0,omegai0)
# Constants and references
Tr=1173.0  #Reference temperature (K)
Pr=0.2  #Reference pressure (GPa)
gsr=1.34e-5 #ref grain size
deltaB=1.04 #Burgers element strength
alpha=0.274 #Anelastic frequency exponent
tauLo=1e-3  #reference lower High Temperature Background period (s)
tauHo=1.0e7  #reference higher High Temperature Background period (s)
tauMo=3.02e7 #reference Maxwell period (s)
ma=1.31 #Anelastic grain size exponent
mv=3.0 #Viscous grain size exponent
#EB=3.75e5 #Activation energy
#AV=6E-6 #425e-5 #Activation volume
EB=3.75e5#3.6e5 #Activation energy (J/mol)
AV= 1.2300e-5 #6e-6 # # Activation volume (m3/mol)
tauPo=3.98e-4 #reference peak period
deltaP=0.057 #Peak hight
sig=4.0 #Peak width
R=8.31446261815324 #Gas universal cte

# Basic calculations
iTr=1/Tr
PT=Pr/Tr
AVR=AV/R
ER=EB/R
gr=gsi0/gsr #Normalized grain size

# Keep on going
#cp=deltaP*(2*pi)^((-0.5)/sig)
cp=deltaP*(2*pi)^(-0.5)/sig
taut=exp((ER)*(1.0/Ti0-iTr))*exp(AVR*(presi0/Ti0-PT))
  if taut > 1e+26
      taut=1e+26
  end  # this is for numerical reasons only (larger taut's give error)
  tauH=tauHo*gr^ma*taut
  tauL=tauLo*gr^ma*taut
  tauP=tauPo*gr^ma*taut
  tauM=tauMo*gr^mv*taut

# Go for the integrals
#Initialize results
on=1;

# original values with deifned tolerance
# ij1, err1 = quadgk(x -> J1anel(x,omegai0), tauL, tauH, rtol=1.0e-10,atol=1.0e-16, order=1)
# ij2, err2 = quadgk(x -> J2anel(x,omegai0), tauL, tauH, rtol=1.0e-10,atol=1.0e-16, order=1)
# ip1, err3=quadgk(x ->  J1p(x,omegai0,tauP),0.0,Inf,atol=1.e-5,rtol=1.e-4)
# ip2, err4=quadgk(x ->  J2p(x,omegai0,tauP),0.0,Inf,atol=1.e-5,rtol=1.e-4)

ij1, err1 = quadgk(x -> J1anel(x,omegai0), tauL, tauH)
ij2, err2 = quadgk(x -> J2anel(x,omegai0), tauL, tauH)
ip1, err3=quadgk(x ->  J1p(x,omegai0,tauP),0.0,Inf)
ip2, err4=quadgk(x ->  J2p(x,omegai0,tauP),0.0,Inf)

#qagi(J2p,0.0,1,1.e-5,1.e-4,ip2,abserr,neval,ier,omegai0,tauP)
#
# Go for the J1 and J2
Jb1=alpha*deltaB*ij1/(tauH^alpha-tauL^alpha)
Jb2=omegai0*alpha*deltaB*ij2/(tauH^alpha-tauL^alpha)
Jp1=cp*ip1
Jp2=cp*omegai0*ip2
J1=1+Jp1+Jb1
J2=(Jb2+Jp2)+1.0/(omegai0*tauM)

Qs=J1/J2

# End of function, return my values!
      return Qs,J1,J2
end

## The next few functions define the Burguers compliance terms
# as callable functiones that we need for the integrals.
function J1anel(tau,omega)
    alpha=0.274
    J1anel=(tau^(alpha-1))/(1+(omega^2)*(tau^2))
    return J1anel
end

function J2anel(tau,omega)
       alpha=0.274
       J2anel=(tau^alpha)/(1+(omega^2)*(tau^2))
    return J2anel
end

function J1p(tau,omega,tauP)
    sig=4.0
    J1p=(1.0/tau)*exp(-0.5*(log(tau/tauP)/sig)^2)/(1+(omega*tau)^2)
   return J1p
end

function J2p(tau,omega,tauP)
       sig=4.0
       J2p=exp(-0.5*(log(tau/tauP)/sig)^2)/(1+(omega*tau)^2)
    return J2p
end
