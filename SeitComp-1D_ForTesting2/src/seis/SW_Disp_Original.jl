"SW_Disp: Surface Wave dispersion for Julia
This file is a small package to perform forward modelling of Surface Wave Dispersion.
It is based on the codes by Haney and Tsai (2015,2017,2019) as well as on many comments
Prof. Robert Herrmann and his Computer Programs in Seismology (version 3.30).
Use with care.
This package was developed by Mariano Arnaiz (marianoarnaiz@gmail.com)
at the Universidad Complutense de Madrid
in spring 2021 as part of the WINTER_C project with Dr. Javier Fullea."
## Begin by loading the modules used
# All this are necessary to run SW_DISP properly
# You need to istall "brew install arpack"!!!!!
using PolynomialRoots, LinearAlgebra, Arpack, SparseArrays, Interpolations

## CALLING FUNCTION
"This function only CALLS ALL THE OTHER FUNCTIONS. It is intended to ease the input/output of the data"
##
function SW_Disp(Model,h_index,Obs_T,Nmodes,Earth,Elastic,Anisotropy)

# Print Message for user
#println("Runnning SW_Disp 1.0")
#println("Computing Dispersion for: ",Nmodes," Surface Wave Modes")
#println("Considering:")
# if Earth==1
#     #println(" - Earth's Sphericity")
# else
#     println(" - Flat Earth")
# end
# if Elastic==1
#     println(" - Elastic Wave Propagation")
# else
#     println(" - Anaelastic Wave Propagation (Qs is expected)")
# end
# if Anisotropy==0
#     println(" - Isotropic Earth")
# else
#     println(" - Anisotropic Earth (Radial Anisotropy is expected)")
# end

# Prepare Inputs
fks,Mesh,vpv,vsv,vsv_v,vsv_h,rhov,qsv,qR,Nsolid,hsolid,Nfluid,hfluid,vpfluid,vpfv,rhofv = fix_SWF_input(Obs_T,h_index,Model,Elastic,Anisotropy);
# Compute de Dispersion for R Wave
kkR, cR, UR, evR = Rayleigh_Forwardsp(Nsolid,vsv_v,vpv,rhov,fks,hsolid,Nmodes,Nfluid,vpfv,rhofv,hfluid,qsv,Earth,Elastic);
#println("Running: R wave dispersion")
# Compute de Dispersion for L Wave
kkL, cL, UL, evL = Love_Forwardsp(Nsolid,vsv_h,rhov,fks,hsolid,Nmodes,qsv,Earth,Elastic);
#println("Running: L wave dispersion")
# Organize the Outputs
T, SWV,EG_FUNC_R_Z, EG_FUNC_R_R, EG_FUNC_L_T = fix_SWF_output(Nmodes,cR, UR, evR, cL, UL, evL, fks);
# println("Done! Wait for outputs to be written:")
# println(" - T: Periods (s)")
# println(" - SWV: c and U (km/s) -> SWV[:, 1(cR) 2(UR) 3(cL) 4(UL), Mode]")
# println(" - EG_FUNC_R_Z: R wave Z Eigenfunction")
# println(" - EG_FUNC_R_R: R wave R Eigenfunction")
# println(" - EG_FUNC_L_T: L wave T Eigenfunction")
return T, SWV,EG_FUNC_R_Z, EG_FUNC_R_R, EG_FUNC_L_T
end

## FUNTION FOR INPUT Data
"This function is intended to fix the inputs for both of the codes that require a bunch of inputs"
##
function fix_SWF_input(Obs_T,h_index,Model,Elastic,Anisotropy)

if h_index == 0
    h=GOLD_MESH;
    println(" ")
    printstyled("   - SW_Disp Says: Using Gold Mesh",color=:white)
    println(" ")
elseif h_index == 1
    h=Make_1D_Mesh(minimum(Obs_T));
    println(" ")
    printstyled("   - SW_Disp Says: Computing Mesh for the Current Inputs (T)",color=:white)
    println(" ")
elseif h_index == 2
    h=Make_1D_Mesh(minimum(Obs_T), Depths./1000);
    println(" ")
    printstyled("   - SW_Disp Says: Computing Mesh for the Current Inputs (T,Z)",color=:white)
    println(" ")
end

fks=1 ./collect(reverse(Obs_T)); # vector of  frequencies  at which the velocities are measured (Hz) Change T (s) range
Step=1000; #step for the fluid
Mesh=cumsum(h,dims=1); #from special h to the FEM mesh
Mesh=Mesh[1:end-1]

# Iterpolate all the properties
iVp=LinearInterpolation(Model[:,1]*1000,  Model[:,2]*1000, extrapolation_bc=Flat());
iVs=LinearInterpolation(Model[:,1]*1000,  Model[:,3]*1000, extrapolation_bc=Flat());
iρ=LinearInterpolation(Model[:,1]*1000,  Model[:,4]*1000, extrapolation_bc=Flat());
vpv=iVp(Mesh*1000);
vsv=iVs(Mesh*1000);
rhov=iρ(Mesh*1000);

if Elastic==0
    iQs=LinearInterpolation(Model[:,1]*1000,  Model[:,5], extrapolation_bc=Flat());
    qsv=iQs(Mesh*1000);
else
    qsv=zeros(size(Mesh,1)); #Dummy vaariable
end

if Anisotropy==1
    iR=LinearInterpolation(Model[:,1]*1000,  Model[:,6], extrapolation_bc=Flat());
    qR=iR(Mesh*1000)./100;
    vsv_v=((3 .-qR)./3).*vsv; #Vertical velocity in realtion to Radial Anisitropy
    vsv_h=(1 .+ (2/3).*qR).*vsv; #Horizontal velocity equal to total velocity

else
    qR=zeros(size(Mesh,1)); #Dummy vaariable
    vsv_v=vsv; #Vertical velocity equal to total velocity
    vsv_h=vsv; #Horizontal velocity equal to total velocity
end

# FEM MODEL
# The model begins with the number of elements
# SOLID PART (the code can consider the water layer)
Nsolid = size(vpv,1); # number of elements in solid
hsolid =[diff(Mesh); Model[end,1]-Mesh[end,1]]*1000;
# Fluid Part
Nfluid = 0; # number of elements in fluid
hfluid = Step*ones(1,Nfluid); # grid spacing of mesh (meters)
vpfluid = 1500; rhofluid = 1030;
vpfv = vpfluid*ones(1,Nfluid);
rhofv = rhofluid*ones(1,Nfluid);

return fks,Mesh,vpv,vsv,vsv_v,vsv_h,rhov,qsv,qR,Nsolid,hsolid,Nfluid,hfluid,vpfluid,vpfv,rhofv
end

## FUNCTION TO FIX Outputs
" This function cleans a bit the outut of the 2 main functions"
##
function fix_SWF_output(Nmodes,cR, UR, evR, cL, UL, evL, fks)
# Make Periods vectors
T=1 ./fks; # Periods in s
# Organize the Velocity Dispersion Output
SWV=zeros(size(fks,1),4,Nmodes);
SWV[:,1,1:Nmodes]=cR[1:Nmodes,1,:]'#[cR[1,1,:] cR[2,1,:] cR[3,1,:]];  #safe Rayleigh phase velocity (m/s)
SWV[:,2,1:Nmodes]=UR[1:Nmodes,1,:]'#[UR[1,1,:] UR[2,1,:] UR[3,1,:]];  #safe Rayleigh group veloecity (m/s)
SWV[:,3,1:Nmodes]=cL[1:Nmodes,1,:]'#[cL[1,1,:] cL[2,1,:] cL[3,1,:]];  #safe Love phase velocity (m/s)
SWV[:,4,1:Nmodes]=UL[1:Nmodes,1,:]'#[UL[1,1,:] UL[2,1,:] UL[3,1,:]];  #safe Love group veloecity (m/s)
SWV=SWV ./1000; # Velocities in Km/s
# Get R wave Vertical Eigenfunction
EG_FUNC_R_Z=evR[2:2:end,:,1:end];
# Get R wave Radial Eigenfunction
EG_FUNC_R_R=evR[1:2:end,:,1:end];
# Get L wave Tangial Eigenfunction
EG_FUNC_L_T=evL
return T, SWV,EG_FUNC_R_Z, EG_FUNC_R_R, EG_FUNC_L_T
end

## Compute Stoneley waves velocities (Explenation at the end)
" Compute Stoneley waves velocities"
##
function Stoneley(a,b,c,f,s)
# Compute all the polynomial coefficients
c16 = 256*(b^16) - (512*(b^18))/(a^2) + (256*(b^20))/(a^4);
c14 = -768*(b^14) + (1280*(b^16))/(a^2) - (512*(b^18))/(a^4) - (512*(b^16))/(c^2) + (1024*(b^18))/((a^2)*(c^2)) - (512*(b^20))/((a^4)*(c^2));
c12 = 832*(b^12) - (1024*(b^14))/(a^2) + (256*(b^16))/(a^4) + (256*(b^16))/(c^4) - (512*(b^18))/((a^2)*(c^4)) + (256*(b^20))/((a^4)*(c^4)) + (1536*(b^14))/(c^2) - (2560*(b^16))/((a^2)*(c^2)) + (1024*(b^18))/((a^4)*(c^2)) - (64*(b^12)*(f^2))/(s^2);
c10 = -416*(b^10) + (288*(b^12))/(a^2) - (768*(b^14))/(c^4) + (1280*(b^16))/((a^2)*(c^4)) - (512*(b^18))/((a^4)*(c^4)) - (1664*(b^12))/(c^2) + (2048*(b^14))/((a^2)*(c^2)) - (512*(b^16))/((a^4)*(c^2)) + (96*(b^10)*(f^2))/(s^2) + (96*(b^12)*(f^2))/((a^2)*(s^2)) + (64*(b^12)*(f^2))/((c^2)*(s^2));
c8 = 112*(b^8) - (32*(b^10))/(a^2) + (832*(b^12))/(c^4) - (1024*(b^14))/((a^2)*(c^4)) + (256*(b^16))/((a^4)*(c^4)) + (832*(b^10))/(c^2) - (576*(b^12))/((a^2)*(c^2)) - (48*(b^8)*(f^2))/(s^2) - (128*(b^10)*(f^2))/((a^2)*(s^2)) - (32*(b^12)*(f^2))/((a^4)*(s^2)) - (96*(b^10)*(f^2))/((c^2)*(s^2)) - (96*(b^12)*(f^2))/((a^2)*(c^2)*(s^2));
c6 = -16*(b^6) - (416*(b^10))/(c^4) + (288*(b^12))/((a^2)*(c^4)) - (224*(b^8))/(c^2) + (64*(b^10))/((a^2)*(c^2)) + (16*(b^6)*(f^2))/(s^2) + (48*(b^8)*(f^2))/((a^2)*(s^2)) + (32*(b^10)*(f^2))/((a^4)*(s^2)) + (48*(b^8)*(f^2))/((c^2)*(s^2)) + (128*(b^10)*(f^2))/((a^2)*(c^2)*(s^2)) + (32*(b^12)*(f^2))/((a^4)*(c^2)*(s^2));
c4 = (b^4) + (112*(b^8))/(c^4) - (32*(b^10))/((a^2)*(c^4)) + (32*(b^6))/(c^2) + ((b^4)*(f^4))/(s^4) - (2*(b^4)*(f^2))/(s^2) - (16*(b^6)*(f^2))/((a^2)*(s^2)) - (16*(b^6)*(f^2))/((c^2)*(s^2)) - (48*(b^8)*(f^2))/((a^2)*(c^2)*(s^2)) - (32*(b^10)*(f^2))/((a^4)*(c^2)*(s^2));
c2 = -((16*(b^6))/(c^4)) - (2*(b^4))/(c^2) - (2*(b^4)*(f^4))/((a^2)*(s^4)) + (2*(b^4)*(f^2))/((a^2)*(s^2)) + (2*(b^4)*(f^2))/((c^2)*(s^2)) + (16*(b^6)*(f^2))/((a^2)*(c^2)*(s^2));
c0 = (b^4)/(c^4) + ((b^4)*(f^4))/((a^4)*(s^4)) - (2*(b^4)*(f^2))/((a^2)*(c^2)*(s^2));
# Find the roots of the polynomial and tirn to velocity
Candidates=sqrt.(1 ./real(PolynomialRoots.roots([c0,c2,c4,c6,c8,c10,c12,c14,c16])));
# Delete all the values larger than the velocity in the Fluid
deleteat!(Candidates, Candidates .> c)
# test which root best satisfies stoneley equation
seq=zeros(length(Candidates));
for ii=1:length(Candidates)
    v = Candidates[ii];
    seq[ii] = sqrt((b/v)^2 - (b/c)^2)*((1 - 2*(b/v)^2)^2 - (4*(b/v)^2)*sqrt((b/v)^2 - (b/a)^2)*sqrt((b/v)^2 - 1)) + (f/s)*sqrt((b/v)^2 - (b/a)^2);
end
#
V_Stoneley=Candidates[argmin(seq)]*1000;
return V_Stoneley
end
    # calculate the fundamental
    # mode velocity of the guided wave for a model of a halfspace of water over
    # a halfspace of an elastic solid. This is called a Stoneley wave since its
    # velocity is less than the water velocity (i.e. it is trapped in both
    # directions, up and down). In contrast, a Scholte wave occurs for a finite
    # water depth when the guided wave velocity is greater than in water. The
    # Stoneley wave velocity is the solution of an eighth order polynomial.

## LOVE DISPERSION
" THIS IS THE MAIN LOVE WAVE DISPERSION FUNCTION"
##
function Love_Forwardsp(Nsolid,vsv,rhov,f,hsolid,Nmodes,qsv,Earth::Int64=0,Elastic::Int64=1)
fnum=size(f,1);
# Spherical Correction:
    if Earth==1
        ER=6371;
        HH=[0; hsolid]
        SHELLS=[cumsum(HH[1:end-1]./1000) cumsum(HH[2:end]./1000)]
        RAD=ER .-SHELLS
        CORRV=(2*ER) ./(RAD[:,1]+RAD[:,2]);
        CORRRHO_L=CORRV.^-5;

        vsv=vsv.*CORRV
        rhov=rhov.*CORRRHO_L
    end
# Some moduli
muv = rhov.*vsv.*vsv;
# Angular frecuency
# Calculate make angular frequency
ω = 2*pi*f;
## Initialize some matrixes
# initialize some local matrices
L1 = spzeros(2,2);
L3 = spzeros(2,2);
M1 = spzeros(2,2);
# initialize the global matrix
Ka1 = spzeros(Nsolid,Nsolid);
Ka3 = spzeros(Nsolid,Nsolid);
M = spzeros(Nsolid,Nsolid);
# For Solid Part of the Model
# Loop for the solid part of the Model
for ii=1:Nsolid
    # grab grid interval of current element
    h = hsolid[ii];
    #grab material properties of current element
    mu = muv[ii];
    # make elemental mass matrix
    M1 = spzeros(2,2);
    M1[1,1] = h*rhov[ii]/2;
    M1[2,2] = h*rhov[ii]/2;
    # make elemental stiffness matrices
    L1 = spzeros(2,2);
    L3 = spzeros(2,2);
    # some alternate variables from Lysmer
    alph = mu/6;
    bet = mu/6;
    # the 16 entries of the 4x4 elemental stiffness matrices of Lysmer
    L1[1,1] = 2*alph*h;
    L3[1,1] = (6*bet/h);
    L1[1,2] = alph*h;
    L3[1,2] = -(6*bet/h);
    L1[2,1] = L1[1,2];
    L3[2,1] = L3[1,2];
    L1[2,2] = L1[1,1];
    L3[2,2] = L3[1,1];
    # assemble mass and stiffness matrices from elemental matrices
    if (ii == Nsolid)
        M[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] = M[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] + M1[1:1,1:1];
        Ka1[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] = Ka1[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] + L1[1:1,1:1];
        Ka3[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] = Ka3[(1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)] + L3[1:1,1:1];
    else
        M[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = M[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + M1;
        Ka1[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = Ka1[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + L1;
        Ka3[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = Ka3[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + L3;
    end
end
# Set a lower bound on the Love wave speed
lspd = minimum(vsv);
# find the eigenvalue closest to the upper-bound eigenvalue
x=zeros(size(M,1),Nmodes,fnum);
d=zeros(Nmodes,1,fnum);
for o=1:fnum
    evm1=(ω[o]*ω[o]*M)-Ka3;
    evm2=Ka1;
    @fastmath @inbounds dp,xp=eigs(evm1,evm2,nev=Nmodes,sigma=(ω[o]/lspd)^2,tol=0.0, ncv=20, maxiter=300, which=:LM);
    @fastmath @inbounds dp=real.(dp);
    @fastmath @inbounds xp=real.(xp)
    # Test for artificial numerical modes and remove them
    pmi=zeros(Int64,Nmodes,1);
    for mm=1:Nmodes
        # test if artificial or not
        if (maximum(diff(xp[:,mm])/maximum(abs.(xp[:,mm]))) > 0.5)
            pmi[mm] = 1;
        else
            pmi[mm]= 0;
        end
    end
    if sum(pmi)> 0
        println("LW at $o Some artificial Modes have been found, please densify the grid and repeat the computation")
    end
    x[:,:,o] = xp;
    d[:,:,o] = dp;
end
# Make sure none positive are not accounted for:
d[findall(d.<0)] .= NaN
# Normalize the eigenfunction
ev=zeros(Nsolid,Nmodes,fnum)
for o=1:fnum
iev=zeros(Nsolid,Nmodes)
fctr=zeros(Nmodes,1)
    for i=1:Nmodes
        N=1 ./(x[1:Nsolid,i,o]'*M*x[1:Nsolid,i,o])
        fctr[i] = N[1];
        iev[:,i] = x[1:Nsolid,i,o].*sqrt.(fctr[i]);
    end
    ev[:,:,o]=iev
end
# NOW! Some results!
# The wavenumber used is d
# The  Computed phase velocity (c)
vpk=zeros(Nmodes,1,fnum);
for o=1:fnum
    vpk[:,:,o] = ω[o]./sqrt.(d[:,:,o]);
end
# The Computed group velocity (U)
vgk=zeros(Nmodes,1,fnum);
for o=1:fnum
    ivgk=zeros(size(vpk,1),1);
    for i=1:size(vpk,1)
        ivgk[i] = (transpose(x[1:Nsolid,i,o])*(2*sqrt.(d[i,:,o]).*Ka1)*x[1:Nsolid,i,o]/(2*ω[o]))*(1/((transpose(x[1:Nsolid,i,o]))*M*x[1:Nsolid,i,o]));
    end
    vgk[:,:,o]=ivgk
end
# Atenuation Aproximation
# Working with the eigenfunctions
if Elastic == 0
for o=1:fnum
# The eigenfunctions
V3=ev[:,:,o] # SH-wave eigenfunction
# Components of the Lagrangian
I0=sum(rhov .*(V3.*V3) .*hsolid,dims=1);
#δcδvs
A=(vsv.*rhov)/(I0'.*vgk[:,:,o]);
TERM2=((([diff(V3,dims=1); zeros(1,Nmodes)])./hsolid).*(1 ./sqrt.(d[:,:,o]))').^2;
TERM1=(V3).^2;
B=sum((TERM1+TERM2).*hsolid,dims=1);
δcδvs=A.*B;
#δcδvp: Love wave do not depend on Vp velocity
δcδvp=0;
# Attenuation Term
Qp=Inf;
#Qs=85;
Qterm=sum((δcδvs.*(vsv./1000)).*(1 ./qsv),dims=1)
# 𝛄 gamma
TREF=50; # Reference period 50 s
ωr=2*pi/TREF # T= 50 s as a reference
𝛄=(ω[o] ./(2*(vpk[:,:,o].*vpk[:,:,o]))).*Qterm';
# Anaelastic Phase velocity
c=vpk[:,:,o].+((1/pi*log(ω[o]/ωr))*Qterm)';
# Anaelastic Group velocity
A=2 .-(vgk[:,:,o]./vpk[:,:,o]);
B=(c-vpk[:,:,o])./vpk[:,:,o];
C=2*𝛄.*vgk[:,:,o]/pi*ω[o];
U=vgk[:,:,o].*(1 .+(A.*B) .+C);
# Safe the Values
vpk[:,:,o]=c;
vgk[:,:,o]=U;
end
end
return d, vpk, vgk, ev
end

## RAYLEIGH DISPERSION
" THIS IS THE MAIN RAYLEIGH WAVE DISPERSION FUNCTION"
##
function Rayleigh_Forwardsp(Nsolid,vsv,vpv,rhov,f,hsolid,Nmodes,Nfluid,vpfv,rhofv,hfv,qsv,Earth::Int64=0,Elastic::Int64=1)
fnum=size(f,1);
# Spherical Correction:
if Earth==1
    ER=6371;
    HH=[0; hsolid]
    SHELLS=[cumsum(HH[1:end-1]./1000) cumsum(HH[2:end]./1000)]
    RAD=ER .-SHELLS
    CORRV=(2*ER) ./(RAD[:,1]+RAD[:,2]);
    CORRRHO_R=CORRV.^-2.275;
    vpv=vpv.*CORRV
    vsv=vsv.*CORRV
    rhov=rhov.*CORRRHO_R
end
# Fluid part of the model
# Check for the number of nodes in the fluid, based on the number of elements
if (Nfluid > 0)
    Nnfo = Nfluid + 1;
else
    Nnfo = 0;
end
# make fluid portion of model
# make kappa of the fluid , the modulus
kappafv = rhofv.*vpfv.*vpfv;
# Angular frecuency
# Calculate make angular frequency
ω = 2*pi*f;
# Initialize some matrixes
# initialize some local matrices
L1 = spzeros(2,2);
L3 = spzeros(2,2);
M1 = spzeros(2,2);
# initialize the global matrix
Ka1 = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
Ka3 = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
M = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
# MATRIXES IN THE FLUID PART OF THE Model
if Nfluid>0
# for all elements
    for ii=1:Nfluid
        # grab grid interval of current element
        h = hfluid[ii];
        # grab material properties of current element
        rhof = rhofv[ii];
        kappaf = kappafv[ii];
        # make elemental mass matrix
        M1 = spzeros(2,2);
        M1[1,1] = h/(2*kappaf);
        M1[2,2] = h/(2*kappaf);
        # make elemental stiffness matrices
        L1 = spzeros(2,2);
        L3 = spzeros(2,2);
        # some alternate variables from Lysmer
        alph = 1/(6*rhof);
        bet = 1/(6*rhof);
        # the 4 entries of the 2x2 elemental stiffness matrices of Lysmer
        L1[1,1] = 2*alph*h;
        L3[1,1] = (6*bet/h);
        L1[1,2] = alph*h;
        L3[1,2] = -(6*bet/h);
        L1[2,1] = L1[1,2];
        L3[2,1] = L3[1,2];
        L1[2,2] = L1[1,1];
        L3[2,2] = L3[1,1];
        # assemble mass and stiffness matrices from elemental matrices
        M[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = M[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + M1;
        Ka1[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = Ka1[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + L1;
        Ka3[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] = Ka3[(1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))] + L3;

    end
end
# Correct before solid part of Model
M[1,1] = M[1,1]*2;
Ka1[1,1] = Ka1[1,1]*2;
Ka3[1,1] = Ka3[1,1]*2;
# For Solid Part of the Model
# make solid portion of model
# make mu and lambda
muv = rhov.*vsv.*vsv;
lamdav = rhov.*vpv.*vpv - 2*muv;
#initialize some matrices
Ka2 = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
L1 = spzeros(4,4);
L2 = spzeros(4,4);
L3 = spzeros(4,4);
M1 = spzeros(4,4);
# Loop for the solid part of the Model
 for ii=1:Nsolid
    # grab grid interval of current element
    h = hsolid[ii];
    #grab material properties of current element
    mu = muv[ii];
    lamda = lamdav[ii];
    # make elemental mass matrix
    M1 = spzeros(4,4);
    M1[1,1] = h*rhov[ii]/2;
    M1[2,2] = h*rhov[ii]/2;
    M1[3,3] = h*rhov[ii]/2;
    M1[4,4] = h*rhov[ii]/2;
    # make elemental stiffness matrices
    L1 = spzeros(4,4);
    L2 = spzeros(4,4);
    L3 = spzeros(4,4);
    # some alternate variables from Lysmer
    alph = ((2*mu)+lamda)/6;
    bet = mu/6;
    theta = (mu+lamda)/4;
    psi = (mu-lamda)/4;
    # the 16 entries of the 4x4 elemental stiffness matrices of Lysmer
    L1[1,1] = 2*alph*h;
    L3[1,1] = (6*bet/h);
    L2[1,2] = 2*psi;
    L1[1,3] = alph*h;
    L3[1,3] = -(6*bet/h);
    L2[1,4] = 2*theta;
    L2[2,1] = L2[1,2];
    L1[2,2] = 2*bet*h;
    L3[2,2] = (6*alph/h);
    L2[2,3] = -2*theta;
    L1[2,4] = bet*h;
    L3[2,4] = -(6*alph/h);
    L1[3,1] = L1[1,3];
    L3[3,1] = L3[1,3];
    L2[3,2] = L2[2,3];
    L1[3,3] = L1[1,1];
    L3[3,3] = L3[1,1];
    L2[3,4] = -2*psi;
    L2[4,1] = L2[1,4];
    L1[4,2] = L1[2,4];
    L3[4,2] = L3[2,4];
    L2[4,3] = L2[3,4];
    L1[4,4] = L1[2,2];
    L3[4,4] = L3[2,2];
    # assemble mass and stiffness matrices from elemental matrices
    if ii == Nsolid
        M[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] = M[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] + M1[1:2,1:2];
        Ka1[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] = Ka1[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] + L1[1:2,1:2];
        Ka2[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] = Ka2[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] + L2[1:2,1:2];
        Ka3[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] = Ka3[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))] + L3[1:2,1:2];
    else
        M[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] = M[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] + M1;
        Ka1[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] = Ka1[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] + L1;
        Ka2[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] = Ka2[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] + L2;
        Ka3[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] = Ka3[(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))] + L3;
    end
end
# Construct the coupling matrix
if (Nfluid > 0)
Cm = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
Cm[Nnfo,Nnfo+2] = 1;
Cm[Nnfo+2,Nnfo] = 1;
else
    Cm = spzeros(Nnfo+(2*Nsolid),Nnfo+(2*Nsolid));
end
# find the rayleigh/scholte wave speed which would exist if the solid model
# were a halfspace with the minimum model velocity
# this is a lower bound on the velocity that can be passed to EIGS, based
# on ARPACK
if (Nfluid > 0)
    msval  = minimum(vsv);
    msloc = argmin(vsv);
    mfval  = minimum(vpfv);
    mfloc = argmin(vpfv);
    vsmay = msval;
    vpmay = vpv[msloc];
    vpfmay = mfval;
    rhofmay = rhofv[mfloc];
    rhomay = rhov[msloc];
    rspd = Stoneley(vpmay/1000,vsmay/1000,vpfmay/1000,rhofmay/1000,rhomay/1000);
else
    msval  = minimum(vsv);
    msloc = argmin(vsv);
    vsmay = msval;
    vpmay = vpv[msloc];
    # coefficients of rayleigh's polynomial
    t1 = 1/(vsmay^6);
    t2 = -8/(vsmay^4);
    t3 = ((24/(vsmay^2))-(16/(vpmay^2)));
    t4 = -16*(1-((vsmay/vpmay)^2));
    # rayleigh wave speed
    rspd = sqrt(minimum(real(PolynomialRoots.roots([t4, t3, t2, t1]))));
end
# Find the eigenvalue closest to the upper-bound eigenvalue
x=zeros(2*size(M,1),Nmodes,fnum);
d=zeros(Nmodes,1,fnum);
for o=1:fnum
    evm1=[spzeros((Nnfo+(2*Nsolid)),(Nnfo+(2*Nsolid))) sparse(I,(Nnfo+(2*Nsolid)),(Nnfo+(2*Nsolid))); ((ω[o]*ω[o]*M)-Ka3-(ω[o]*Cm)) Ka2];
    evm2=[sparse(I,(Nnfo+(2*Nsolid)),(Nnfo+(2*Nsolid))) spzeros((Nnfo+(2*Nsolid)),(Nnfo+(2*Nsolid))); spzeros((Nnfo+(2*Nsolid)),(Nnfo+(2*Nsolid))) Ka1];
    @fastmath @inbounds dp,xp=eigs(evm1,evm2,nev=Nmodes,sigma=ω[o]/rspd,tol=0.0, ncv=20, maxiter=300, which=:LM);
    @fastmath @inbounds x[:,:,o] = real.(xp);
    @fastmath @inbounds d[:,:,o] = real.(dp);
end
# Normalize the eigenfunctions
ev=zeros(size((Nnfo+1):(Nnfo+(2*Nsolid)),1),Nmodes,fnum);
for o=1:fnum
    iev=zeros(size((Nnfo+1):(Nnfo+(2*Nsolid)),1),Nmodes);
    for e=1:Nmodes
        fctr = 1/ (((x[1:1:(Nnfo+(2*Nsolid)),e,o]')*M*x[1:1:(Nnfo+(2*Nsolid)),e,o])-(((x[1:1:(Nnfo+(2*Nsolid)),e,o]')*Cm*x[1:1:(Nnfo+(2*Nsolid)),e,o])/2*ω[o]));
        evp = x[1:1:(Nnfo+(2*Nsolid)),e,o]*sqrt(fctr)*sign(x[Nnfo+1,e,o]);
        # # return only the eigenvector in the solid
        iev[:,e] = evp[(Nnfo+1):(Nnfo+(2*Nsolid))];
    end
    ev[:,:,o]=iev
end
# NOW! Some results!
# Make sure none positive are not accounted for:
d[findall(d.<0)] .= NaN
# The wavenumber used is d
# The  Computed phase velocity (c)
vpk=zeros(Nmodes,1,fnum);
for o=1:fnum
    vpk[:,:,o] = ω[o]./d[:,:,o];
end
# The Computed group velocity (U)
vgk=zeros(Nmodes,1,fnum);
for o=1:fnum
    ivgk=zeros(size(vpk,1),1);
    for i=1:size(vpk,1)
        ivgk[i] = (((x[1:1:(Nnfo+(2*Nsolid)),i,o]')*((2*d[i,1,o]*Ka1)-Ka2)*x[1:1:(Nnfo+(2*Nsolid)),i,o])/((2*ω[o]*((x[1:1:(Nnfo+(2*Nsolid)),i,o]')*M*x[1:1:(Nnfo+(2*Nsolid)),i,o]))-((x[1:1:(Nnfo+(2*Nsolid)),i,o]')*Cm*x[1:1:(Nnfo+(2*Nsolid)),i,o])));
    end
    vgk[:,:,o]=ivgk
end
# Atenuation Aproximation
# Working with the eigenfunctions
if Elastic == 0
for o=1:fnum
    # The eigenfunctions
V1=ev[2:2:end,:,o] # vertical Eigenfcuntions
V2=ev[1:2:end,:,o] # radial eigenfunctions
# Components of the Lagrangian
I0=sum(rhov.*(V1.*V1 .+ V2.*V2).*hsolid,dims=1);
#δc/δvp
A=(vpv.*rhov)/(I0'.*vgk[:,:,o]);
derivativeV1=(([diff(V1,dims=1); zeros(1,Nmodes)])./hsolid) .*(1 ./d[:,:,o])'
B=sum(((V2-derivativeV1).^2).*hsolid,dims=1)
δcδvp=A.*B;
#δc/δvs
A=(vsv.*rhov)/(I0'.*vgk[:,:,o]);
derivativeV1=(([diff(V1,dims=1); zeros(1,Nmodes)])./hsolid) .*(4 ./d[:,:,o])'
derivativeV2=(([diff(V2,dims=1); zeros(1,Nmodes)])./hsolid) .*(1 ./d[:,:,o])'
TERM1=(V1+derivativeV2).^2
TERM2=V2.*derivativeV1
B=sum((TERM1+TERM2).*hsolid,dims=1)
δcδvs=A.*B;
# Attenuation Term
Qp=Inf;
Qterm=sum((δcδvp.*(vpv./1000)).*(1/Qp)+(δcδvs.*(vsv./1000)).*(1 ./qsv),dims=1)
# 𝛄 gamma
TREF=50; # Reference period 50 s
ωr=2*pi/TREF # T= 50 s as a reference
𝛄=(ω[o] ./(2*(vpk[:,:,o].*vpk[:,:,o]))).*Qterm';
# Anaelastic Phase velocity
c=vpk[:,:,o]+((1/pi*log(ω[o]/ωr))*Qterm)'
# Anaelastic Group velocity
A=2 .-(vgk[:,:,o]./vpk[:,:,o]);
B=(c-vpk[:,:,o])./vpk[:,:,o];
C=2*𝛄.*vgk[:,:,o]/pi*ω[o];
U=vgk[:,:,o].*(1 .+(A.*B) .+C);
# Safe the Values
vpk[:,:,o]=c;
vgk[:,:,o]=U;
end
end
# Return Values
return d, vpk, vgk, ev
end

## SHOWME
"shome me is a simple function to see a full vector in Julia Repel"
##
function showme(A)
show(IOContext(stdout, :limit=>false), MIME"text/plain"(), A)
return
end

## Make_1D_Mesh


function Make_1D_Mesh(Tmin,DepthsM=nothing)
## Inputs
#Tmin=5; # Read min T from file
#Depths=[0 5 20 40 200  2889]; # Read Depths from the user input()

## Constants
if DepthsM!=nothing
    zmax=DepthsM[end]; # Read CMB from file
else
    zmax=2890; # Read CMB from file
end

#if Tmin >= 10
    vmin=1; # min velocity of the R wave c
    n=40; #According to Haney
    a=0.5; #or 0.63
    nodes= 5; # Number of nodes for each frequency
# else
#     vmin=0.1; # min velocity of the R wave c
#     n=40; #According tho Haney
#     a=0.5; #or 0.63
#     nodes= 30; # Number of nodes for each frequency
# end

## Compute some values & get the initial Depths of the mesh()
zmin=(vmin*Tmin*a)/nodes; # Zmin is the smalles depth of a node in the model for the required computation

# if zmin < 0.01
#     zmin=0.01;
# end

Nmax=ceil(n*a*log(zmax/zmin)); # Max bumber of layers
@fastmath @inbounds Z=zmin.*exp.((0:1:Nmax)./(n*a)); # Deoth vector that grows exponentially

## Polish the mesh a little bit
Z[end]=zmax; # Make the end of the vector the same as the input of the user

# Make sure that the User's Input depths appear in the mesh by moving the
# closest point to the user's value
if DepthsM!=nothing
    for i=2:size(Depths,2)-1
        Ind = argmin((abs.(Z.-DepthsM[i]/1000)))
        Z[Ind]=DepthsM[i]/1000;
    end
end


## From Depths to distance between nodes
FEM=[zmin; diff(Z,dims=1)]; # Maje the FEM mesh()

## Delete really small values & replace with the zmin
nfix=ceil(Int,sum(FEM[FEM.<zmin])/zmin);
FEM_fixed=[Z[1]; ones(nfix)*zmin; FEM[FEM.>zmin]]
return FEM_fixed

end

## GOLDEN Mesh
"I call this the golden mesh. It is the perfect FEM mesh for the problem according to Haney and Tsai
it is rounded uo to be able to be used and covers the entire mantle"
##
GOLD_MESH=[


0.5
  0.5
  0.5
  0.5
  0.5
  0.5
  0.5
  0.5
  0.5
  0.5
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.0
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  1.5
  2.0
  2.0
  2.0
  2.0
  2.0
  2.0
  2.0
  2.0
  2.0
  2.0
  2.0
  2.0
  2.0
  2.0
  2.5
  2.5
  2.5
  2.5
  2.5
  2.5
  2.5
  2.5
  2.5
  2.5
  2.5
  2.6
  2.6
  2.7
  2.7
  2.8
  2.9
  2.9
  3.0
  3.0
  3.1
  3.2
  3.2
  3.3
  3.4
  3.4
  3.5
  3.6
  3.6
  3.7
  3.8
  3.9
  3.9
  4.0
  4.1
  4.2
  4.3
  4.3
  4.4
  4.5
  4.6
  4.7
  4.8
  4.9
  5.0
  5.1
  5.2
  5.3
  5.4
  5.5
  5.6
  5.7
  5.9
  6.0
  6.1
  6.2
  6.4
  6.5
  6.6
  6.7
  6.9
  7.0
  7.2
  7.3
  7.5
  7.6
  7.8
  7.9
  8.1
  8.2
  8.4
  8.6
  8.7
  8.9
  9.1
  9.3
  9.5
  9.7
  9.9
 10.1
 10.3
 10.5
 10.7
 10.9
 11.1
 11.3
 11.6
 11.8
 12.0
 12.3
 12.5
 12.8
 13.0
 13.3
 13.6
 13.8
 14.1
 14.4
 14.7
 15.0
 15.3
 15.6
 15.9
 16.2
 16.6
 16.9
 17.2
 17.6
 18.0
 18.3
 18.7
 19.1
 19.4
 19.8
 20.2
 20.6
 21.1
 21.5
 21.9
 22.4
 22.8
 23.3
 23.7
 24.2
 24.7
 25.2
 25.7
 26.2
 26.8
 27.3
 27.9
 28.4
 29.0
 29.6
 30.2
 30.8
 31.4
 32.0
 32.7
 33.3
 34.0
 34.7
 35.4
 36.1
 36.9
 37.6
 38.4
 39.1
 39.9
 40.7
 41.5
 42.4
 43.2
 44.1
 45.0
 45.9
 46.8
 47.8
 48.7
 49.7
 50.7
 51.8
 52.8
 53.9
 55.0
 56.1
 57.2
 47.7];
