#By Mariano Arnaiz
"This file is a Julia interpretation of the Matlab file rf_forward2.m
Computation of receiver functions and displacements at the surface
for plane P-wave impingent from below.
 The stack of homogeneous, isotropic layers is completed by elastic air
 rho = 1.2 kg/m3
 VP  = 0.34 km/s
 VS  = 0.02 km/s i.e. VERY low shear modulus
INPUTS:
 zs [km]     : Surface at z=0,
               zs[1]   = Bottom of first layer
               zs[end] = Top of halfspace
 VPs [km(s]  : VPs[1] between 0 and zs(1)
               VPs[i] between zs(i-1) and zs(i)
               VPs[end] below zs(end)
               hence length(VPs) = length(zs)+1
 VSs [km/s]  : default VSs = VPs/sqrt(3)
 rhos [kg/m3]: default rhos = 0.32*(VPs*1000) + 770; %Birch medium
 p [s/deg]   : wave slowness, default 6.4 s/deg
 ts [s]      : Time axis, assumed periodic. Take care of wrap around. Note
               that it is internally extended to a nice length for fft.
 Gaussian_factor : With omegas = 2*pi*frequency the waveform spectrum is
               H(omegas) = k*exp(-omegas.^2/(4*Gaussian_factor^2))
 rot_ang [deg]: Rotation of (r,up) before RF-computation.
               If rot_ang is NaN or not specified then default is the
               theoretical angle, i_app, where P-energy is cancelled on Q.
 OUTPUTS:
% Q_rf : Radial reciever function
% L_rf : Receiver function for vertical component (positive upwards)
% ts_nice: time samples, modified to get efficient fft."

using LinearAlgebra, StructArrays, FFTW, Statistics
include("Get_p.jl")

function P_RF_forward(Model,depth,dist,Gaussian_factor,Clean_Model::Int64=1,rotate::Int64=0)

## Compute surface response in Fourier domain
# for plane impingent delta P-wave wavefront
# Transmitted P-wave impulse arrives at surface at ts=0

## ts: we will allways compute the entire receiver function!
ts = -10:0.05:90;
#ts = -120:0.05:120;
##Compute ray parameter p
p=Get_p(depth,dist,"P")
## Clean model 1 run the clean mode, 0 use the model as it is (May cause problems)
if Clean_Model ==1
    zs,VPs,VSs,rhos=Fix_RF_Model(Model)
    println("Cleaning model")
    Used_model=[[zs Model[end,1]*1000]' VPs' VSs' rhos'];
else
Model=Model[searchsortednearest(Model[:,1],0.5):searchsortednearest(Model[:,1],1000),:];
zs=Model[2:end,1];
VPs=Model[:,2];
VSs=Model[:,3];
rhos=Model[:,4];
VPs = 1000*VPs';
VSs = 1000*VSs';
zs = 1000*zs';
rhos = 1000*rhos';
zs = [0 zs]; #soil at sea level implicit; force a row
air_fac = 0.1; #air factor; 1=sea level; 0.1=stratosphere# Density is assumed to obey Birch's-relation:
VPs  = [340*air_fac VPs];
rhos = [1.2*air_fac rhos];
VSs  = [(0.1*VPs[1]/sqrt(3)) VSs]; #%Slightly unphysical elastic air
Used_model=[[zs Model[end,1]*1000]' VPs' VSs' rhos'];
#println("Using Model as it is")
end

## rotation
if rotate == 0
    rot_ang = 0; #theoretical apparent incidence angle
    #println("No rotation")
else
    rot_ang = 2*asind((VSs[2]/1000)*(p/111.111)); #theoretical apparent incidence angle
    #println("Considering Rotation, rot_ang = $rot_ang ")
end

## Other parameters
p_sm = p/111111; #1 degree per 111111 meter
ns = 2;
N_layers = length(VPs)-2; #the number of finite layers, typically layers inside the crust
ts = ts';
##times
dt=ts[2]-ts[1];
Nt = size(ts,2);
T_max = ts[Nt];
N2=nextpow(2,Nt);
N3=round(N2/1.5);
p2 = 2 .^(1:N2);#we demand nice_N to be even
p3 = 3 .^(0:N3);
p22,p33=meshgrid(p2,p3);
p2233 = p22.*p33;
p2233sort = sort(p2233[:]);
N_fft = Int64(minimum(p2233sort[p2233sort .>= Nt]));
ts_fft = T_max .- dt .*((N_fft-1):-1:0); #pad before zero
T_min_fft = ts_fft[1];
df = 1/(N_fft*dt);
fs = df .*[0:N_fft/2-1;-N_fft/2:-1];
#
# ##fixed times
# dt=ts[2]-ts[1];
# Nt = 2001;
# T_max = 90.0;
# N2=2048;
# N3=1365.0;
# # p2 = 2 .^(1:N2);#we demand nice_N to be even
# # p3 = 3 .^(0:N3);
# # p22,p33=meshgrid(p2,p3);
# # p2233 = p22.*p33;
# # p2233sort = sort(p2233[:]);
# N_fft = 2048#Int64(minimum(p2233sort[p2233sort .>= Nt]));
# ts_fft = T_max .- dt .*((N_fft-1):-1:0); #pad before zero
# T_min_fft = ts_fft[1];
# df = 1/(N_fft*dt);
# fs = df .*[0:N_fft/2-1;-N_fft/2:-1];

##All Wave conversions
Td,Ru,Tu,Rd,S=vrho2tr(p_sm,VPs,VSs,rhos);
M=zeros(4,4,N_layers+1);
#Theb do
for i=1:N_layers+1
    invTd=inv(Td[:,:,i]);
    M[:,:,i] = [invTd -invTd*Ru[:,:,i];Rd[:,:,i]*invTd (Tu[:,:,i]-Rd[:,:,i]*invTd*Ru[:,:,i])];
end

## %Wave inputs:
Ad = zeros(2,N_fft); #downgoing wavefields above earth surface; no input from above.
#Bu = [ones(1,N_fft);zeros(1,N_fft)]; %upgoing wavefields below lowermost interface; impulsive P-wave at t=T_min_fft and no S-wave
P_shift = sum(sqrt.((1 ./VPs[2:N_layers+1].^2) .-p_sm^2) .*diff(zs',dims=1)); #apparent vertical traveltime of direct P-wave from bottom to surface.
Bu = [shift(fs,-T_min_fft - P_shift)' ; zeros(1,N_fft)]; #upgoing wavefields below lowermost interface; impulsive P-wave at t=T_min_fft and no S-wave
#Wave responses (just prealocations:
global Au = 9999*im*ones(2,N_fft); #upgoing wavefields above earth surface
global Bd = 9999*im*ones(2,N_fft); #downgoing wavefields below lowermost interface
global B_intern = zeros(ComplexF64,N_fft,4,length(ns)); #accumulate Fouriertransforms of wavefields beneath indexed interfaces.

global C = zeros(ComplexF64,4,4,N_layers+1); #Tdelay effects, one for each homogeneous region
global C[1,1,1] = im; #just to force the array size to be complex
##
#Can we make do?  med N_fft/2+1???
global FBCM = zeros(ComplexF64,4,4,N_layers+2);
for j=1:Int64((N_fft/2)+1)

  global C[1,1,2:N_layers+1] = exp.(+im*2*pi*fs[j]*sqrt.((1 ./VPs[2:N_layers+1].^2) .-p_sm^2) .*diff(zs',dims=1))';
  global C[2,2,2:N_layers+1] = exp.(+im*2*pi*fs[j]*sqrt.((1 ./VSs[2:N_layers+1].^2) .-p_sm^2) .*diff(zs',dims=1))';
  global C[3,3,:] = conj(C[1,1,:]);
  global C[4,4,:] = conj(C[2,2,:]);
  global FBCM = zeros(ComplexF64,4,4,N_layers+2);
  global FBCM[:,:,N_layers+2]=I(4);
  for i=(N_layers+1):-1:2
      global FBCM[:,:,i]=C[:,:,i]*M[:,:,i]*FBCM[:,:,i+1];
  end
  global FBCM[:,:,1]=M[:,:,1]*FBCM[:,:,2];
  D = FBCM[:,:,1];
  Ddd = D[1:2,1:2];
  Dud = D[1:2,3:4];
  Ddu = D[3:4,1:2];
  Duu = D[3:4,3:4];

  invDdd = inv(Ddd);
  global Bd[:,j] = Ddd\(Ad[:,j]-Dud*Bu[:,j]);
  global Au[:,j] = Ddu*Bd[:,j] + Duu*Bu[:,j];
  global B_bottom = [Bd[:,j];Bu[:,j]];
  for inn=1:length(ns)
     global B_intern[j,:,inn] = FBCM[:,:,ns[inn]]*B_bottom;
  end
end

# Now consider the remaining negative frequencies:
for k=Int64(N_fft/2+2):N_fft
    j=Int64(N_fft-k+2);
    global Bd[:,k] = conj(Bd[:,j]);
    global Au[:,k] = conj(Au[:,j]);
    for in=1:length(ns)
         global B_intern[k,:,in] = conj(B_intern[j,:,in]);
    end
end


Pu = Au[1,:]; #P-wave in elastic air
Su = Au[2,:]; #S-wave in elastic air

##
#--->>> BHJ change May 2008:
#VPs(in) changed to VPs(ns(in)) etc.
#effect is that previously we converted P and S below surface with
# a PS2rz computed in air. The effect was small; a few percent!
#This error is not present in code used for inversion, hence results in
#BSSA-paper "Enhanced ..." are valid.

Bpd2rz=zeros(2,length(ns));
Bsd2rz=zeros(2,length(ns));
Bpu2rz=zeros(2,length(ns));
Bsu2rz=zeros(2,length(ns));
PS2rz=zeros(ns,4,length(ns));
for inn=1:length(ns)
   Bpd2rz[:,inn] = [VPs[ns[inn]]*p_sm;+sqrt(1-(VPs[ns[inn]]*p_sm)^2)];
   Bsd2rz[:,inn] = [+sqrt(1-(VSs[ns[inn]]*p_sm)^2);-VSs[ns[inn]]*p_sm];
   Bpu2rz[:,inn] = [VPs[ns[inn]]*p_sm;-sqrt(1-(VPs[ns[inn]]*p_sm)^2)];
   Bsu2rz[:,inn] = [+sqrt(1-(VSs[ns[inn]]*p_sm)^2);VSs[ns[inn]]*p_sm];
   PS2rz[:,:,inn] =  [Bpd2rz[:,inn] Bsd2rz[:,inn] Bpu2rz[:,inn] Bsu2rz[:,inn]];
end
##
#Compute displacements R(fs) and UP(fs) beneath surface and associated receiver function RF(fs)
inn = 1; #ns(1) should be 2 in order to make values relate to the position just below the surface
R_pd =  Bpd2rz[1,inn].*B_intern[:,1,inn];
UP_pd = -Bpd2rz[2,inn].*B_intern[:,1,inn];

R_sd =  Bsd2rz[1,inn].*B_intern[:,2,inn];
UP_sd = -Bsd2rz[2,inn].*B_intern[:,2,inn];

R_pu =  Bpu2rz[1,inn].*B_intern[:,3,inn];
UP_pu = -Bpu2rz[2,inn].*B_intern[:,3,inn];

R_su =  Bsu2rz[1,inn].*B_intern[:,4,inn];
UP_su = -Bsu2rz[2,inn].*B_intern[:,4,inn];
# Get R and Z components
R  = R_pd + R_sd + R_pu + R_su;
UP = UP_pd + UP_sd + UP_pu + UP_su;

# Begin the Lump array, this array holds a lot of useful and silly information
# we do not save as an struct array to avoid complications
lump_RF_fd = R./UP;
lump_R_fd = R;
lump_UP_fd = UP;

lump_FBCM = FBCM;
lump_B_bottom = B_bottom;
lump_B_intern = B_intern;
lump_PS2rz = PS2rz;
# # Transform back to timedomain:
 R = real(ifft(R));
 UP = real(ifft(UP));
##

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Compute the time domain displacements
#% as well as the receiver function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omegas = 2*pi*fs[:];
H = exp.(-omegas.^2/(4*Gaussian_factor^2));
h = real(ifft(H));
k = 1/maximum(h);
H = k*H; #normalized to amplitude 1.0 in time domain

# H2 = exp.(-omegas.^2/(4*(0.5*Gaussian_factor)^2));
# h2 = real(ifft(H2));
# k2 = 1/maximum(h2);
# H2 = k2*H2; #normalized to amplitude 1.0 in time domain


r = real(ifft(H.*lump_R_fd));
up = real(ifft(H.*lump_UP_fd));

shift_fd = shift(fs[:],-minimum(ts_fft));
r_rf = real(ifft((H.*lump_R_fd)./(lump_UP_fd.*shift_fd))); #note that the spectrum of UP normally has no zero crossings
up_rf = real(ifft(H.*1.0.*shift_fd)); #%UP./UP = 1.0

#now do the rotation
L_fd = ( cosd(rot_ang)*lump_UP_fd + sind(rot_ang)*lump_R_fd);#.*shift_fd;
Q_fd = (-sind(rot_ang)*lump_UP_fd + cosd(rot_ang)*lump_R_fd);#.*shift_fd;

L = real(ifft(H.*L_fd));
Q = real(ifft(H.*Q_fd));
L_rf = real(ifft(H.*1.0.*shift_fd));
Q_rf = real(ifft((H.*Q_fd)./(L_fd.*shift_fd)));

#safe more stuff
lump_H = H;
lump_h = h;
lump_r = r;
lump_up = up;
lump_r_rf = r_rf;
lump_up_rf = up_rf;
lump_rot_ang = rot_ang;

# Requested outpus
R_rf=reverse(r_rf); #Unrotated Radial Receiver Function
Z_rf=up_rf; #Unrotated Receiver function for vertical component (positive upwards)
Q_rf = reverse(Q_rf); #Rotated Radial Receiver Function
L_rf = L_rf; # Rotated Receiver function for vertical component (positive upwards)
time = ts_fft; #time samples, modified to get efficient fft.

# return time, R_rf, Z_rf, Q_rf, L_rf, Used_model
return time, Q_rf, Used_model

end

##
"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file = vrho2tr.jl
% Derived from vrho2tr.m that is derived from ps_conv.m (update: 17/12-2002; after AkiRich.m 5/11-02 BHJ & Sonny Krag)
% Purpose: returns conversioncoefficients at several boundaries
% S: structure with 16 fields, each representing a conversion type at boundary n
%              [PdPu  SdPu  PuPu  SuPu ]
%		       [PdSu  SdSu  PuSu  SuSu ]
%         S =  [PdPd  SdPd  PuPd  SuPd ] = [Rd Tu]
%			   [PdSd  SdSd  PuSd  SuSd ]   [Td Ru]
% PdPu: a matrix of classical P reflection coefficients for incidence from above.
%       indexed as PdPu(ips,ins)
% :
% SuSd: a vector of S-wave refelction coefficients for incidence from below.
% Transmission of down-going P and S above interface, Ad = [AdP;AdS]
% into down-going P and S below interface, Bd = [BdP;BdS] = Td*Ad, where
% Td = [PdPd SdPd ; PdSd SdSd];
% Likewise
% Au = Ru*Ad; Ru = [PdPu SdPu;PdSu SdSu];
% Au = Tu*Bu; Tu = [PuPu SuPu;PuSu SuSu];
% Bd = Ru*Bu; Ru = [PuPd SuPd;PuSd SuSd];
% INPUTS:
% ns : boundary numbers for which to compute the conversions
% p  : ray parameter sin(incidence angle)/alpha [s/m]
% alpha : P-velocities in layers [m/s]
% beta  : S-hastighed for alle lag [m/s]
% rho :  Densities in layers [kg/m3]
% OUTPUT:
% Td = [PdPd SdPd ; PdSd SdSd] = S(3:4,1:2)
% Ru = [PdPu SdPu ; PdSu SdSu] = S(1:2,1:2)
% Tu = [PuPu SuPu ; PuSu SuSu] = S(1:2,3:4)
% Rd = [PuPd SuPd ; PuSd SuSd] = S(2:4,3:4)
%
% CALL:  Td,Ru,Tu,Rd,S=vrho2tr(p,VPs,VSs,rho)"

using StructArrays

function vrho2tr(p,alpha,beta,rho)

#initialize some variables
nn=length(alpha)-1;
A=zeros(4,4,nn);
Td=zeros(2,2,nn);
Rd=zeros(2,2,nn);
Tu=zeros(2,2,nn);
Ru=zeros(2,2,nn);
# Weird S array
S = StructArray((PdPu = zeros(1,nn), SdPu = zeros(1,nn), PuPu = zeros(1,nn), SuPu = zeros(1,nn), PdSu = zeros(1,nn), SdSu = zeros(1,nn), PuSu = zeros(1,nn), SuSu = zeros(1,nn), PdPd = zeros(1,nn), SdPd = zeros(1,nn), PuPd = zeros(1,nn), SuPd = zeros(1,nn), PdSd = zeros(1,nn), SdSd = zeros(1,nn), PuSd = zeros(1,nn), SuSd = zeros(1,nn)))

#Main loop, for each interface (consider layer above and layer bellow)
for n=1:length(alpha)-1

    N_lay = length(alpha);
    if (n>=N_lay)|(n<0)
        println("n illegal. Number of interfaces is $N_lay")
    end

  #VELOCITIES IN THE LAYERS UNDER CONSIDERATION
  alpha1  =  alpha[n];
  alpha2  =  alpha[n+1];
  beta1  =  beta[n];
  beta2  =  beta[n+1];
  rho1   =  rho[n];
  rho2   =  rho[n+1];

  i1  =  asin(p*alpha1);
  j1  =  asin(p*beta1);
  i2  =  asin(p*alpha2);
  j2  =  asin(p*beta2);

  #THE MATRIXES
  M=zeros(4,4);
  N=zeros(4,4);
  M[1,1] = -alpha1*p;
  M[1,2] = -cos(j1);
  M[1,3] = alpha2*p;
  M[1,4] = cos(j2);

  M[2,1] = cos(i1);
  M[2,2] = -beta1*p;
  M[2,3] = cos(i2);
  M[2,4] = -beta2*p;


  M[3,1] = 2*rho1*beta1^2*p*cos(i1);
  M[3,2] = rho1*beta1*(1-2*beta1^2*p^2);
  M[3,3] = 2*rho2*beta2^2*p*cos(i2);
  M[3,4] = rho2*beta2*(1-2*beta2^2*p^2);

  M[4,1] = -rho1*alpha1*(1-2*beta1^2*p^2);
  M[4,2] = 2*rho1*beta1^2*p*cos(j1);
  M[4,3] = rho2*alpha2*(1-2*beta2^2*p^2);
  M[4,4] = -2*rho2*beta2^2*p*cos(j2);


  N[1,1] = alpha1*p;
  N[1,2] = cos(j1);
  N[1,3] = -alpha2*p;
  N[1,4] = -cos(j2);

  N[2,1] = cos(i1);
  N[2,2] = -beta1*p;
  N[2,3] = cos(i2);
  N[2,4] = -beta2*p;

  N[3,1] = 2*rho1*beta1^2*p*cos(i1);
  N[3,2] = rho1*beta1*(1-2*beta1^2*p^2);
  N[3,3] = 2*rho2*beta2^2*p*cos(i2);
  N[3,4] = rho2*beta2*(1-2*beta2^2*p^2);

  N[4,1] = rho1*alpha1*(1-2*beta1^2*p^2);
  N[4,2] = -2*rho1*beta1^2*p*cos(j1);
  N[4,3] = -rho2*alpha2*(1-2*beta2^2*p^2);
  N[4,4] = 2*rho2*beta2^2*p*cos(j2);


  A[:,:,n]   =  inv(M)*N;
  Td[:,:,n]  = A[3:4,1:2,n];
  Rd[:,:,n]  = A[1:2,1:2,n];
  Tu[:,:,n]  = A[1:2,3:4,n];
  Ru[:,:,n]  = A[3:4,3:4,n];

  ip=1; #only one ray parameter
  S.PdPu[ip,n] = A[1,1,n];
  S.SdPu[ip,n] = A[1,2,n];
  S.PuPu[ip,n] = A[1,3,n];
  S.SuPu[ip,n] = A[1,4,n];
  S.PdSu[ip,n] = A[2,1,n];
  S.SdSu[ip,n] = A[2,2,n];
  S.PuSu[ip,n] = A[2,3,n];
  S.SuSu[ip,n] = A[2,4,n];
  S.PdPd[ip,n] = A[3,1,n];
  S.SdPd[ip,n] = A[3,2,n];
  S.PuPd[ip,n] = A[3,3,n];
  S.SuPd[ip,n] = A[3,4,n];
  S.PdSd[ip,n] = A[4,1,n];
  S.SdSd[ip,n] = A[4,2,n];
  S.PuSd[ip,n] = A[4,3,n];
  S.SuSd[ip,n] = A[4,4,n];

end
return Td,Ru,Tu,Rd,S
end

## Shift
"werid funtion called shift"
function shift(fs,t_shift);
N=length(fs);
exp_fac=exp.(-im*2*pi*fs*t_shift);
exp_fac[Int64(N/2)+1]=real.(exp_fac[Int64(N/2)+1]);
return exp_fac
end

##
"""
duplicate of matlab meshgrid function
"""
function meshgrid(xin,yin)
nx=length(xin)
ny=length(yin)
xout=zeros(ny,nx)
yout=zeros(ny,nx)
for jx=1:nx
    for ix=1:ny
        xout[ix,jx]=xin[jx]
        yout[ix,jx]=yin[ix]
    end
end
return (x=xout, y=yout)
end

##
" a function to get the closest to the bottom of the model that we need"
function searchsortednearest(a,x)
   idx = searchsortedfirst(a,x)
   if (idx==1); return idx; end
   if (idx>length(a)); return length(a); end
   if (a[idx]==x); return idx; end
   if (abs(a[idx]-x) < abs(a[idx-1]-x))
      return idx
   else
      return idx-1
   end
end


##
" This function returns a clean verion of the Velocity model for the RF
computation. This includes:
- Removing the top layer (if given), actually anything above 500 m
- Removing the base of the Mantle (beyond 1000 km)
- serching for interfaces ( dvp/dz > 0.1 km/s) and avereging in needed
- km/s to m/s
- add elastic air layer for the computation"

function Fix_RF_Model(Model::Array{Float64, 2})

# Get the range of the model that is actually ndeed
Model=Model[searchsortednearest(Model[:,1],0.5):searchsortednearest(Model[:,1],1000),:];

#make columns and all to meters
zs=1000*Model[2:end,1]';
VPs=1000*Model[:,2]';
VSs=1000*Model[:,3]';
rhos=1000*Model[:,4]';

# find the seismic interfaces in the model
A=abs.(diff(round.(VPs'))) .> 200
inter_index=myfind(A);
# Keep the depths of this interfaces
zsClean=zs[inter_index]

VPsClean=zeros(size(inter_index,1)+1,size(inter_index,2))
VSsClean=zeros(size(inter_index,1)+1,size(inter_index,2))
rhosClean=zeros(size(inter_index,1)+1,size(inter_index,2))

for i=1:size(inter_index,1)+1
    if i==1
        VPsClean[i]=mean(VPs[1:inter_index[i]]);
        VSsClean[i]=mean(VSs[1:inter_index[i]]);
        rhosClean[i]=mean(rhos[1:inter_index[i]]);
    elseif i > size(inter_index,1)
        VPsClean[i]=mean(VPs[inter_index[i-1]+1:end]);
        VSsClean[i]=mean(VSs[inter_index[i-1]+1:end-1]);
        rhosClean[i]=mean(rhos[inter_index[i-1]+1:end]);

    else
        VPsClean[i]=mean(VPs[inter_index[i-1]+1:inter_index[i]]);
        VSsClean[i]=mean(VSs[inter_index[i-1]+1:inter_index[i]]);
        rhosClean[i]=mean(rhos[inter_index[i-1]+1:inter_index[i]]);
    end
end


# Add elastic air to the model on top
air_fac = 0.1; #air factor; 1=sea level; 0.1=stratosphere# Density is assumed to obey Birch's-relation:
zs = [0 zsClean']; #soil at sea level implicit; force a row
VPs  = [340*air_fac VPsClean'];
rhos = [1.2*air_fac rhosClean'];
VSs  = [(0.1*VPs[1]/sqrt(3)) VSsClean']; #%Slightly unphysical elastic air

      return zs,VPs,VSs,rhos
end

##
" Find none zeros in array"
function myfind(c)
    a = similar(c, Int)
    count = 1
    @inbounds for i in eachindex(c)
        a[count] = i
        count += (c[i] != zero(eltype(c)))
    end
    return resize!(a, count-1)
end
