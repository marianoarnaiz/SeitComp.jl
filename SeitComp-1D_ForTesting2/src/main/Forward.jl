"This is the master file of C&M1D. This file will create the Geophysical Model
from the temperature and composition prescrived by the user and used it to
compute all the synthetics. Please do not modify this file unless you know what you are doing"

## Compute SWF Dispersion
#SWF
T, SWV,EG_FUNC_R_Z, EG_FUNC_R_R, EG_FUNC_L_T = SW_Disp(Model,h_index,Obs_T,Nmodes,Anisotropy);

## Compute P wave RF.
#P RF
println(" ")
printstyled("   -  🚩  Computing Receiver Functions",color=:white)
println(" ")
P_matrix=zeros(2048,size(dist_P,1)+1);
for i=1:size(dist_P,1)
time_P, P_rf, Used_model = P_RF_forward(Model,depth_P[i],dist_P[i],Gaussian_factor,Clean_Model,rotate_P);
P_matrix[:,1]=time_P;
P_matrix[:,i+1]=P_rf;
end

## Compute S wave RF.
S_matrix=zeros(2048,size(dist_S,1)+1);
try
for i=1:size(dist_S,1)
time_S, S_rf, Used_model = S_RF_forward(Model,depth_S[i],dist_S[i],Gaussian_factor,Clean_Model);
S_matrix[:,1]=time_S;
S_matrix[:,i+1]=S_rf;
end
catch e
    println(" ")
    printstyled(" S to P RF are not compatible with the Model.",color=:yellow)
    println(" ")
    printstyled(" You might need to change lines 47 and 52",color=:yellow)
    println(" ")
    printstyled(" in src/seis/S_RF_forward.jl to make it compatible.",color=:yellow)
    println(" ")
    printstyled(" Here a cutoff depth is imposed (250 km)",color=:yellow)
    println(" ")
end


## Compute SKS wave RF.
#SKS RF
SKS_matrix=zeros(4096,size(dist_SKS,1)+1);
for i=1:size(dist_SKS,1)
time_SKS, SKS_rf, Used_model = SKS_RF_forward(Model,depth_SKS[i],dist_SKS[i],Gaussian_factor,Clean_Model);
SKS_matrix[:,1]=time_SKS;
SKS_matrix[:,i+1]=SKS_rf;
end


## Time to stack stuff

## Reference values
pref=6.4;
cdeg2km=111.12;

##IASP91 velocity model
#depth    vp      vs     n
I91=[   0.00  5.800 3.360 0
   0.00  5.800 3.360 0
  10.00  5.800 3.360 4
  20.00  5.800 3.360 4
  20.00  6.500 3.750 0
  35.00  6.500 3.750 5
  35.00  8.040 4.470 0
  71.00  8.044 4.483 10
 120.00  8.050 4.500 10
 171.00  8.192 4.510 10
 210.00  8.300 4.518 10
 210.00  8.300 4.522 0
 271.00  8.523 4.628 10
 371.00  8.888 4.802 20
 410.00  9.030 4.870 10
 410.00  9.360 5.070 0
 471.00  9.565 5.199 10
 571.00  9.901 5.411 15
 660.00 10.200 5.600 13
 660.00 10.790 5.950 0
 671.00 10.819 5.979 3
 760.00 11.056 6.209 10
 760.00 11.056 6.210 0
 771.00 11.076 6.218 3
 821.00 11.164 6.256 5
 871.00 11.251 6.293 5
 921.00 11.335 6.329 5
 971.00 11.417 6.364 5
1021.00 11.498 6.397 5
1071.00 11.576 6.430 5
1121.00 11.653 6.462 5
1171.00 11.728 6.493 5
1221.00 11.801 6.524 5
1271.00 11.873 6.553 5
1321.00 11.944 6.582 5
1371.00 12.013 6.610 5
1421.00 12.080 6.637 5
1471.00 12.147 6.664 5
1521.00 12.212 6.691 5
1571.00 12.276 6.716 5
1621.00 12.340 6.742 5
1671.00 12.402 6.766 5
1721.00 12.463 6.791 5
1771.00 12.524 6.815 5];

## Fix the model as required
z=[0; diff(I91[:,1])];
vp=I91[:,2];
vs=I91[:,3];

## Stack the P-to-S RF
#  Load  the precalculated data
DATA = round.([[0 dist_P'] ; [0 depth_P'] ; P_matrix],digits=5);
Trf = DATA[3:end, 1];
c =sum(Trf.<0); ## Count the number of samples on negative times
Trf=Trf[Trf.>=0];
RFs= DATA[3:end, 2:end];
RFs=RFs[c+1:end,:];

# Get the ray parameter for each RF
p=zeros(size(RFs,2));
for i=1:size(RFs,2)
    p[i]=Get_p(DATA[2,i+1],DATA[1,i+1],"P");
end

# Compute synth times for the reference model and values
tref=cumsum((sqrt.((vs.^-2) .- (pref/cdeg2km)^2) .- sqrt.((vp.^-2) .- (pref/cdeg2km)^2)) .*z);
tnum=zeros(size(tref,1),size(RFs,2))
for i=1:size(RFs,2)
    tnum[:,i]=cumsum((sqrt.((vs.^-2) .- (p[i]/cdeg2km)^2) .- sqrt.((vp.^-2) .- (p[i]/cdeg2km)^2)) .*z);
end

# Correct the time of each RF
corr_tnum=zeros(size(Trf,1),size(RFs,2))
for i=1:size(RFs,2)
#interpolate each time to the tref
    interp_tnum=LinearInterpolation(tnum[:,i], tref);
#evaluate each time in the new tref
    corr_tnum[:,i]=interp_tnum(Trf);
end

# Make sure no NaNs
corr_tnum[isnan.(corr_tnum)] .= 0;
# To stack each RF we interpolate the result to the same time
time_stack=0:mean(diff(Trf)):floor(minimum(corr_tnum[end,:]))
move_outed_RFs=zeros(size(time_stack,1),size(RFs,2))
for i=1:size(RFs,2)
#interpolate each rf in the move out time
    interp_RFsnum=LinearInterpolation(corr_tnum[:,i], RFs[:,i]);
#evaluate each rf for the time stack vector
    move_outed_RFs[:,i]=interp_RFsnum(time_stack);
end
# Add up for stack RF
Stacked_PRF=[time_stack sum(move_outed_RFs,dims=2)];

## Stack the SKS-to-P RF
DATA = round.([[0 dist_SKS'] ; [0 depth_SKS'] ; SKS_matrix],digits=5);
Trf = DATA[3:end, 1];
c =sum(Trf.<0); ## Count the number of samples on negative times
Trf=Trf[Trf.>=0];
RFs= DATA[3:end, 2:end];
RFs=RFs[c+1:end,:];

# Get the ray parameter for each RF
p=zeros(size(RFs,2));
for i=1:size(RFs,2)
    p[i]=Get_p(DATA[2,i+1],DATA[1,i+1],"SKS");
end

# Compute synth times for the reference model and values
tref=cumsum((sqrt.((vs.^-2) .- (pref/cdeg2km)^2) .- sqrt.((vp.^-2) .- (pref/cdeg2km)^2)) .*z);
tnum=zeros(size(tref,1),size(RFs,2))
for i=1:size(RFs,2)
    tnum[:,i]=cumsum((sqrt.((vs.^-2) .- (p[i]/cdeg2km)^2) .- sqrt.((vp.^-2) .- (p[i]/cdeg2km)^2)) .*z);
end

# Correct the time of each RF
corr_tnum=zeros(size(Trf,1),size(RFs,2))
for i=1:size(RFs,2)
#interpolate each time to the tref
    interp_tnum=LinearInterpolation(tnum[:,i], tref);
#evaluate each time in the new tref
    corr_tnum[:,i]=interp_tnum(Trf);
end

# Make sure no NaNs
corr_tnum[isnan.(corr_tnum)] .= 0;
# To stack each RF we interpolate the result to the same time
time_stack=0:mean(diff(Trf)):floor(minimum(corr_tnum[end,:]))
move_outed_RFs=zeros(size(time_stack,1),size(RFs,2))
for i=1:size(RFs,2)
#interpolate each rf in the move out time
    interp_RFsnum=LinearInterpolation(corr_tnum[:,i], RFs[:,i]);
#evaluate each rf for the time stack vector
    move_outed_RFs[:,i]=interp_RFsnum(time_stack);
end
# Add up for stack RF
Stacked_SKSRF=[time_stack sum(move_outed_RFs,dims=2)];

## Stack the S-to-P RF
DATA = round.([[0 dist_S'] ; [0 depth_S'] ; S_matrix],digits=5);

#IASP91 velocity model
#depth    vp      vs     n
I91=[   0.00  5.800 3.360 0
   0.00  5.800 3.360 0
  10.00  5.800 3.360 4
  20.00  5.800 3.360 4
  20.00  6.500 3.750 0
  35.00  6.500 3.750 5
  35.00  8.040 4.470 0
  71.00  8.044 4.483 10
 120.00  8.050 4.500 10
 171.00  8.192 4.510 10
 210.00  8.300 4.518 10
 210.00  8.300 4.522 0];

# Fix the model as required
z=[0; diff(I91[:,1])];
vp=I91[:,2];
vs=I91[:,3];

# Use the 3 receiver functions for this test. Only keed the positive times. The Rest we do not need
Trf = DATA[3:end, 1];
c =sum(Trf.<0); ## Count the number of samples on negative times
Trf=Trf[Trf.>=0];

RFs= DATA[3:end, 2:end];
RFs=RFs[c+1:end,:];

# Get the ray parameter for each RF
p=zeros(size(RFs,2));

for i=1:size(RFs,2)
    p[i]=Get_p(DATA[2,i+1],DATA[1,i+1],"S");
end

# Compute synth times for the reference model and values
tref=cumsum((sqrt.((vs.^-2) .- (pref/cdeg2km)^2) .- sqrt.((vp.^-2) .- (pref/cdeg2km)^2)) .*z);
tnum=zeros(size(tref,1),size(RFs,2))
for i=1:size(RFs,2)
    tnum[:,i]=cumsum((sqrt.((vs.^-2) .- (p[i]/cdeg2km)^2) .- sqrt.((vp.^-2) .- (p[i]/cdeg2km)^2)) .*z);
end

# Correct the time of each RF
corr_tnum=zeros(size(Trf,1),size(RFs,2))

for i=1:size(RFs,2)
#interpolate each time to the tref
    interp_tnum=LinearInterpolation(tnum[:,i], tref,extrapolation_bc=Line());
#evaluate each time in the new tref
    corr_tnum[:,i]=interp_tnum(Trf);
end

# Make sure no NaNs
corr_tnum[isnan.(corr_tnum)] .= 0;
# To stack each RF we interpolate the result to the same time
time_stack=0:mean(diff(Trf)):floor(minimum(corr_tnum[end,:]))

move_outed_RFs=zeros(size(time_stack,1),size(RFs,2))

for i=1:size(RFs,2)
#interpolate each rf in the move out time
    interp_RFsnum=LinearInterpolation(corr_tnum[:,i], RFs[:,i]);
#evaluate each rf for the time stack vector
    move_outed_RFs[:,i]=interp_RFsnum(time_stack);
end

# Add up for stack RF

Stacked_SRF=[time_stack sum(move_outed_RFs,dims=2)];
