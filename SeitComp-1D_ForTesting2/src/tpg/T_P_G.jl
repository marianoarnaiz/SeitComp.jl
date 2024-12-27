"This file computes the temperature in the Earth by solving the heat equation
and the mantle adiabatic equation. Extra functions compute the Preassure and the
Gravity for a Density vs Depth table"

function Make_Model(Temperature_Nodes,Temp_Depth_Nodes,Depths,LAB_index,Mantle_Composition)

## Some important variables
global z=Temp_Depth_Nodes[1]:1000:Temp_Depth_Nodes[end];
global basementinx=size(Sediments_Composition,1)+1;
global mohoinx=size(Sediments_Composition,1)+size(Crust_Composition,1)+1;

## First we get an Initial Thermal model with all the constants
t0heat, kt0, rhocp, hrmv, LAB_z_index, p0, km, rhov0 = Make_Initial_Thermal_Model(z,Temperature_Nodes,Temp_Depth_Nodes,Depths,LAB_index)
## Then, we solve the Heat Equation to get an initial idea of T
# here we set TA to 0, so there is no thermal anomaly
t0,qv = Solve_Heat_Equation(t0heat,z,TA_z,zeros(size(TA_z,1)),kt0,rhocp,hrmv,LAB_z_index,Temperature_Nodes)
## Get TPC matrix
TPC=Build_TPC(z,Depths,t0,p0,Sediments_Composition,Crust_Composition,Mantle_Composition);
## Get the Density for all the layers
S=TPC[findall( x -> x== "S" , TPC[:,7]),1:6];
S[:,2]=S[:,2]*1e-4;
C=TPC[findall( x -> x== "C" , TPC[:,7]),1:6];
C[:,2]=C[:,2]*1e-4;
M=TPC[findall( x -> x== "M" , TPC[:,7]),1:4];
M[:,1]=M[:,1].-273.15; #temp in K to C
## New density model
New_Density1=[ Sediments_TPC_2_VpVsRho(S)[:,3] ; Igneous_TPC_2_VpVsRho(C)[:,3] ; Mantle_TPC_2_VpVsRho(M)[:,3] ]*1000;
## Update the Thermal model
t0heat, kt1, rhocp, hrmv, LAB_z_index, p1, km, rhov1=Make_Thermal_Model(z,Temperature_Nodes,Temp_Depth_Nodes,Depths,LAB_index,New_Density1)
## Then, we solve the Heat Equation to get an initial idea of T
# here we set TA to 0, so there is no thermal anomaly
t1,qv = Solve_Heat_Equation(t0heat,z,TA_z,zeros(size(TA_z,1)),kt1,rhocp,hrmv,LAB_z_index,Temperature_Nodes)
## Get TPC matrix
TPC=Build_TPC(z,Depths,t1,p1,Sediments_Composition,Crust_Composition,Mantle_Composition);
## Get the Density for all the layers
S=TPC[findall( x -> x== "S" , TPC[:,7]),1:6];
S[:,2]=S[:,2]*1e-4;
C=TPC[findall( x -> x== "C" , TPC[:,7]),1:6];
C[:,2]=C[:,2]*1e-4;
M=TPC[findall( x -> x== "M" , TPC[:,7]),1:4];
M[:,1]=M[:,1].-273.15; #temp in K to C
## New density model
New_Density2=[ Sediments_TPC_2_VpVsRho(S)[:,3] ; Igneous_TPC_2_VpVsRho(C)[:,3] ; Mantle_TPC_2_VpVsRho(M)[:,3] ]*1000;

## Update the Thermal model
t0heat, kt2, rhocp, hrmv, LAB_z_index, p2, km, rhov2=Make_Thermal_Model(z,Temperature_Nodes,Temp_Depth_Nodes,Depths,LAB_index,New_Density2)
## Then, we solve the Heat Equation to get an initial idea of T
# here we DO consider the thermal anomaly propoused
t2,qv2 = Solve_Heat_Equation(t0heat,z,TA_z,TA,kt2,rhocp,hrmv,LAB_z_index,Temperature_Nodes)
## Get TPC matrix
TPC=Build_TPC(z,Depths,t2,p2,Sediments_Composition,Crust_Composition,Mantle_Composition);
## Get the Density for all the layers
S=TPC[findall( x -> x== "S" , TPC[:,7]),1:6];
#S[:,2]=S[:,2]*1e-4;
S[:,2]=S[:,2]*1e-4;
C=TPC[findall( x -> x== "C" , TPC[:,7]),1:6];
C[:,2]=C[:,2]*1e-4;
#C[:,2]=C[:,2]*1e-4;
M=TPC[findall( x -> x== "M" , TPC[:,7]),1:4];
M[:,1]=M[:,1].-273.15; #temp in K to C
## Arrenge Qs
zQs=zeros(size(z));
## Arrenge the Radial Anisotropy vector to the considerted depths
z_R_Ani=[[Depths[1:end-1]'/1000 z[argmin(abs.(TPC[:,2].-140000))]/1000 z[argmin(abs.(TPC[:,2].-240000))]/1000 Depths[end]/1000]' [R_Ani[1]; R_Ani']];
int_Ani=LinearInterpolation(z_R_Ani[:,1], z_R_Ani[:,2]);
vR_Ani=int_Ani.(z/1000);

## Arrenge the Outputs
#New_Density3=[ Sediments_TPC_2_VpVsRho(S)[:,3] ; Igneous_TPC_2_VpVsRho(C)[:,3] ; Mantle_TPC_2_VpVsRho(M)[:,3] ]*1000;
Model=[z/1000 [ Sediments_TPC_2_VpVsRho(S) ; Igneous_TPC_2_VpVsRho(C) ; Mantle_TPC_2_VpVsRho(M)] zQs vR_Ani];

#update pressure one last time
gend=Gravitiy_Inside_Earth(z,Model[:,4]*1000)
#pend=Pressure_Inside_Earth(z,Model[:,4]*1000,g2)
pend=Pressure_Inside_Earth(z,Model[:,4]*1000,gend)
TPC=Build_TPC(z,Depths,t2,pend,Sediments_Composition,Crust_Composition,Mantle_Composition);

zTPC=[collect(z/1000) TPC];
# ## Arrenge the Outputs
#
# zTPC=[collect(z/1000) TPC];
# Rho=New_Density3/1000;
# Vp=[Rho2Vp.(Rho[1:Int(Depths[4]/2000)]);Mantle_VpVsRho[:,1]];
# Vs=[Vp2Vs.(Vp[1:Int(Depths[4]/2000)]);Mantle_VpVsRho[:,2]];
# Model=[collect(z/1000) Vp Vs Rho]
    return zTPC,Model,qv2
end

##
function Make_Initial_Thermal_Model(z,Temperature_Nodes,Temp_Depth_Nodes,Depths,LAB_index)

## First define the initial properites to be considered for each layer

## Radiogenic heat production, Thermal Conductivity, Heat Capacity, Density
hrm_Sediments=1.0e-6; # Radiogenic heat production, W/m^3
km_Sediments = 4.0; # Thermal conductivity, W/m/K
cpd_Sediments=775; # Heat capacity, J/kg/K
rho_0_Sediments=Sediments_TPC_2_VpVsRho([Temperature_Nodes[1].*ones(size(Sediments_Composition,1)) zeros(size(Sediments_Composition,1)) Sediments_Composition])[:,3] #initial density for the sediments



hrm_Upper_Crust=2.0e-7; # Radiogenic heat production, W/m^3
km_Upper_Crust = 2.0; # Thermal conductivity, W/m/K
cpd_Upper_Crust=1172; # Heat capacity, J/kg/K
rho_0_Crust=Igneous_TPC_2_VpVsRho([Temperature_Nodes[1].*ones(size(Crust_Composition,1)) mean(rho_0_Sediments)*Depths[basementinx]*1e-9*ones(size(Crust_Composition,1)) Crust_Composition])[:,3]; #initial density for the upper crust

hrm_Lower_Crust=2.0e-7; # Radiogenic heat production, W/m^3
km_Lower_Crust = 2.0; # Thermal conductivity, W/m/K
cpd_Lower_Crust=1140; # Heat capacity, J/kg/K
#rho_0_Lower_Crust=Crust_TPC_2_VpVsRho([13000 Crust_Composition[2,1] Crust_Composition[2,2] Crust_Composition[2,3] Depths[4]],3); #initial density for the lower crust

hrm_Lithospheric_Mantle=2e-08; # Radiogenic heat production, W/m^3
#km_Lithospheric_Mantle = Hofmeister_Mantle_Thermal_Conductivity(T,P); # Thermal conductivity, W/m/K
cpd_Lithospheric_Mantle=1005; # Heat capacity, J/kg/K
rho_0_Lithospheric_Mantle=Mantle_TPC_2_VpVsRho([1500 25000 Mantle_Composition[1,:]'],3); #initial density for the lithospheric mantle

hrm_Mantle=7e-11; # Radiogenic heat production, W/m^3
#km_Mantle = Hofmeister_Mantle_Thermal_Conductivity(T,P); # Thermal conductivity, W/m/K
cpd_Mantle=1005; # Heat capacity, J/kg/K
rho_0_Mantle=Mantle_TPC_2_VpVsRho([2500 500000 round.(mean(Mantle_Composition[2:end,:],dims=1))],3); #initial density for the mantle

## Make vectors out of these numbers
Density=[rho_0_Sediments' rho_0_Crust' rho_0_Lithospheric_Mantle rho_0_Mantle]*1000;
hrm=[hrm_Sediments*ones(basementinx-1)' hrm_Upper_Crust*ones(Int(floor(size(Crust_Composition,1)/2)))' hrm_Lower_Crust*ones(Int(ceil(size(Crust_Composition,1)/2)))' hrm_Lithospheric_Mantle hrm_Mantle];
cpd=[cpd_Sediments*ones(basementinx-1)' cpd_Upper_Crust*ones(Int(floor(size(Crust_Composition,1)/2)))' cpd_Lower_Crust*ones(Int(ceil(size(Crust_Composition,1)/2)))' cpd_Lithospheric_Mantle cpd_Mantle];
km=[km_Sediments*ones(basementinx-1)' km_Upper_Crust*ones(Int(floor(size(Crust_Composition,1)/2)))' km_Lower_Crust*ones(Int(ceil(size(Crust_Composition,1)/2)))' 0.0 0.0];

## Create a full model of properties
rhov=zeros(1,size(z,1));
hrmv=zeros(1,size(z,1));
kt=zeros(1,size(z,1));
cp0=zeros(1,size(z,1));

for j=2:size(Depths,2)
    for i=1:size(z,1)
        if z[i]<=Depths[j] && z[i]>=Depths[j-1]
            rhov[i]=Density[j-1]; # density
            hrmv[i]=hrm[j-1]; # heat production conductivity
            cp0[i]=cpd[j-1]; # Heat capacity, J/kg/K
        end
end
end

## Compute gravity and Preassure for this model
g=Gravitiy_Inside_Earth(z,rhov);
p=Pressure_Inside_Earth(z,rhov',g);

## General T gradient in Earth's Mantle
Upper_mantle_gradient=0.45; #k/km
Lower_mantle_gradient=0.25; #k/km
#We change the Gradient when the preassure jumps from upper to lower mantle
Z_Grad_Change=argmin(abs.((p*0.0001).-24));


## Make a initial Temperature Distribution
#zsize=Temp_Depth_Nodes[end]; # Model size, m
zstp=1000; # Grid step. Every 2 km. This works well for the solution of the euqation
#Creating vector for nodal point positions
#z=Temp_Depth_Nodes[1]:zstp:Temp_Depth_Nodes[end];
znum=size(z,1); # Number of nodes
LAB_z_index=argmin(abs.(z.-Temp_Depth_Nodes[LAB_index])); #get the index of the LAB in the depth (z) vector
Dpp=znum-(round(Int,250/(zstp/1000))-1); # get the index of the  D'' index (last 250 km)

# Initial temperature distribution in the Lithosphere. It begins HOT, exept
# for the first element that is COLD
global t0heat=zeros(znum,1); # Initialize the temp range
# Top and bottom
global t0heat[1]=Temperature_Nodes[1]; # the top node is cold
global t0heat[2]=Temperature_Nodes[LAB_index]; # The next node has the lithpshere temperature
global t0heat[end]=Temperature_Nodes[end]; # the bottom node is very hot

#From the third node to the Upper Mantle - Lower Mantle Boundary we add the Upper_mantle_gradient
for i=3:Z_Grad_Change
    global t0heat[i]=t0heat[i-1]+Upper_mantle_gradient*(zstp/1000);
end

#From the the Upper Mantle - Lower Mantle Boundary we add the Lower_mantle_gradient
for i=Z_Grad_Change+1:znum-1
    global t0heat[i]=t0heat[i-1]+Lower_mantle_gradient*(zstp/1000);
end

## Compute Kt for the Model
kt=Kt_Model(z,t0heat,rhov,Depths,km,p)

kappa= kt./(rhov.*cp0); # Thermal diffusivity, m^2/s
rhocp=rhov.*cp0; # volumetric heat capacity, J/K/m^3

    return t0heat, kt, rhocp, hrmv, LAB_z_index, p, km, rhov
end


##
function Make_Thermal_Model(z,Temperature_Nodes,Temp_Depth_Nodes,Depths,LAB_index,rhov)

    ## First define the initial properites to be considered for each layer

    ## Radiogenic heat production, Thermal Conductivity, Heat Capacity, Density
    hrm_Sediments=1.0e-6; # Radiogenic heat production, W/m^3
    km_Sediments = 4.0; # Thermal conductivity, W/m/K
    cpd_Sediments=775; # Heat capacity, J/kg/K



    hrm_Upper_Crust=2.0e-7; # Radiogenic heat production, W/m^3
    km_Upper_Crust = 2.0; # Thermal conductivity, W/m/K
    cpd_Upper_Crust=1172; # Heat capacity, J/kg/K

    hrm_Lower_Crust=2.0e-7; # Radiogenic heat production, W/m^3
    km_Lower_Crust = 2.0; # Thermal conductivity, W/m/K
    cpd_Lower_Crust=1140; # Heat capacity, J/kg/K

    hrm_Lithospheric_Mantle=2e-08; # Radiogenic heat production, W/m^3
    #km_Lithospheric_Mantle = Hofmeister_Mantle_Thermal_Conductivity(T,P); # Thermal conductivity, W/m/K
    cpd_Lithospheric_Mantle=1005; # Heat capacity, J/kg/K

    hrm_Mantle=7e-11; # Radiogenic heat production, W/m^3
    #km_Mantle = Hofmeister_Mantle_Thermal_Conductivity(T,P); # Thermal conductivity, W/m/K
    cpd_Mantle=1005; # Heat capacity, J/kg/K

## Make vectors out of these numbers
hrm=[hrm_Sediments*ones(basementinx-1)' hrm_Upper_Crust*ones(Int(floor(size(Crust_Composition,1)/2)))' hrm_Lower_Crust*ones(Int(ceil(size(Crust_Composition,1)/2)))' hrm_Lithospheric_Mantle hrm_Mantle];
cpd=[cpd_Sediments*ones(basementinx-1)' cpd_Upper_Crust*ones(Int(floor(size(Crust_Composition,1)/2)))' cpd_Lower_Crust*ones(Int(ceil(size(Crust_Composition,1)/2)))' cpd_Lithospheric_Mantle cpd_Mantle];
km=[km_Sediments*ones(basementinx-1)' km_Upper_Crust*ones(Int(floor(size(Crust_Composition,1)/2)))' km_Lower_Crust*ones(Int(ceil(size(Crust_Composition,1)/2)))' 0.0 0.0];

## Create a full model of properties
#rhov=zeros(1,size(z,1));
hrmv=zeros(1,size(z,1));
kt=zeros(1,size(z,1));
cp0=zeros(1,size(z,1));

for j=2:size(Depths,2)
    for i=1:size(z,1)
        if z[i]<=Depths[j] && z[i]>=Depths[j-1]
            hrmv[i]=hrm[j-1]; # heat production conductivity
            cp0[i]=cpd[j-1]; # Heat capacity, J/kg/K
        end
end
end
## Compute gravity and Preassure for this model
g=Gravitiy_Inside_Earth(z,rhov);
p=Pressure_Inside_Earth(z,rhov',g);

## General T gradient in Earth's Mantle
Upper_mantle_gradient=0.45; #k/km
Lower_mantle_gradient=0.25; #k/km
#We change the Gradient when the preassure jumps from upper to lower mantle
Z_Grad_Change=argmin(abs.((p*0.0001).-24));


## Make a initial Temperature Distribution
#zsize=Temp_Depth_Nodes[end]; # Model size, m
zstp=1000; # Grid step. Every 2 km. This works well for the solution of the euqation
#Creating vector for nodal point positions
#z=Temp_Depth_Nodes[1]:zstp:Temp_Depth_Nodes[end];
znum=size(z,1); # Number of nodes
LAB_z_index=argmin(abs.(z.-Temp_Depth_Nodes[LAB_index])); #get the index of the LAB in the depth (z) vector
Dpp=znum-(round(Int,250/(zstp/1000))-1); # get the index of the  D'' index (last 250 km)

# Initial temperature distribution in the Lithosphere. It begins HOT, exept
# for the first element that is COLD
global t0heat=zeros(znum,1); # Initialize the temp range
# Top and bottom
global t0heat[1]=Temperature_Nodes[1]; # the top node is cold
global t0heat[2]=Temperature_Nodes[LAB_index]; # The next node has the lithpshere temperature
global t0heat[end]=Temperature_Nodes[end]; # the bottom node is very hot

#From the third node to the Upper Mantle - Lower Mantle Boundary we add the Upper_mantle_gradient
for i=3:Z_Grad_Change
    global t0heat[i]=t0heat[i-1]+Upper_mantle_gradient*(zstp/1000);
end

#From the the Upper Mantle - Lower Mantle Boundary we add the Lower_mantle_gradient
for i=Z_Grad_Change+1:znum-1
    global t0heat[i]=t0heat[i-1]+Lower_mantle_gradient*(zstp/1000);
end

## Compute Kt for the Model
kt=Kt_Model(z,t0heat,rhov,Depths,km,p)

kappa= kt./(rhov.*cp0); # Thermal diffusivity, m^2/s
rhocp=rhov.*cp0; # volumetric heat capacity, J/K/m^3

    return t0heat, kt, rhocp, hrmv, LAB_z_index, p, km, rhov
end

## Solve the Heat Equation
function Solve_Heat_Equation(t0heat,z,TA_z,TA,kt,rhocp,hrmv,LAB_z_index,Temperature_Nodes)

#always solve for this time vector
dt=500*(1e+6*365.25*24*3600); # dt in Ma
#dt=[500*(1e+6*365.25*24*3600)]; # dt in Ma
zstp=z[2]-z[1];
tnum=length(dt)
znum=length(z)
t1heat=zeros(length(z))
timesum=0;
qv=0.0;
for t=1:1:tnum
        # Matrix of coefficients initialization for implicit solving
        Lheat=spzeros(znum,znum);
        # Vector of right part initialization for implicit solving
        Rheat=zeros(znum,1);
        # Implicit solving of 1D temperature equation:
        # RHO*Cp*dT/dt=d(k*dT/dx)/dx
            for i=1:1:znum
                # Global index
                # Boundary nodes
                if i==1 || i==znum || i==LAB_z_index
                    # Upper boundary
                    if i==1
                        # Constant temperature: T(i,j)=ttop
                        Lheat[i,i]=1;
                        Rheat[i,1]=Temperature_Nodes[1];
                    end
                    # Lower boundary
                    if i==znum
                        # Constant temperature: T(i,j)=ttop
                        Lheat[i,i]=1;
                        Rheat[i,1]=Temperature_Nodes[end];
                    end
                    # LAB
                    if i==LAB_z_index
                        # Constant temperature: T(i,j)=ttop
                        Lheat[i,i]=1;
                        Rheat[i,1]=Temperature_Nodes[2];
                    end
                # Internal nodes
                else
                    # RHO*Cp*dT/dt=d(k*dT/dz)/dz
                    # -d(k*dT/dz)/dz
                    Lheat[i,i-1]=-(kt[i-1]+kt[i])/2/zstp^2; # coefficient for T(i,j-1)
                    Lheat[i,i+1]=-(kt[i+1]+kt[i])/2/zstp^2; # coefficient for T(i,j+1)
                    Lheat[i,i]=(kt[i-1]+kt[i])/2/zstp^2+(kt[i+1]+kt[i])/2/zstp^2; # coefficient for T(i,j+1)
                    # RHO*Cp*(dT/dt)
                    Lheat[i,i]=Lheat[i,i]+rhocp[i]/dt; # ADD coefficient for T(i,j)
                    # Right part
                    Rheat[i,1]=(rhocp[i]*t0heat[i]/dt)+ hrmv[i];
                end
            end
        # Obtaining solution
         t1heat=Lheat\Rheat;
        # Compute surface Heat Flow: DT/Dz*k*1000
        # temperature of the 2 uppermost solid layers (not considering the air temperature in the first one)
        # and the average of the associated kt. (*1000 to get mW/m^2)
         qv=((t1heat[3]-t1heat[2])/zstp)*((kt[3]+kt[2])/2)*1000;
        # Add time counter
         timesum=timesum+dt;
        # Reassign temperature profiles for the next step
         t0heat=t1heat;
end

# We need to solve the equation once more in order to account for any thermal anomalies
dt=1*(1e+6*365.25*24*3600); # dt in Ma
Lheat=spzeros(znum,znum);
Rheat=zeros(znum,1);
#Add thermal anomaly

for j=1:size(TA,1)
    for i=1:size(t0heat,1)
        if z[i]>=TA_z[j,1] && z[i]<=TA_z[j,2]
             t0heat[i]=t0heat[i]+TA[j]; # Add temperature anomaly
        end
     end
end

# Implicit solving of 1D temperature equation:
    for i=1:1:znum
        if i==1 || i==znum || i==LAB_z_index
            if i==1
                Lheat[i,i]=1;
                Rheat[i,1]=Temperature_Nodes[1];
            end
            if i==znum
                Lheat[i,i]=1;
                Rheat[i,1]=Temperature_Nodes[end];
            end
            if i==LAB_z_index
                # Lheat[i,i]=1;
                # Rheat[i,1]=Temperature_Nodes[2];
                Lheat[i,i-1]=-(kt[i-1]+kt[i])/2/zstp^2; # coefficient for T(i,j-1)
                Lheat[i,i+1]=-(kt[i+1]+kt[i])/2/zstp^2; # coefficient for T(i,j+1)
                Lheat[i,i]=(kt[i-1]+kt[i])/2/zstp^2+(kt[i+1]+kt[i])/2/zstp^2; # coefficient for T(i,j+1)
                Lheat[i,i]=Lheat[i,i]+rhocp[i]/dt; # ADD coefficient for T(i,j)
                Rheat[i,1]=(rhocp[i]*t0heat[i]/dt)+ hrmv[i];
            end
        else
            Lheat[i,i-1]=-(kt[i-1]+kt[i])/2/zstp^2; # coefficient for T(i,j-1)
            Lheat[i,i+1]=-(kt[i+1]+kt[i])/2/zstp^2; # coefficient for T(i,j+1)
            Lheat[i,i]=(kt[i-1]+kt[i])/2/zstp^2+(kt[i+1]+kt[i])/2/zstp^2; # coefficient for T(i,j+1)
            Lheat[i,i]=Lheat[i,i]+rhocp[i]/dt; # ADD coefficient for T(i,j)
            Rheat[i,1]=(rhocp[i]*t0heat[i]/dt)+ hrmv[i];
        end
    end
# Obtaining solution
 t1heat=Lheat\Rheat;
# Compute surface Heat Flow: DT/Dz*k*1000
 qv=((t1heat[3]-t1heat[2])/zstp)*((kt[3]+kt[2])/2)*1000;
 timesum=timesum+dt;

    return t1heat, qv
end

## Compute Thermal_Conductivity as a function of T and P. Based on Hofmeister 1999
function Hofmeister_Mantle_Thermal_Conductivity(T,P)
# The formula:
#k(T,P)= ko * (298/T)^a * exp((-4*ɣ +1/3) * ∫ α(T)dT) * (1 + Ko'*P/Kt) + krad(T)
## P is required in Gpa and the code works in kbar
P=P*0.0001; #P in Gpa
## Parameters according to Fullea et al 2009
ɣ = 1.25; #thermodynamic/thermal Gru ̈neisen parameter
#α T-dependent coefficient of thermal expansion
## Mineralogy Dependant parameters
if P <= 14 #Gpa. This for alpha Olivine in the lithosphere and asthenosphere
    a = 0.25; #fitting parameter
    ko=5.0; #average ambient thermal conductivity
    Kt=129;#isothermal bulk modulus
    K0prime= 5.37; #isothermal bulk modulus pressure derivative
    a0=0.2854e-4; a1=1.0080e-8; a2=-0.3842;
elseif P > 14 && P <= 17.5 ## Gpa.  This for Garnet in the Upper Tranzition zone (410 to 550)
    a = 0.45; #fitting parameter
    ko=3.3#5.03; #average ambient thermal conductivity of olivine
    Kt=170#174;#isothermal bulk modulus
    K0prime= 3.945; #isothermal bulk modulus pressure derivative
    a0=0.2843e-4; a1=0.5772e-8; a2=-0.8901;
elseif P > 17.5 && P <= 24 ## Gpa.  This for Spinel (Wadsleyite) in the Lower Tranzition zone (550 to 660)
    a = 0.45; #fitting parameter
    ko=9.48; #average ambient thermal conductivity of Spinel
    Kt=183;#isothermal bulk modulus
    K0prime= 4.3; #isothermal bulk modulus pressure derivative
    a0=0.2497e-4; a1=0.3639e-8; a2=-0.6531;
elseif P > 24 # Gpa This for Perovskites  in the Lower Mantle (z > 660)
    a = 0.33; #fitting parameter
    ko=1.5; #average ambient thermal conductivity of Perovskites
    Kt=261;#isothermal bulk modulus
    K0prime= 4; #isothermal bulk modulus pressure derivative
    a0=0.3156e-4; a1=0.9421e-8; a2=-0.3271;
end

#krad(T)=radiation contribution to k
## Terms of the equation
#Term1: ko * (298/T)^a
Term1= ko * (298/T)^a;
#Term2: (-4*ɣ +1/3)
Term2=-4*ɣ+(1/3);
#Term3: ∫ α(T)dT
#according to Fei(1995) α(T)=ao+a1*T+a2*T^-2
# Solving the integral ∫ α(T)dT with boundaries 298 and T
#∫ α(T)dT=a0*(T - 298) - 44402*a1 + (T^2*a1)/2 - a2*(1/T - 1/298)
Term3=a0*(T - 298) - 44402*a1 + (T^2*a1)/2 - a2*(1/T - 1/298);
#Term4: (1 + Ko'*P/Kt)
Term4=(1 + K0prime*P/Kt);
#Term5: krad(T)
Term5=0.01753 - 0.00010365*T + 2.2451T*T*1e-7 - 3.407T*T*T*1e-11;

## Equation turns into
#k(T,P)= Term1 * exp(Term2 * Term3) * Term4 + Term5
k= Term1 * exp(Term2 * Term3) * Term4 + Term5;
return k
end

## Compute Thermal conductivity in a model
function Kt_Model(z,t0heat,rhov,Depths,km,p)
# Compute kt
kt=zeros(1,size(z,1));
for j=2:size(Depths,2)
    for i=1:size(z,1)
        if j<=4 && z[i]<=Depths[j] && z[i]>=Depths[j-1]
            kt[i]=km[j-1]; # thermal conductivity
        elseif j > 4 && z[i]<=Depths[j] && z[i]>=Depths[j-1]
            kt[i]=Hofmeister_Mantle_Thermal_Conductivity(t0heat[i],p[i]); # thermal conductivity

        end
end
end
return kt
end

## Compute Gratity inside the Earth
function Gravitiy_Inside_Earth(Depth,Density)
# Some constants

core_mass=1.932e24#1.895e24; #1.98e24 # mass of the earth core in kg/1.932
earth_radius=6371008.8; #earth radius in m
zsize=length(Depth);
# initiate the gravity vector (g) in m/s
g=zeros(zsize);

# Loop. For each density ring in the Earth
# first compute the mass of the ring and them the gravity force
for i=zsize:-1:1
    # For the case of the top most ring
    if i==1
         core_mass=core_mass+(Density[i])*4*pi*((earth_radius-Depth[i])^3-(earth_radius-Depth[i+1])^3)/3;
         g[i]=core_mass*6.6743e-11/(earth_radius-Depth[i])^2;
    #For all the other rings
    else
     core_mass=core_mass+(Density[i])*4*pi*((earth_radius-Depth[i-1])^3-(earth_radius-Depth[i])^3)/3;
     g[i]=core_mass*6.6743e-11/(earth_radius-Depth[i])^2;
    end
end
return g
end


## Compute Pressure Inside Earth
# Apply formula
function Pressure_Inside_Earth(Depth,Density,Gravity)
P=zeros(size(z))
H=[1000 ;diff(Depth,dims=1)];
for i=1:length(z)
    P[i]=Density[i]*Gravity[i]*H[i];
end
Pressure=cumsum(P)*1e-5 # in kbar

#Pressure=cumsum((Density.*Gravity.*[2000 ;diff(Depth,dims=1)])*1e-5,dims=1); # in kbar
return Pressure
end
## This function arreages the Temperature, Preassure and composition as needed for
# Mantle_TPC_2_VpVsRho#

# basementinx
#mohoinx
function Build_TPC(z,Depths,t1,p,Sediments_Composition,Crust_Composition,Mantle_Composition)

    TPC= Array{Any}(undef, size(t1,1),7)
    TPC[:,1:5].=0.0;
        for j=2:size(Depths,2)
            for i=1:size(t1,1)
                if j<=basementinx  && z[i]<=Depths[j] && z[i]>=Depths[j-1] #if the layer is sedimentary
                        TPC[i,1]=t1[i] # T
                        TPC[i,2]=p[i] # P
                        TPC[i,3]=Sediments_Composition[j-1,1]; # Qz
                        TPC[i,4]=Sediments_Composition[j-1,2]; # CaCO3
                        TPC[i,5]=Sediments_Composition[j-1,3]; # Clay
                        TPC[i,6]=Sediments_Composition[j-1,4]; # Porosity
                        TPC[i,7]="S"; # type of layer
                    # else
                    #     TPC[i,1]=t1[i] # T
                    #     TPC[i,2]=p[i] # P
                    #     TPC[i,3]=0.5*(TPC[i,3]+TPC[i,3]+Sediments_Composition[j-1,1]); # Qz
                    #     TPC[i,4]=0.5*(TPC[i,4]+Sediments_Composition[j-1,2]); # CaCO3
                    #     TPC[i,5]=0.5*(TPC[i,5]+Sediments_Composition[j-1,3]); # Clay
                    #     TPC[i,6]=0.5*(TPC[i,6]+Sediments_Composition[j-1,4]); # Porosity
                    #     TPC[i,7]="S"; # type of layer
                    #end
                elseif j>basementinx && j<=mohoinx  && z[i]<=Depths[j] && z[i]>=Depths[j-1] #if the layer is in the upper crust
                    #println("crust")
                    TPC[i,1]=t1[i] # T
                    TPC[i,2]=p[i] # P
                    TPC[i,3]=Crust_Composition[j-basementinx,1]; # Mineral1
                    TPC[i,4]=Crust_Composition[j-basementinx,2]; # Mineral2
                    TPC[i,5]=Crust_Composition[j-basementinx,3]; # Mineral3
                    TPC[i,6]=Crust_Composition[j-basementinx,4]; # Type of Rock
                    TPC[i,7]="C"; # type of layer
                elseif j > mohoinx && z[i]>Depths[mohoinx]
                    if t1[i] <= Temperature_Nodes[2] # if the layer is within the lithospheric mantle
                        #println("litmant")
                        TPC[i,1]=t1[i] # T
                        TPC[i,2]=p[i] # P
                        TPC[i,3]=Mantle_Composition[1,1]; # Al203
                        TPC[i,4]=Mantle_Composition[1,2]; # FeO
                        TPC[i,7]="M"; # type of layer
                    elseif t1[i] > Temperature_Nodes[2] && p[i] <= 140000 # if the layer is within the asthenosphere
                            #println("astmant")
                            TPC[i,1]=t1[i] # T
                            TPC[i,2]=p[i] # P
                            TPC[i,3]=Mantle_Composition[2,1]; # Al203
                            TPC[i,4]=Mantle_Composition[2,2]; # FeO
                            TPC[i,7]="M"; # type of layer
                    elseif p[i] > 140000 && p[i] <= 240000 # if the layer is within the Transition zone
                        #println("midmant")
                            TPC[i,1]=t1[i] # T
                            TPC[i,2]=p[i] # P
                            TPC[i,3]=Mantle_Composition[3,1]; # Al203
                            TPC[i,4]=Mantle_Composition[3,2]; # FeO
                            TPC[i,7]="M"; # type of layer
                    elseif p[i] > 240000 # if the layer is in the lower mantle
                        #println("lowmant")
                            TPC[i,1]=t1[i] # T
                            TPC[i,2]=p[i] # P
                            TPC[i,3]=Mantle_Composition[4,1]; # Al203
                            TPC[i,4]=Mantle_Composition[4,2]; # FeO
                            TPC[i,7]="M"; # type of layer
                        end
                    end

                end
       end

       #Add Compositional anomaly
       for j=1:size(MantleCA,1)
           for i=1:size(t0heat,1)
               if z[i]>=MantleCA_z[j,1] && z[i]<=MantleCA_z[j,2]
                    TPC[i,3]=TPC[i,3]+MantleCA[j,1]; # Add temperature anomaly
                    TPC[i,4]=TPC[i,4]+MantleCA[j,2];
               end
            end
       end

    return TPC
end



## Vp to Vs: From Brotcher 2005
# All in km/s
function Vp2Vs(Vp)
Vs= 0.7858 - 1.2344*Vp + 0.7949*Vp^2 - 0.1238*Vp^3 + 0.0064*Vp^4
    return Vs
end

## Vp to Rho
function Vp2Rho(Vp)
    Rho = 1.6612*Vp - 0.4721*Vp^2 + 0.0671*Vp^3 - 0.0043*Vp^4 + 0.000106*Vp^5;
    return Rho
end
## Vp of Igneous rocks based in composition. The inputs are given in %wt and p in kbar
#
# function Crust_Vp(p,SiO2,MgO,CaO,z)
# # Coeficients of the equation
# Coef=[2.92e+01	2.40e+02	2.90e+02
# 1.54e-02	2.00e-02	1.67e-02
# 2.28e-02	9.83e-03	1.35e-02
# -2.66e+00	-2.87e+01	-3.43e+01
# 1.31e-01	1.39e+00	1.65e+00
# -3.28e-03	-3.45e-02	-4.05e-02
# 4.40e-05	4.62e-04	5.42e-04
# -3.02e-07	-3.21e-06	-3.76e-06
# 8.40e-10	9.05e-09	1.06e-08
# 3.66e-01	2.48e+00	2.71e+00
# -4.96e-02	-1.71e-01	-3.50e-01
# 7.43e-03	6.20e-03	1.18e-02
# -3.07e-04	-5.35e-04	-5.52e-04
# -8.73e-03	-7.86e-02	-8.33e-02
# 3.92e-05	5.83e-04	6.72e-04
# -1.34e-03	4.39e-03	1.01e-02
# 2.37e-05	-2.82e-05	-9.34e-05];
#
# # Na2O(z)
# Na2O= 5.522e-5*((z/1000)^3) -0.0061*((z/1000)^2) + 0.2007*(z/1000) + 1.8111;
#
# # Compute Vp from the composition and account for p and Na2O if necessary
# if p < 12000
#     Vp = 6.90 - 0.011*SiO2 + 0.037*MgO+ 0.045*CaO;
# elseif p >= 12000 # && p < 15000
#     Vp = Coef[1,1] +Coef[2,1]*MgO+Coef[3,1]*CaO+Coef[4,1]*SiO2+Coef[5,1]*SiO2^2 +Coef[6,1]*SiO2^3 +Coef[7,1]*SiO2^4 + Coef[8,1]*SiO2^5 + Coef[9,1]*SiO2^6 + Coef[10,1]*Na2O + Coef[11,1]*Na2O^2 + Coef[12,1]*Na2O^3 + Coef[13,1]*Na2O^4 + Coef[14,1]*SiO2*Na2O + Coef[15,1]*(SiO2^2)*Na2O + Coef[16,1]*SiO2*(Na2O^2) + Coef[17,1]*(SiO2^2)*(Na2O^2);
# # elseif p >= 15000 &&  p < 20000
# #      Vp = Coef[1,2] +Coef[2,2]*MgO+Coef[3,2]*CaO+Coef[4,2]*SiO2+Coef[5,2]*SiO2^2 +Coef[6,2]*SiO2^3 +Coef[7,2]*SiO2^4 + Coef[8,2]*SiO2^5 + Coef[9,2]*SiO2^6 + Coef[10,2]*Na2O + Coef[11,2]*Na2O^2 + Coef[12,2]*Na2O^3 + Coef[13,2]*Na2O^4 + Coef[14,2]*SiO2*Na2O + Coef[15,2]*(SiO2^2)*Na2O + Coef[16,2]*SiO2*(Na2O^2) + Coef[17,2]*(SiO2^2)*(Na2O^2);
# # elseif p >= 20000
# #      Vp = Coef[1,3] +Coef[2,3]*MgO+Coef[3,3]*CaO+Coef[4,3]*SiO2+Coef[5,3]*SiO2^2 +Coef[6,3]*SiO2^3 +Coef[7,3]*SiO2^4 + Coef[8,3]*SiO2^5 + Coef[9,3]*SiO2^6 + Coef[10,3]*Na2O + Coef[11,3]*Na2O^2 + Coef[12,3]*Na2O^3 + Coef[13,3]*Na2O^4 + Coef[14,3]*SiO2*Na2O + Coef[15,3]*(SiO2^2)*Na2O + Coef[16,3]*SiO2*(Na2O^2) + Coef[17,3]*(SiO2^2)*(Na2O^2);
# end
#     return Vp
# end

## Old way to solve heat equation
# ## Solve the Heat Equation
# function Solve_Heat_Equation(t0heat,z,TA_z,TA,kt,rhocp,hrmv,LAB_z_index,Temperature_Nodes)
#
# #always solve for this time vector
# dt=10*(1e+6*365.25*24*3600); # dt in Ma
# zstp=z[2]-z[1];
# tnum=length(dt)
# znum=length(z)
# t1heat=zeros(length(z))
# timesum=0;
# qv=0.0;
# while t0heat[LAB_z_index]-Temperature_Nodes[2]>0
#         # Matrix of coefficients initialization for implicit solving
#         Lheat=spzeros(znum,znum);
#         # Vector of right part initialization for implicit solving
#         Rheat=zeros(znum,1);
#         # Implicit solving of 1D temperature equation:
#         # RHO*Cp*dT/dt=d(k*dT/dx)/dx
#             for i=1:1:znum
#                 # Global index
#                 # Boundary nodes
#                 if i==1 || i==znum || i==LAB_z_index
#                     # Upper boundary
#                     if i==1
#                         # Constant temperature: T(i,j)=ttop
#                         Lheat[i,i]=1;
#                         Rheat[i,1]=Temperature_Nodes[1];
#                     end
#                     # Lower boundary
#                     if i==znum
#                         # Constant temperature: T(i,j)=ttop
#                         Lheat[i,i]=1;
#                         Rheat[i,1]=Temperature_Nodes[end];
#                     end
#                     # LAB
#                     if i==LAB_z_index
#                         # Constant temperature: T(i,j)=ttop
#                         Lheat[i,i]=1;
#                         Rheat[i,1]=Temperature_Nodes[2];
#                     end
#                 # Internal nodes
#                 else
#                     # RHO*Cp*dT/dt=d(k*dT/dz)/dz
#                     # -d(k*dT/dz)/dz
#                     Lheat[i,i-1]=-(kt[i-1]+kt[i])/2/zstp^2; # coefficient for T(i,j-1)
#                     Lheat[i,i+1]=-(kt[i+1]+kt[i])/2/zstp^2; # coefficient for T(i,j+1)
#                     Lheat[i,i]=(kt[i-1]+kt[i])/2/zstp^2+(kt[i+1]+kt[i])/2/zstp^2; # coefficient for T(i,j+1)
#                     # RHO*Cp*(dT/dt)
#                     Lheat[i,i]=Lheat[i,i]+rhocp[i]/dt; # ADD coefficient for T(i,j)
#                     # Right part
#                     Rheat[i,1]=(rhocp[i]*t0heat[i]/dt)+ hrmv[i];
#                 end
#             end
#         # Obtaining solution
#          t1heat=Lheat\Rheat;
#         # Compute surface Heat Flow: DT/Dz*k*1000
#         # temperature of the 2 uppermost solid layers (not considering the air temperature in the first one)
#         # and the average of the associated kt. (*1000 to get mW/m^2)
#          qv=((t1heat[3]-t1heat[2])/zstp)*((kt[3]+kt[2])/2)*1000;
#         # Add time counter
#          timesum=timesum+dt;
#         # Reassign temperature profiles for the next step
#          t0heat=t1heat;
# end
#
# # We need to solve the equation once more in order to account for any thermal anomalies
# dt=1*(1e+6*365.25*24*3600); # dt in Ma
# Lheat=spzeros(znum,znum);
# Rheat=zeros(znum,1);
# #Add thermal anomaly
#
# for j=1:size(TA,1)
#     for i=1:size(t0heat,1)
#         if z[i]>=TA_z[j,1] && z[i]<=TA_z[j,2]
#              t0heat[i]=t0heat[i]+TA[j]; # Add temperature anomaly
#         end
#      end
# end
#
# # Implicit solving of 1D temperature equation:
#     for i=1:1:znum
#         if i==1 || i==znum || i==LAB_z_index
#             if i==1
#                 Lheat[i,i]=1;
#                 Rheat[i,1]=Temperature_Nodes[1];
#             end
#             if i==znum
#                 Lheat[i,i]=1;
#                 Rheat[i,1]=Temperature_Nodes[end];
#             end
#             if i==LAB_z_index
#                 Lheat[i,i]=1;
#                 Rheat[i,1]=Temperature_Nodes[2];
#             end
#         else
#             Lheat[i,i-1]=-(kt[i-1]+kt[i])/2/zstp^2; # coefficient for T(i,j-1)
#             Lheat[i,i+1]=-(kt[i+1]+kt[i])/2/zstp^2; # coefficient for T(i,j+1)
#             Lheat[i,i]=(kt[i-1]+kt[i])/2/zstp^2+(kt[i+1]+kt[i])/2/zstp^2; # coefficient for T(i,j+1)
#             Lheat[i,i]=Lheat[i,i]+rhocp[i]/dt; # ADD coefficient for T(i,j)
#             Rheat[i,1]=(rhocp[i]*t0heat[i]/dt)+ hrmv[i];
#         end
#     end
# # Obtaining solution
#  t1heat=Lheat\Rheat;
# # Compute surface Heat Flow: DT/Dz*k*1000
#  qv=((t1heat[3]-t1heat[2])/zstp)*((kt[3]+kt[2])/2)*1000;
#  timesum=timesum+dt;
#
#     return t1heat, qv
# end
