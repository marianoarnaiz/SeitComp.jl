"This is a srction Helps set up all the enviroment. JUST RUN ONCE or TWICE!"

using DelimitedFiles, GMT, Plots, Interpolations

function Set_Up()
    ## Clear the console
    clearconsole()
    ## Make sure Julia is in the right path
    println(" ")
    println(" ")
    printstyled("⨁****************************************************⨁",color=:light_cyan)
    println(" ")
    printstyled("*                  🌎 SeiTCom-1D 🌍                    *",color=:light_cyan,bold=true)
    println(" ")
    printstyled("*                       V 1.0                        *",color=:light_cyan,bold=true)
    println(" ")
    printstyled("*                Arnaiz-R and Fullea                 *",color=:light_cyan)
    println(" ")
    printstyled("*                       (2022)                       *",color=:light_cyan)
    println(" ")
    printstyled("⨁****************************************************⨁",color=:light_cyan)
    println(" ")
    printstyled(" Welcome to SeiTCom-1D Ver. 1.0.",color=:light_cyan)
    println(" ")
    printstyled(" This is a program for modelling Geophysical data",color=:light_blue)
    println(" ")
    printstyled(" from Geochemical composition and Temperature .",color=:light_blue)
    println(" ")
    printstyled(" We are setting up the modelling enviroment.",color=:light_cyan)
    println(" ")
    printstyled(" Please, wait ...",color=:red)
    println(" ")
    println(" ")
    printstyled(" Testing all the modules ... ⏳",color=:light_cyan)
    println(" ")

## Check the working Directory
    Directory=pwd()
    cd(Directory)
## Make sure Interpolations package is the right version!!!!!!
# we need to fix this is in a efficent way
#Pkg.add(Pkg.PackageSpec(;name="Interpolations", version="0.13.1"));
#add GMT#master
## Load all the necessary modules
# Temperature Preassure Gravity
    include("src/tpg/T_P_G.jl"); # Temperature, Preassure and Gravity from top to CMB.
    include("src/tpg/Melt_Fraction.jl"); # Compute melt fraction at each depth.
# Seismology
    include("src/seis/SW_Disp.jl"); # SW Dispersion Forward Modelling
    include("src/seis/P_RF_forward.jl"); # P-to-S Receiver srction (RF) Forward Modelling
    include("src/seis/S_RF_forward.jl"); # S-to-P Receiver srction (RF) Forward Modelling
    include("src/seis/SKS_RF_forward.jl"); # SKS-to-P Receiver srction (RF) Forward Modelling
    include("src/seis/Get_p.jl"); # Compute ray parameter
    include("src/seis/Qs_Burgers.jl"); # Compute Qs according to Burgers' formulation
# Composition to properties
    include("src/composition2vpvsrho/sediments/Sediments_TPC_2_VpVsRho.jl") # Compute Vp,Vs and Rho for Sediments from T,P and Mineral Composition and Porosity
    include("src/composition2vpvsrho/crust/Crust_TPC_2_VpVsRho.jl") # Compute Vp,Vs and Rho for Sediments from T,P and Mineral Composition
    include("src/composition2vpvsrho/crust/Igneous_TPC_2_VpVsRho.jl") # # Compute Vp,Vs and Rho for Sediments from T,P and Mineral Composition and Rock Type
    include("src/composition2vpvsrho/mantle/Mantle_TPC_2_VpVsRho.jl") # Compute Vp,Vs and Rho for the Mantle from T,P and wh% Composition
# Commands
    include("src/main/Run.jl") # Forward Modelling Master srction
    include("src/main/Draw_Plots.jl") # Plotting script
    #include("Input_Template.jl") # Master input file
    include("src/main/Load_Input.jl") # Load Master input file info to memory
    include("src/main/Load_GM.jl") # Load a Model file to memory
    include("src/main/Write_Outs.jl") # Write files to ascii in results/
    include("src/main/Build_Model.jl") # Only write to memotu the geophysical model
    include("src/main/Make_Triangle.jl") # Print ternary diagrams for Vp,Vs and Rho for a S or I rock
    include("src/main/FindCz.jl") # Print ternary diagrams for specific values at some depth conditions
    include("src/main/Run_All.jl"); # Run all the commands at once.
# Others
    include("src/ineptias/tools.jl") # Load Master input file info to memory

## Load Reference models
    global PREM=readdlm("data/PREM.txt",Float64,comments=true,comment_char='#'); #This is a reference model of the PREM
    global IASP91=readdlm("data/IASP91.txt",Float64,comments=true,comment_char='#'); #This is a reference model of the IASP91
    global IASP91wLAB=readdlm("data/IASP91wLAB.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB
    global AK135=readdlm("data/AK135.txt",Float64,comments=true,comment_char='#'); #This is a reference model of the AK135

## Try to run The srctions once to force precompilation
# this may cause a problem and require to rerun Set_Up()
try
    #SWF
    #SW_Disp(Model,FEM_MESH,T0,Tstep,TF,Nmodes,Earth,Elastic,Anisotropy)
    Obs_T=10:50:300::Int64; #Shortes period (s)
    Nmodes=3::Int64; #Number of modes (Currently 1 to 3)
    Earth=1::Int64; #Earth: Flat (0) or Spherical (1)
    Elastic=1::Int64; #Elasticity: Elastic (1) or Anaelastic (0)
    Anisotropy=0::Int64; #Anisotropy: Isotropic (0) or Anisotropic (1)
    h_index=0; # Index that determines the model to be used in the SWF Disperison
    SW_Disp(AK135,h_index,Obs_T,Nmodes,Anisotropy);

    ## Inputs for Receiver srctions
    Gaussian_factor = 2.5; #see Ammon 1991
    Clean_Model = 0 # Clean Model 1= yes! 0 = no!
    rotate_P= 1 # 1 = Rotate according to Vs value, 0 no rotaton
    depth=100.0;
    dist_P=90.0;
    dist_S=85.0;
    dist_SKS=120.0;

    #RF
    println(" ")
    printstyled("   -  🚩  Computing Receiver Functions",color=:white)
    println(" ")
    # P Wave RF
    P_RF_forward(AK135,depth,dist_P,Gaussian_factor,Clean_Model,rotate_P);
    # S Wave RF
    S_RF_forward(AK135,depth,dist_S,Gaussian_factor,Clean_Model);
    # SKS wave RF
    SKS_RF_forward(AK135,depth,dist_SKS,Gaussian_factor,Clean_Model);

    ## Declare Global Variables
    #global T, SWV, P_Matrix, S_Matrix, SKS_Matrix

        println(" ")
        printstyled(" Set_Up() Report:",color=:cyan)
        println(" ")
        printstyled("   -  All the functions have been included.",color=:light_cyan)
        println(" ")
        printstyled("   -  Julia's path is in working Directory.",color=:light_cyan)
        println(" ")
        printstyled("   -  Standard Models are in Memory.",color=:light_cyan)
        println(" ")
        printstyled("   -  All the functions have been successfully precompiled.",color=:light_cyan)
        println(" ")
        printstyled("   -  All global variables declared.",color=:light_cyan)
        println(" ")
        println(" ")
        printstyled(" You are ready to begin! 👍 ",color=:green,bold=true)
        println(" ")
        println(" ")
catch e
    println(" ")
    printstyled(" An error ocurred.",color=:red)
    println(" ")
    println(" ")
    printstyled("Please run Set_Up() again. 🔁",color=:red)
    println(" ")
    sleep(1)
    #return
    #Set_Up()
end
return
end
