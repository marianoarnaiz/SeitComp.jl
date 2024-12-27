"Input_Template.jl: This file is a Template for the SeiTCom1D input. It is basically a file that holds
all the inputs that the user is required to provided. It also holds a lot of tips
to run the code. Please, tread lightly.
- SECTION 1: TEMPERATURE
- SECTION 2: COMPOSITION
- SECTION 3: DEPTHS
- SECTION 4: RADIAL ANISOTROPY
- SECTION 5: SURFACE WAVE DISPERSION
- SECTION 6: RECEIVER FUNCTIONS"

## SECTION 1: TEMPERATURE
## 1.a Tempertaure nodes
Temperature_Nodes=[25 1345 3630].+273.15; # Temperature nodes +273.15. Includes surfaces and bottom in C and converted to K.
Temp_Depth_Nodes= [0  200  2889].*1000; # Depths of the temperature Nodes  in km (but *1000 so in m)
LAB_index=2; # This index is which is the LAB in the previous vectors

## 1.b Temperature anomalies
# Depth range for the temperature anomaly in km but *1000 so in m
TA_z=[30 50;
      300  500;
      600 900].*1000;
# Temperature of the anomaly in C or K... its the same thing!
TA  = [-000;
       -000;
       +000];

## SECTION 2: COMPOSITION
## 2.1 Sediments composition
#                      Qz(%)  Carbo(%)  Clays(%) Porosity(%)
Sediments_Composition=[30      10         60       5];

## 2.b Crust composition  a
#Felsic rocks: Qz(%) Kfeldespar(%) Plagioclase(%)
#Mafic rocks: An(%) Cpx(%) Opc(%)
#Ultramafic rocks: Ol(%) Cpx(%) Opc(%)
Crust_Composition=[ 10  90  0  "Felsic"; #Upper Crust 1
                    20  60  20  "Mafic"]; # Lower Crust

## 2.c Mantle composition
# in %w             Al2O3 FeO
Mantle_Composition=[3.6  8.0; # Lithospheric Mantle #Fullea
                    3.6  8.0; # Asthernospheric Mantle #Fullea
                    3.6  8.0; # Transition zone Mantle % 4.5 2.0
                    3.6  8.0]; # Lower zone Mantle (Wood and Rubie,1996) & Boujibar (2016)
#P_SLAB=zTPC[:,3];
                    # Mantle_Composition=[3.75  8.0; # Lithospheric Mantle #Fullea
                    #                     4.55  8.1; # Asthernospheric Mantle #Fullea
                    #                     3.0   8.1; # Transition zone Mantle % 4.5 2.0
                    #                     3.0   8.0]; # Lower zone Mantle (Wood and Rubie,1996) & Boujibar (2016)

## 2.d Mantle Composition anomalies
# Depth range for the composition anomaly in km but *1000 so in m
MantleCA_z=[250  400;
            420  580;
            2000 2200].*1000;
# Composition anomaly in %w:
#   Al2O3 FeO
MantleCA=[-0.0 +0.0;
          -0.0 -0.0;
          +0.0 -0.0];
## SECTION 3: DEPTHS
## Depths to every interfaces. Note that this represent the depth to the base of every layer in the input!
# Crustal interphases (sed and igneous), LAB and CMB are required.
#       {Topography}  {1 Sed. layers } {4 Crustal Layers}  {LAB} {CMB}
#Depths=[     0            2             10 20 30  40        Temp_Depth_Nodes[2]    2889].*1000;
Depths=[Temp_Depth_Nodes[1]/1000  2   20   40   Temp_Depth_Nodes[2]/1000  Temp_Depth_Nodes[3]/1000].*1000;

## SECTION 4: RADIAL ANISOTROPY & Qs
## 4.1 Radial anisotropy (%) of each layer.
# Radial anisotripy is required for every layer in the crust and the mantle.
# Note that the layers in the mantle are fixed to 4 layers
#       {1 Sed. layers } {2 Crustal Layers} {4 Mantle Layers}
R_Ani=[     0.0             0.0 0.0            0.0 0.0 0.0 0.0 ]; # Radial Anisotropy at each layer

## 4.2 Grain Size Model for Qs
# Choose between the:
# - Grain size according to Schierjott 2020 ("S")
# - Grain size according to Dannberg 2017 ("D")
# - Grain size is constant 10 mm in all the mantle ("C") -> Sometimes simple is better!
Grain_Size_Model="C"

## 4.3 Melt Switch

Melt_Crust="ON" # Turn "ON" or "OFF"
Melt_Mantle="ON" # Turn "ON" or "OFF"

## SECTION 5: SURFACE WAVE DISPERSION
## 5.a Observed Surface Wave Dispersion File
Obs_SWD=readdlm("data/Obs_SWD.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB

## 5.b Inputs for Surface Wave Dispersion Forward Modelling
Anisotropy=1::Int64; #Anisotropy: Isotropic (0) or Anisotropic (1)
h_index=0; # Index that determines the model to be used in the SWF Disperison. Gold_Mesh (0), Compute (1)


## SECTION 6: RECEIVER FUNCTIONS"
## 6.a Observed Receiver Function Files
Obs_PRF=readdlm("data/No_PRF.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB
Obs_SRF=readdlm("data/No_SRF.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB
Obs_SKSRF=readdlm("data/No_SKSRF.txt",Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB

## 6.b General Inputs for Receiver functions
Gaussian_factor = 2.5; #see Ammon 1991
Clean_Model = 0 # Clean Model 1= yes! (Experimental feature, not recommended) 0 = no!
## 6.c Inputs for P wave Receiver functions
rotate_P= 1 # 1 = Rotate according to Vs value, 0 no rotaton

## end of this file
