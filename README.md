# SeitComp.jl
Julia software to model geophysical parameters with petrological constrains

Manual for SEITCOMP V α.0
Mariano S. Arnaiz-Rodríguez (mararnai@ucm.es or mararnai@ipgp.fr)
Javier Fullea
Universidad Complutense de Madrid
2021.

1. Introduction
Welcome to SEITCOM1D!

This package was desing to be able to model geophysical data from geochemical (peetrological) and geothermal (temperature) constrains. All the work was made during 2021 by Mariano Arnaiz-R and Javier Fullea as part of the 2018-T1/AMB/11493 Integrated Geophysical-Petrological Modelling of the Crust and Mantle at global Scale at Universidad Complutense de Madrid.

The idea of this code (so far) is to be able to model in 1D:

- Dispersion of Surface Waves, considering group [U] and phase [c] velocities for the fundamental mode and for some overtones.
- Receiver Functions (of P-to-S, S-to-P, SKS-to-P waves) to model the crust and mantle structure.
- Topography
- Heat flow data.

The input of the modelling is made of petrological and thermal data (no geophysical data is input, although we do give a option for this).

2. Needed Packages.
I will admit that it is invonient that this particular package is bound to some very particular versions of some packages but, then again, Julia is still in version 1 so there is a lot of room for improvement and testing. Be sure the install the next list of versions before using the code:

- Arpacy_jll must be the newer version. The older one has problems with large eigenvalues problem.
Pkg.add(Pkg.PackageSpec(;name="Arpack_jll", version="3.8"));
- The above library (as of today) requires the old version of Arpack as it is more stable
Pkg.add(Pkg.PackageSpec(;name="Arpack", version="0.5.0"));
- Interpolations is a great package but in version 0.13.3 they deleted the option to use repeater nodes, so... we use the 0.13.1
Pkg.add(Pkg.PackageSpec(;name="Interpolations", version="0.13.1"));

2. Inside the package:

Whn you open the directory you will find:
	data/: Directory containing all the input files
	figs/: Directory containing plots and such.
	results/: Directory containing written outputs
	src/: Directory containing all the codes (DO NOT MOVE THINGS AROUND!)
	Input_Template.jl: Input file. Modify as needed.
	Set_Up.jl: Main calling function.
	Manual.txt: This File.


3. Standard Running Sequence:

Before making new files and changing things around we will discuss how to run the package.
I have tried to make this as simple as possible for both Julians and none-Julians.

a. Make sure that in Atom you are have selected the SEITCOM1D folder. If you are new to julia know that:
	To change directories use cd("path")
	To know your current path use pwd()
	To list the contect of a directory use readdir()
	To clean the Julia terminal use clearconsole()

b. Run the line: include("Set_Up.jl"). This will load Set_Up.jl to your memory.
c. Run the command: Set_Up(). You might need to run this twice (Depending on your OS). This command will load the entire enviroment and try to force the precompilation of the functions to be used. In some OS precompilation does not work, so it the code is slow at first it will faster the second time running it.
d. Run the command Run_All(). This command runs the full sequence of SEITCOM1D.

3.1 Running by sections
-  Make sure to run "a", "b", and "c" first and then:
- Run the command Load_Input("FILE.jl"). This will load to memory all the parameters in Input_Template. which holds all the input parameters and files paths.
- Run Build_Model(), this will create a model from you inputs. Alternatively you can run Load_GM("data/AK135.txt") to load a geophysical model of your own (in this case Herrmann's AK135).
- Run the command Run(). This will perform all the forward computations.
- To view all results run Draw_Plots()
- To have all the results run Write_Outs()
- If you need to see a ternary diagram for a T,P,Rock use:
General Use:	Make_Triangle(T,P,TypeofRock,SecondInput)
For Mafic Igneous Rocks:	Make_Triangle(800,0.15,"I","Mafic")
For Felsic Igneous Rocks:	Make_Triangle(800,0.15,"I","Felsic")
For Ultramafic Igneous Rocks:	Make_Triangle(800,0.15,"I","Ultramafic")
For Sedimentary Rocks (second input is porosity [%]):	Make_Triangle(800,0.1,"S",5.0)
For Mantle Rocks:Make_Triangle(1850,24,"M")

10. FAQ

- Can I change SEITCOM1D?
	Yes! You can change anything you'd like. Be mindfull of your changes because the package might not run afterwards.
- Can you change SEITCOM1D for me?
	We might need a grant for that.
- Do you have a 2D/3D version of the code?
	Not yet. Maybe one day.
- How do I change the plots created?
	Most of the plots have been created in the GMT Julia wrapper. If you know how to used GMT you can just go to /src/main/Draw_Plots.jl and change most things there. A couple of them are created with Plots. If you need more info on how to use Plots.jl go here: http://docs.juliaplots.org/latest/
- There is a problem with two used libraries Plots and GMT.
	If you are having problems specify to Julia which package you are trying to run (e.g. Plots.plot(x,y,show=true) or GMT.plot(x,y,show=true)).
- How do I install GMT to be able to used SEITCOM1D?
	This is always changing. Go here for more information: https://github.com/GenericMappingTools/GMT.jl
- In which variable are the Surface Wave Dispersions Stored?
	Write T for the periods in seconds and SWV for the velocities. Each matrix is a mode (from fundamental SWV[:,:,1] to n higher overtone SWV[:,:,n+1]). Each column is cR UR cL UL.
- Where are my RF stored?
	There are 3 matrixes P_matrix, S_matrix, SKS_matrix. Each holds all your RF. If in doubt go to "results/Cal_*RF.txt"
- The RF look very complex! Is there any way to make them simpler?
	All (I really mean all body wave conversions) are computed by the code, so RF (particularly the P-to-S one) can have a lot of spikes. You can used the variable     Clean_Model = 0 (Clean Model 1= yes! 0 = no!) to automalically simplify your velocity structure and reduce this effect, but, a word of warning, it might change the arrival times a bit too.
