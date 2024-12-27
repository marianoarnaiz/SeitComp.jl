"This file makes FEM meshes for the SWD computation"

function Make_1D_Mesh(Tmin,Depths)
## Inputs
#Tmin=5; # Read min T from file
#Depths=[0 5 20 40 200  2889]; # Read Depths from the user input()

## Constants
zmax=Depths[end]; # Read CMB from file
vmin=1; # min velocity of the R wave c
#zmin=0.5
n=40; #According to Haney
a=0.5; #or 0.63
nodes= 5; # Number of nodes for each frequency

## Compute some values & get the initial Depths of the mesh()
zmin=(vmin*Tmin*a)/nodes; # Zmin is the smalles depth of a node in the model for the required computation
Nmax=ceil(n*a*log(zmax/zmin)); # Max bumber of layers
@fastmath @inbounds Z=round.(zmin.*exp.((0:1:Nmax)./(n*a)),digits=2); # Deoth vector that grows exponentially

## Polish the mesh a little bit
Z[end]=zmax; # Make the end of the vector the same as the input of the user

# Make sure that the User's Input depths appear in the mesh by moving the
# closest point to the user's value
for i=2:size(Depths,2)-1
    Ind = argmin((abs.(Z.-Depths[i])))
    Z[Ind]=Depths[i];
end

## From Depths to distance between nodes
FEM=[zmin; diff(Z,dims=1)]; # Maje the FEM mesh()

## Delete really small values & replace with the zmin
nfix=ceil(Int,sum(FEM[FEM.<zmin])/zmin);
FEM_fixed=[Z[1]; ones(nfix)*zmin; FEM[FEM.>zmin]]
return FEM_fixed

end
