#By Mariano Arnaiz
"Melt Fraction:Compute the Melt Fraction.
This function compute melt fraction (xmelt) and latent heat (hlat)
at given pressure (P), temperature (T)and rock type (C).
This package was developed by Mariano Arnaiz (marianoarnaiz@gmail.com)
at the Universidad Complutense de Madrid
in spring 2022 as part of the WINTER_C project with Dr. Javier Fullea."

function Melt_fraction(zTPC,Melt_Crust,Melt_Mantle)

## Initiate variables
Melt=zeros(size(zTPC,1));
Latent_Heat=zeros(size(zTPC,1)); #latent heat
T_Solidus=zeros(size(zTPC,1))
T_Liquidus=zeros(size(zTPC,1))
for i=1:size(zTPC,1)
	# Calculate melt fraction using marker type
	P=zTPC[i,3]*1e-1; # Pressure from Bar to MPa
	global tl=0; # Liquidus temperature
	global ts=0; # Solidus temperature


    # Sediments: latent heat 300 kJ/kg
	if zTPC[i,8] == "S" && Melt_Crust == "ON"
	 	# Solidus Temperature
    	if (P<1200)
        global 	ts=889+17900/(P+54)+20200/(P+54)^2;
    	else
        global 	ts=831+0.06*P;
    	end
    	# Liquidus temperature
    	global tl=1262+0.09*P;
    	# Latent heat
    	HL=300000;
	end

	# Felsic Rocks: latent heat 300 kJ/kg
    if zTPC[i,8] == "C" && zTPC[i,7] == "Felsic"  && Melt_Crust == "ON"
    	# Solidus Temperature
    	if (P<1200)
        global 	ts=889+17900/(P+54)+20200/(P+54)^2;
    	else
        global 	ts=831+0.06*P;
    	end
    	# Liquidus temperature
    	global tl=1262+0.09*P;
    	# Latent heat
    	HL=300000;
	end

    # Mafic rock: e.g. Gabbro and Basalt: latent heat 380 kJ/kg
    if zTPC[i,8] == "C" && zTPC[i,7] == "Mafic"  && Melt_Crust == "ON"
    	# Solidus Temperature
    	if (P<1600)
        global 	ts=973-70400/(P+354)+77800000/(P+354)^2;
    	else
        global 	ts=935+0.0035*P+0.0000062*P^2;
    	end
    	# Liquidus temperature
    	global tl=1423+0.105*P;
    	# Latent heat
    	HL=380000;
	end

	# Ultramafic rock: e.g.  latent heat 380 kJ/kg (Same as mafic for now)
    if zTPC[i,8] == "C" && zTPC[i,7] == "Ultramafic"  && Melt_Crust == "ON"
    	# Solidus Temperature
    	if (P<1600)
        	global ts=973-70400/(P+354)+77800000/(P+354)^2;
    	else
        	global ts=935+0.0035*P+0.0000062*P^2;
    	end
    	# Liquidus temperature
    	global tl=1423+0.105*P;
    	# Latent heat
    	HL=380000;
	end

    # 4: THE ENTIRE MANTLE!: latent heat 400 kJ/kg
    if zTPC[i,8] == "M"  && Melt_Mantle == "ON"
		P=P*1e+1*1e-4; #Turn to GPa
		if P<10
		     # Solidus Temperature Katz
		     ts=(1085.7+132.9*P-5.1*P^2)+273.15;
		     # Liquidus temperature Latisov
		     tl=1968 + 67.1*P + -4.35*P^2 + 0.17*P^3 + -2.25e-03*P^4
		     # Correct for higher presures P > 10 as Gerya does, to avoid comeback of the function
		elseif  P >= 10.0 && P <=25
		     # Solidus Temperature HertzbergGerya
			 ts=(1086-5.7*P+390*log(P))+273.15
		     #ts=2178+0.030819*((P*1000)-10000);
		     # Liquidus temperature Liquidus temperature Latisov
		     #T_Liquidus_Matrix[i,1]=(1780+45*P-2*P^2)+273.15
			 tl=1968 + 67.1*P + -4.35*P^2 + 0.17*P^3 + -2.25e-03*P^4
		# In the Lower Mantle: Fu et al 2018
		elseif P >= 25.0
		  # Solidus Temperature Fu
		  ts=2328+14.6*P-0.0412*P^2;
		  # Liquidus temperature Fu
		  tl=2876+15.7*P-0.0421*P^2;
		end
		# Latent heat
		HL=400000;
	end

# Melt fraction and latent heat calc, check
T_Solidus[i]=ts;
T_Liquidus[i]=tl;
xmelt=0;
hlat=0;
mtk=zTPC[i,2]
#println(" ts: $ts tl: $tl")
if tl>0
	# Solidus and liquidus must not entersect
    # in the extrapolation region
    if (ts>tl-5)
        ts=tl-100
    end
    # Melt fraction
    xmelt=100*(mtk-ts)/(tl-ts);
    if (xmelt<0)
        xmelt=0;
    end
	if (xmelt>0) && (xmelt<0.1)
        xmelt=0.1;
    end
    if (xmelt>2.0) # The max melt is 2.0 %
        xmelt=2.0;
    end
	# Latent heat calc
	hlat=HL*xmelt/100;
end

Melt[i]=xmelt;
Latent_Heat[i]=hlat;

end
return  Melt, Latent_Heat, T_Solidus, T_Liquidus
end

## Correct For Melt
# In this function we correct the global Geophysical model for the Melt EfFects
# if any. Correction After Chantel et al., 2016 Science
function Melt_Model(zTPC,Model,Melt)
	# For each layer
	for i = 1:size(Model,1)
		## Correct Vp and Vs for all nodes
			Vp_melt=Model[i,2]+0.07*(Melt[i]^2)-0.5566*Melt[i]; # Vp
            Vs_melt=Model[i,3]+0.065*(Melt[i]^2)-0.5565*Melt[i]; # Vs
		## Correct Density for Melt. THis is after Gerya 2007
		# The density correction depends on the material been molten!!!!
		 if zTPC[i,8] == "S" # In case of sediments
			 Rho0solid=2.7; Rho0molten=2.4;
		 elseif zTPC[i,8] == "C" && zTPC[i,7] == "Felsic"  # For felsic rocks
			 Rho0solid=2.7; Rho0molten=2.4;
		 elseif zTPC[i,8] == "C" && zTPC[i,7] == "Mafic"  # For mafic rocks
			 Rho0solid=2.9; Rho0molten=2.6;
		 elseif zTPC[i,8] == "C" && zTPC[i,7] == "Ultramafic"  # For mafic rocks
			Rho0solid=3.25; Rho0molten=2.7;
		 elseif zTPC[i,8] == "M" # for  mantle
			Rho0solid=3.3; Rho0molten=2.7;
		 end
		Rho_melt=Model[i,4]*(1-Melt[i]+Melt[i]*(Rho0molten/Rho0solid))
		# Estimate Qs from the metl fraction

            if Melt[i] > 0.0
				Qs_melt=100/(2.4063*log(Melt[i])+5.9284) # Qs
			else
				Qs_melt=1000
			end
		## Correct the model accordingly
		Model[i,2]=Vp_melt; #Vp
		Model[i,3]=Vs_melt; #Vs
		Model[i,4]=Rho_melt; #Rho
		Model[i,5]=minimum([Model[i,5],Qs_melt])
	end
	return Model
end
