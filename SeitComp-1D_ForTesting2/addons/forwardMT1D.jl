
function forwardMT1D(ρ, H, ℐ)

# Compute some parameters
mu = 4*pi*1e-7; #Magnetic Permeability (H/m)
ω = 2 * pi * ℐ; #Angular Frequency (Radians);
n=size(ρ,2); #Number of Layers

#initialize impedances vector
impedances = zeros(ComplexF64,n,1);
#Layering in this format
#  Layer     j
# Layer 1    1
# Layer 2    2
# Layer 3    3
# Layer 4    4
# Halfspahce   5
#

# Steps for modelling (for each geoelectric model and frequency)
# 1. Compute basement impedance Zn using sqrt((w * mu * resistivity))
# 2. Iterate from bottom layer to top(not the basement)
    # 2.1. Calculate induction parameters
    # 2.2. Calculate Exponential factor from intrinsic impedance
    # 2.3 Calculate reflection coeficient using current layer
    #          intrinsic impedance and the below layer impedance

# 3. Compute apparent resistivity from top layer impedance
        #   apparent resistivity = (Zn^2)/(mu * w)

#Symbols
# Zn - Basement Impedance
# Zi - Layer Impedance
# wi - Intrinsic Impedance
# di - Induction parameter
# ei - Exponential Factor
# ri - Reflection coeficient
# re - Earth R.C.

#Step 1 : Calculate halfspace impedance
Zn = sqrt(im*ω*mu*ρ[n]);
impedances[n] = Zn;

#Iterate through layers starting from layer j=n-1 (i.e. the layer above the halfspace)
for j = n-1:-1:1
    ρ_j = ρ[j];
    H_j = H[j];

    # 3. Compute apparent resistivity from top layer impedance
    #Step 2. Iterate from bottom layer to top(not the basement)
    # Step 2.1 Calculate the intrinsic impedance of current layer
    dj = sqrt(im* (ω * mu * (1/ρ_j)));
    wj = dj * ρ_j;

    # Step 2.2 Calculate Exponential factor from intrinsic impedance
    ej = exp(-2*H_j*dj);

    # Step 2.3 Calculate reflection coeficient using current layer
    #          intrinsic impedance and the below layer impedance
    belowImpedance = impedances[j + 1];
    rj = (wj - belowImpedance)/(wj + belowImpedance);
    re = rj*ej;
    Zj = wj * ((1 - re)/(1 + re));
    impedances[j] = Zj;
end
# Step 3. Compute apparent resistivity from top layer impedance
Z = impedances[1];
absZ = abs(Z);
ρa_i = (absZ * absZ)/(mu * ω); # get apparent resistivity
𝚽_i = atan(imag(Z),real(Z)); # get phase

return ρa_i,𝚽_i

end
