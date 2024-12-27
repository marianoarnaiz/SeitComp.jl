"MT forward modelling?"

using Plots
#theme(:vibrant)
include("forwardMT1D.jl")
## Field Data:

Data=[0.001122         2138          10
      0.0029242        1200          10
      0.0079433        1000          10
      0.015625         850           10
      0.035645         700           10
      0.061376         389.05        10
      0.15996          177.83        10
      0.53456          64.565        10
      4.0738           19.498        10
      14.997           42            10
      44.444           105           10
      100              323.59        10
      302              489.78        10
      1472.3           316.23        10
      4444.4           250           10];

  Obs_F = Data[:,1]; #frequency in Hz
  Obs_ρa = Data[:,2]; #aparent resistivity Ω/m
  Obs_e = Data[:,3]; #data error



## Model
lay = @layout [ a{0.4w} grid(2,1) ]

ρ = [300 2500 0.78 3000 2500]; #resistivity Ω/m
H = [200 400 40 500]; #layer thickness m
ρa = zeros(size(Obs_ρa)); #aparent resistivity
𝚽 = zeros(size(Obs_ρa)); #phase
for i = 1 : size(Obs_F,1)
    F = Obs_F[i]; # selected frequency
    ρa_i,𝚽_i = forwardMT(ρ, H, F);
    ρa[i] = ρa_i;
    𝚽[i] = 𝚽_i;
end
𝚽=rad2deg.(𝚽);

## Function
#Model
p0=plot(Moo[:,2],Moo[:,1]/1000,linetype=:steppre,yflip=true,xaxis=:log, lc=:darkred, label="ρ",title="Resistivity Model",xlabel = "ρ (Ω/m)", ylabel="Depth (km)")

#resistivity
p1=plot(Obs_F, Obs_ρa, xaxis=:log, seriestype = :scatter,yaxis=:log, marker=:circle,title="Apparent Resistivity",label="Obs",xlabel="Frequency (Hz)", ylabel="ρa (Ω/m)",legend=:bottomleft)
plot!(Obs_F, ρa, xaxis=:log, yaxis=:log,label="Cal")
# Phase
p2=plot(Obs_F, 𝚽, xaxis=:log, lc=:black,title="Phase",label="Cal",xlabel="Frequency (Hz)", ylabel="Phase (°)",ylims=(0,90))

plot(p0,p1,p2, layout = lay,size=(700,600),bottom_margin = 0.5Plots.cm, left_margin=0.5Plots.cm,right_margin=1.0Plots.cm,top_margin=0.25Plots.cm)
