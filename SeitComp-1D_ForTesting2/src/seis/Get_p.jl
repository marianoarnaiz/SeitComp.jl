"Compute ray parameter for P,S or SKS phases"

using Interpolations, JLD2

@load "src/seis/Pwave_p.jld"
@load "src/seis/Swave_p.jld"
@load "src/seis/SKSwave_p.jld"

function Get_p(Depth,Distance,Phase::String)

ss=length(Distance);

if ss == 1
    p=NaN;

    ## For P wave
    if Phase == "P"

        for i=1:ss
            if Distance[i] < 30 || Distance[i] > 90
                  println("Distance of event out of range (30 to 90 for P), check line $i")
            elseif Depth[i] < 0 || Depth[i] > 600
                  println("Event depth out of range (0 to 600), check line $i")
            else
                 p= Pwave_p(Depth[i],Distance[i])
            end
        end

    elseif Phase == "S"

        for i=1:ss
            if Distance[i] < 55 || Distance[i] > 85
                  println("Distance of event out of range (55 to 85 for S), check line $i")
            elseif Depth[i] < 0 || Depth[i] > 600
                  println("Event depth out of range (0 to 600), check line $i")
            else
                 p= Swave_p(Depth[i],Distance[i])
            end
        end

    elseif Phase == "SKS"
        for i=1:ss
            if Distance[i] < 85 || Distance[i] > 120
                 println("Distance of event out of range (85 to 120 for SKS), check line $i")
            elseif Depth[i] < 0 || Depth[i] > 600
                 println("Event depth out of range (0 to 600), check line $i")
            else
                 p= SKSwave_p(Depth[i],Distance[i])
            end
        end
    end


else
    p=fill(NaN, ss);
## For P wave
if Phase == "P"

    for i=1:ss
        if Distance[i] < 30 || Distance[i] > 90
              println("Distance of event out of range (30 to 90 for P), check line $i")
        elseif Depth[i] < 0 || Depth[i] > 600
              println("Event depth out of range (0 to 600), check line $i")
        else
             p[i]= Pwave_p(Depth[i],Distance[i])
        end
    end

elseif Phase == "S"

    for i=1:ss
        if Distance[i] < 55 || Distance[i] > 85
              println("Distance of event out of range (55 to 85 for S), check line $i")
        elseif Depth[i] < 0 || Depth[i] > 600
              println("Event depth out of range (0 to 600), check line $i")
        else
             p[i]= Swave_p(Depth[i],Distance[i])
        end
    end

elseif Phase == "SKS"
    for i=1:ss
        if Distance[i] < 85 || Distance[i] > 120
             println("Distance of event out of range (85 to 120 for SKS), check line $i")
        elseif Depth[i] < 0 || Depth[i] > 600
             println("Event depth out of range (0 to 600), check line $i")
        else
             p[i]= SKSwave_p(Depth[i],Distance[i])
        end
    end
end
end
return p

end
