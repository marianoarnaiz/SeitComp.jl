"This is a function Runs all the Main instructions to aid in modelling"

function Run_All(FILE::String)
        try
                printstyled(" Running All Instructions ... ⏳ ",color=:yellow)
                println(" ")
                Load_Input(FILE);
                println(" ")
                Build_Model();
                println(" ")
                Run();
                println(" ")
                Draw_Plots();
                println(" ")
                Write_Outs();
                println(" ")
        catch e2
                println(" ")
                printstyled(" An error ocurred. Please check everything",color=:red)
                println(" ")
        end
return
end
