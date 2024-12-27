"This is a function Helps run all the forward calculations"

function Run()
        try
                printstyled(" Running ... ⏳ ",color=:light_cyan)
                include("src/main/Forward.jl");
                println(" ")
                printstyled(" Forward Calculations completed!",color=:green)
        catch e2
                printstyled(" An error ocurred. Please check Input_Template.jl file",color=:red)
        end
return
end
