"This is a function Helps to load a Geophysical Model to memory. Notice that Load_input()
is necesarry for this option to function!"

function Load_GM(fileandpath::String)
        try
        global Model=readdlm(fileandpath,Float64,comments=true,comment_char='#'); #This is the IASP91 with LAB
        printstyled(" Geophysical Model correctly loaded!",color=:green)
        catch e2
        printstyled(" An error ocurred. Please check your model file",color=:red)
        end
return Model
end
