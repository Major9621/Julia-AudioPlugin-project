module AudioPlugin

using FFTW, Plots, DSP, PlutoUI, WAV, FileIO, WebIO, OffsetArrays

export greet


function greet()
    println("Hello from AudioPlugin!")
end


end # module AudioPlugin
