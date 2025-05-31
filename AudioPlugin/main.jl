using Pkg
Pkg.activate(".")

include("src/AudioPlugin.jl")
using .AudioPlugin


AudioPlugin.greet()