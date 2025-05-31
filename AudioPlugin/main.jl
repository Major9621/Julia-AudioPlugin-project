using Pkg
Pkg.activate(".")

include("src/AudioPlugin.jl")
using .AudioPlugin



audio_name = "guitar.wav"
audio_dir = "audio"
audio_path = joinpath(audio_dir, audio_name)

if !isfile(audio_path)
    error("Audio file not found: $audio_path")
end


#AudioPlugin.visualize_audio(audio_path)

#AudioPlugin.visualize_frequency_spectrum(audio_path)

AudioPlugin.modify_audio(audio_path, 200.0)

AudioPlugin.compare_audio_files(audio_path, joinpath(audio_dir, "guitar_modified.wav"))


