using Pkg
Pkg.activate(".")

include("src/AudioPlugin.jl")
using .AudioPlugin



audio_name = "overdriven riff fxoff 02.06.2025.wav"
audio_name2 = "clean riff fxoff 02.06.2025.wav"
audio_dir = "audio"
audio_path = joinpath(audio_dir, audio_name)

if !isfile(audio_path)
    error("Audio file not found: $audio_path")
end


#AudioPlugin.visualize_audio(audio_path)

#AudioPlugin.visualize_frequency_spectrum(audio_path)

#AudioPlugin.compare_audio_files(audio_path, joinpath(audio_dir, "guitar_modified.wav"))

#AudioPlugin.test_plot_stft(audio_path, 0.6, frame_size=2048, hop_size=1024)


#AudioPlugin.compare_audio_files(audio_path, joinpath(audio_dir, "guitar_equalized.wav"))


#                                        Tone, Level, Gain
AudioPlugin.custom_distortion(audio_path, 10.0, 6.0, 5.0)
