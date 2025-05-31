module AudioPlugin

using FFTW, Plots, DSP, PlutoUI, WAV, FileIO, WebIO, OffsetArrays, Statistics

export visualize_audio, visualize_frequency_spectrum, modify_audio, compare_audio_files


function visualize_audio(file_path::String)
    audio_data = wavread(file_path)
    audio_data[1]

    signal = audio_data[1][:,1] # Tylko pierwszy kanał
    plot(signal)

end

function get_frequency_spectrum(file_path::String)
    audio_data, fs = wavread(file_path)
    signal = ndims(audio_data) == 1 ? audio_data : audio_data[:, 1]

    N = length(signal)
    signal_fft = fft(signal)
    magnitudes = abs.(signal_fft)[1:div(N, 2)]
    freqs = (0:div(N, 2)-1) .* fs / N 

    return freqs, magnitudes
end

function compare_audio_files(path_1::String, path_2::String)
    if !isfile(path_1) || !isfile(path_2)
        error("Oba pliki muszą istnieć")
    end

    freqs1, original = get_frequency_spectrum(path_1)
    freqs2, distorted = get_frequency_spectrum(path_2)

    p1 = plot(freqs1, original, title="Original", xlabel="Frequency (Hz)", ylabel="Magnitude", xlim=(20, 2000))
    p2 = plot(freqs2, distorted, title="Distorted", xlabel="Frequency (Hz)", ylabel="Magnitude", xlim=(20, 2000))

    plot(p1, p2, layout=(2, 1), size=(800, 600))
end



function visualize_frequency_spectrum(file_path::String)
    audio_data = wavread(file_path)
    signal = ndims(audio_data) == 1 ? audio_data : audio_data[:, 1]
    signal_fft = fft(signal)
    plot(abs.(signal_fft[20:2000]),)
end


function apply_distortion(signal::Vector{Float64}; drive::Float64=2.0)
    # Zapamiętaj oryginalny RMS
    rms_original = sqrt(mean(signal .^ 2))

    # Dodaj distortion (np. tanh to popularna funkcja nieliniowa)
    distorted = tanh.(drive .* signal)

    # Oblicz RMS po przekształceniu
    rms_new = sqrt(mean(distorted .^ 2))

    # Znormalizuj do oryginalnego poziomu RMS
    normalized = distorted .* (rms_original / rms_new)

    return normalized
end


function modify_audio(file_path::String, gain::Float64)
    if(file_path[length(file_path)-3:end] != ".wav")
        error("Plik musi być w formacie WAV")
    end


    audio_data, samplerate = wavread(file_path)
    signal = ndims(audio_data) == 1 ? audio_data : audio_data[:, 1]

    #modified_signal = signal .* gain
    #modified_signal = clamp.(modified_signal, -1.0f0, 1.0f0)

    modified_signal = apply_distortion(signal; drive=gain)
    modified_signal = clamp.(modified_signal, -1.0f0, 1.0f0)
    modified_signal = Float32.(modified_signal)

    # Zapis
    output_path = chop(file_path, head=0, tail=4)
    output_path = output_path * "_modified.wav"

    if isfile(output_path)
        rm(output_path) # Usuń plik, jeśli istnieje
    end

    wavwrite(modified_signal, output_path, Fs=samplerate)
    println("Zapisano zmodyfikowaną ścieżkę do: $output_path")
    
end
end # module AudioPlugin
