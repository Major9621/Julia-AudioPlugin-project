module AudioPlugin

using FFTW, Plots, DSP, PlutoUI, WAV, FileIO, WebIO, OffsetArrays, Statistics

export visualize_audio, visualize_frequency_spectrum, modify_audio, compare_audio_files, test_stft


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
    rms_original = sqrt(mean(signal .^ 2))

    distorted = tanh.(drive .* signal)

    rms_new = sqrt(mean(distorted .^ 2))

    normalized = distorted .* (rms_original / rms_new)

    return normalized
end

function new_audio_path(old_path::String, suffix::String)

    if(old_path[length(old_path)-3:end] != ".wav")
        error("Plik musi być w formacie WAV")
    end

    output_path = chop(old_path, head=0, tail=4)
    output_path = output_path * suffix * ".wav"

    return output_path
end


function process_audio_stft(input_path::String, gain_factor=1.0)
    # Wczytaj dane audio
    signal, fs = wavread(input_path)
    signal = ndims(signal) == 1 ? signal : signal[:, 1]  # tylko pierwszy kanał
    
    # Parametry STFT
    frame_size = 1024
    hop_size = 512
    window = hanning(frame_size)

    stft_result = stft(signal, frame_size, hop_size, window=window)

    # Modyfikacja widma – przykład: wzmacniamy wszystkie amplitudy
    #stft_result .*= gain_factor  # lub zastosuj własną funkcję na każdej kolumnie

    # Odwracamy STFT
    reconstructed = istft(stft_result, window, hop_size)

    # Konwersja na Float32 i zapis
    reconstructed = Float32.(reconstructed)
    output_path = new_audio_path(input_path, "_stft")
    wavwrite(reconstructed, fs, output_path)

    println("Zapisano zmodyfikowany plik do: $output_path")
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
    output_path = new_audio_path(file_path, "_modified")

    if isfile(output_path)
        rm(output_path) # Usuń plik, jeśli istnieje
    end

    wavwrite(modified_signal, output_path, Fs=samplerate)
    println("Zapisano zmodyfikowaną ścieżkę do: $output_path")
    
end

function istft(spectrogram::Matrix{ComplexF64}, window::Vector{Float64}, hop_size::Int)
    frame_size, num_frames = size(spectrogram)
    output_length = (num_frames - 1) * hop_size + frame_size
    signal = zeros(Float64, output_length)
    window_sum = zeros(Float64, output_length)

    for i in 0:num_frames-1
        frame = real(ifft(spectrogram[:, i + 1]))
        range = (i * hop_size + 1):(i * hop_size + frame_size)
        signal[range] .+= frame .* window
        window_sum[range] .+= window.^2
    end

    # Normalizacja (żeby nie przesterować tam, gdzie okna się nakładają)
    nonzero = window_sum .> 1e-8
    signal[nonzero] ./= window_sum[nonzero]

    return signal
end






function my_stft(signal::Vector{Float64}, frame_size::Int, hop_size::Int)
    window = hann(frame_size)
    total_length = length(signal)
    n_frames = ceil(Int, (total_length - frame_size) / hop_size) + 1
    padded_length = (n_frames - 1) * hop_size + frame_size
    padded_signal = vcat(signal, zeros(padded_length - total_length))

    stft_matrix = Matrix{ComplexF64}(undef, frame_size, n_frames)
    for i in 0:n_frames-1
        start = i * hop_size + 1
        frame = padded_signal[start:start+frame_size-1] .* window
        stft_matrix[:, i+1] = fft(frame)
    end

    return stft_matrix
end



function my_istft(stft_matrix::Matrix{ComplexF64}, frame_size::Int, hop_size::Int)
    n_frames = size(stft_matrix, 2)
    signal_length = (n_frames - 1) * hop_size + frame_size
    signal = zeros(Float64, signal_length)
    window = hann(frame_size)
    normalization = zeros(Float64, signal_length)

    for i in 0:n_frames-1
        frame = real(ifft(stft_matrix[:, i+1]))
        start = i * hop_size + 1
        signal[start:start+frame_size-1] .+= frame .* window
        normalization[start:start+frame_size-1] .+= window.^2
    end

    # Zapobieganie dzieleniu przez 0
    normalization .= map(x -> x == 0 ? 1 : x, normalization)
    signal ./= normalization
    return signal
end


function test_stft(filepath::String; frame_size=1024, hop_size=512)
    y, fs = wavread(filepath)
    y = vec(y)  # konwersja do 1D (mono)
    y ./= maximum(abs, y)  # normalizacja

    stft_data = my_stft(y, frame_size, hop_size)
    reconstructed = my_istft(stft_data, frame_size, hop_size)

    # Przycinanie lub padding
    reconstructed = reconstructed[1:length(y)]
    reconstructed ./= maximum(abs, reconstructed)  # ponowna normalizacja

    # Zapis pliku
    parts = splitext(filepath)
    outname = parts[1] * "_stft" * parts[2]
    wavwrite(reconstructed, outname, Fs=fs)

    println("Zapisano wynik do: $outname")
end







end # module AudioPlugin
