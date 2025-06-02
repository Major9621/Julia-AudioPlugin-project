module AudioPlugin

using FFTW, Plots, DSP, PlutoUI, WAV, FileIO, WebIO, OffsetArrays, Statistics

export visualize_audio, visualize_frequency_spectrum, compare_audio_files, test_stft, test_compression, test_plot_stft, equalizer_function, equalize_audio


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
    p2 = plot(freqs2, distorted, title="Modified", xlabel="Frequency (Hz)", ylabel="Magnitude", xlim=(20, 2000))

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


function apply_compression(signal::Vector{Float64}, threshold::Float64, ratio::Float64, attack_coeff=0.01, release_coeff=0.001)
    compressed_signal = similar(signal)
    gain = 1.0  # początkowy gain
    for i in eachindex(signal)
        sample = signal[i]
        abs_sample = abs(sample)

        # Oblicz target_gain zależnie od przekroczenia thresholdu
        if abs_sample > threshold
            target_gain = threshold + (abs_sample - threshold) / ratio
            target_gain /= abs_sample  # znormalizowany gain < 1
        else
            target_gain = 1.0  # brak kompresji
        end

        # Wygładzanie gainu
        if target_gain < gain
            # Atak: gain spada szybciej
            gain = (1 - attack_coeff) * gain + attack_coeff * target_gain
        else
            # Release: gain rośnie wolniej
            gain = (1 - release_coeff) * gain + release_coeff * target_gain
        end

        compressed_signal[i] = sample * gain
    end
    compressed_signal = clamp.(compressed_signal, -1.0, 1.0)
    return compressed_signal
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

    return stft_matrix, total_length
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

    normalization .= map(x -> x == 0 ? 1 : x, normalization)
    signal ./= normalization
    return signal
end


function test_stft(filepath::String; frame_size=2048, hop_size=1024) # orygnalnie dwa razy mniejsze, możen się bawić
    y, fs = wavread(filepath)

    # Obsługa stereo i mono
    if ndims(y) == 2 && size(y, 2) == 2
        y = mean(y, dims=2)  # uśrednianie kanałów
        y = vec(y)
    else
        y = vec(y)
    end

    y ./= maximum(abs, y)  # normalizacja

    @info "Oryginalna długość: $(length(y))"

    stft_data, original_len = my_stft(y, frame_size, hop_size)
    reconstructed = my_istft(stft_data, frame_size, hop_size)

    reconstructed = reconstructed[1:original_len]
    reconstructed ./= maximum(abs, reconstructed)

    @info "Zrekonstruowana długość: $(length(reconstructed))"

    parts = splitext(filepath)
    outname = parts[1] * "_stft" * parts[2]
    wavwrite(reconstructed, outname, Fs=fs)

    println("Zapisano wynik do: $outname")
end



function plot_stft_frame(stft_matrix::Matrix{ComplexF64}, frame_index::Int, fs::Real)
    frame = stft_matrix[:, frame_index]
    magnitudes = abs.(frame)  # amplitudy (moduł liczb zespolonych)
    
    # Oś częstotliwości (dla FFT symetrycznej)
    n = length(frame)
    freqs = fs * (0:(n-1)) ./ n  # np. 0:22050 Hz jeśli fs=44100 i n=1024

    p = plot(freqs[1:div(n,2)], magnitudes[1:div(n,2)],
        xlabel="Częstotliwość [Hz]", ylabel="Amplituda",
        title="Widmo ramki $frame_index", legend=false,
        xlims=(0, 2000))
    
    display(p)
end



function test_plot_stft(filepath::String, percent_frame::Float64; frame_size=2048, hop_size=1024)
    y, fs = wavread(filepath)

    # Obsługa stereo i mono
    if ndims(y) == 2 && size(y, 2) == 2
        y = mean(y, dims=2)  # uśrednianie kanałów
        y = vec(y)
    else
        y = vec(y)
    end

    y ./= maximum(abs, y)  # normalizacja

    stft_data, original_len = my_stft(y, frame_size, hop_size)


    n_frames = ceil(Int, (original_len - frame_size) / hop_size) + 1
    frame_index = ceil(Int, percent_frame * n_frames)

    plot_stft_frame(stft_data, frame_index, fs)
    println("Wyświetlono ramkę $frame_index z $(n_frames) całkowitych ramek.")

end



function test_compression(filepath::String; threshold=0.5, ratio=6.0)
    y, fs = wavread(filepath)

    # Obsługa stereo i mono
    if ndims(y) == 2 && size(y, 2) == 2
        y = mean(y, dims=2)  # uśrednianie kanałów
        y = vec(y)
    else
        y = vec(y)
    end

    y ./= maximum(abs, y)  # normalizacja

    compressed_signal = apply_compression(y, threshold, ratio)

    outname = new_audio_path(filepath, "_compressed")
    wavwrite(compressed_signal, outname, Fs=fs)

    println("Zapisano wynik do: $outname")
end


function equalizer_function(frequency::Float64, bass_multiplier::Float64, low_mid_multiplier::Float64, high_mid_multiplier::Float64, treble_multiplier::Float64)
    lowest_freq = 20.0
    center_bassfreq = 80.0
    center_lowmidfreq = 350.0
    center_highmidfreq = 1000.0
    center_trebelfreq = 3500.0
    highest_freq = 20000.0

    if frequency <= lowest_freq
        return 0.0
    elseif frequency >= highest_freq
        return 0.0
    end

    multiplier_list = [1.0, bass_multiplier, low_mid_multiplier, high_mid_multiplier, treble_multiplier, 1.0]
    frequency_list = [lowest_freq, center_bassfreq, center_lowmidfreq, center_highmidfreq, center_trebelfreq, highest_freq]
    
    idx = searchsortedlast(frequency_list, frequency) + 1
    #println("Index: $idx, Frequency: $frequency")
    #println(frequency_list[idx])

    in_between_percentage = (frequency - frequency_list[idx-1]) / (frequency_list[idx] - frequency_list[idx-1])
    after_equalizer = multiplier_list[idx-1] + (in_between_percentage * (multiplier_list[idx] - multiplier_list[idx-1]))
    return after_equalizer
end



function equalize_audio(filepath::String,
                        bass_multiplier::Float64=1.0,
                        low_mid_multiplier::Float64=2.0,
                        high_mid_multiplier::Float64=1.0,
                        treble_multiplier::Float64=1.0,
                        frame_size::Int=2048,
                        hop_size::Int=1024)

    y, fs = wavread(filepath)

    if ndims(y) == 2 && size(y, 2) == 2
        y = mean(y, dims=2)
        y = vec(y)
    else
        y = vec(y)
    end

    y ./= maximum(abs, y)

    stft_data, original_len = my_stft(y, frame_size, hop_size)

    # Oblicz odpowiadające częstotliwości dla każdej linii STFT
    freq_bins = [fs * (k - 1) / frame_size for k in 1:size(stft_data, 1)]

    # Equalizacja — skalowanie magnitudy każdego widma w ramce
    for (i, f) in enumerate(freq_bins)
        gain = equalizer_function(Float64(f),
                        Float64(bass_multiplier),
                        Float64(low_mid_multiplier),
                        Float64(high_mid_multiplier),
                        Float64(treble_multiplier))
        stft_data[i, :] .*= gain
    end

    # ISTFT
    reconstructed = my_istft(stft_data, frame_size, hop_size)
    reconstructed = reconstructed[1:original_len]
    reconstructed ./= maximum(abs, reconstructed)

    # Zapis do pliku
    parts = splitext(filepath)
    outname = parts[1] * "_equalized" * parts[2]
    wavwrite(reconstructed, outname, Fs=fs)

    println("Zapisano wynik do: $outname")
end





end # module AudioPlugin
