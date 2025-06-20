module AudioPlugin

using FFTW, Plots, DSP, PlutoUI, WAV, FileIO, WebIO, OffsetArrays, Statistics

export visualize_audio, visualize_frequency_spectrum, compare_audio_files, test_plot_stft, custom_distortion


function visualize_audio(file_path::String)
    audio_data = wavread(file_path)
    audio_data[1]

    signal = audio_data[1][:,1] # Tylko pierwszy kanał
    plot(signal)

end

function get_frequency_spectrum(file_path::String)
    audio_data, fs = wavread(file_path) # Wczytaj dane audio i częstotliwość próbkowania
    signal = ndims(audio_data) == 1 ? audio_data : audio_data[:, 1]

    N = length(signal)
    signal_fft = fft(signal) # Oblicz FFT sygnału
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
    signal_fft = fft(signal)   # wykres FFT pliku audio
    plot(abs.(signal_fft[20:2000]),)
end

function new_audio_path(old_path::String, suffix::String)
    # Dopisujemy końcówkę do nazwy pliku
    if(old_path[length(old_path)-3:end] != ".wav")
        error("Plik musi być w formacie WAV")
    end

    output_path = chop(old_path, head=0, tail=4)
    output_path = output_path * suffix * ".wav"

    return output_path
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

function plot_stft_frame(stft_matrix::Matrix{ComplexF64}, frame_index::Int, fs::Real)
    frame = stft_matrix[:, frame_index]
    magnitudes = abs.(frame)  # amplitudy (moduł liczb zespolonych)
    
    # Oś częstotliwości (dla FFT symetrycznej)
    n = length(frame)
    freqs = fs * (0:(n-1)) ./ n

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
    # Wykres ramki STFT (częstotliwości w danym oknie czasowym)
    plot_stft_frame(stft_data, frame_index, fs)
    println("Wyświetlono ramkę $frame_index z $(n_frames) całkowitych ramek.")

end




function equalizer_function(frequency::Float64, bass_multiplier::Float64, low_mid_multiplier::Float64, high_mid_multiplier::Float64, treble_multiplier::Float64)
    lowest_freq = 20.0
    center_bassfreq = 80.0
    center_lowmidfreq = 350.0
    center_highmidfreq = 1000.0
    center_trebelfreq = 3500.0
    highest_freq = 20000.0

    # Sprawdzenie, czy częstotliwość jest poza zakresem
    if frequency <= lowest_freq
        return 0.0
    elseif frequency >= highest_freq
        return 0.0
    end

    # Lista częstotliwości i odpowiadających im mnożników
    multiplier_list = [1.0, bass_multiplier, low_mid_multiplier, high_mid_multiplier, treble_multiplier, 1.0]
    frequency_list = [lowest_freq, center_bassfreq, center_lowmidfreq, center_highmidfreq, center_trebelfreq, highest_freq]
    
    # binary search O(nlogn)
    idx = searchsortedlast(frequency_list, frequency) + 1

    # Wyliczamy jak przemnożona powinna być dana częstotliwość w zależności od tego
    # jak się plasuje na funkcji equzlizującej
    in_between_percentage = (frequency - frequency_list[idx-1]) / (frequency_list[idx] - frequency_list[idx-1])
    after_equalizer = multiplier_list[idx-1] + (in_between_percentage * (multiplier_list[idx] - multiplier_list[idx-1]))
    return after_equalizer
end

function equalize_audio_signal(signal::Vector{Float64},
                            fs::Int;
                            bass_multiplier::Float64 = 1.0,
                            low_mid_multiplier::Float64 = 2.0,
                            high_mid_multiplier::Float64 = 1.0,
                            treble_multiplier::Float64 = 1.0,
                            frame_size::Int = 2048,
                            hop_size::Int = 1024)

    # Normalizacja
    signal = copy(signal)  # nie modyfikujemy oryginalnego sygnału
    signal ./= maximum(abs, signal)

    # STFT
    stft_data, original_len = my_stft(signal, frame_size, hop_size)

    # Częstotliwości dla każdego wiersza STFT
    freq_bins = [fs * (k - 1) / frame_size for k in 1:size(stft_data, 1)]

    # Equalizacja
    for (i, f) in enumerate(freq_bins)
        gain = equalizer_function(Float64(f),
                                Float64(bass_multiplier),
                                Float64(low_mid_multiplier),
                                Float64(high_mid_multiplier),
                                Float64(treble_multiplier))
        stft_data[i, :] .*= gain
    end

    # ISTFT i końcowa normalizacja
    reconstructed = my_istft(stft_data, frame_size, hop_size)
    reconstructed = reconstructed[1:original_len]
    reconstructed ./= maximum(abs, reconstructed)

    return reconstructed
end

function add_gain_to_signal(signal::Vector{Float64}; drive::Float64=2.0)
    rms_original = sqrt(mean(signal .^ 2))

    distorted = tanh.(drive .* signal)

    rms_new = sqrt(mean(distorted .^ 2))

    normalized = distorted .* (rms_original / rms_new)

    return normalized
end

function apply_compression_to_signal(signal::Vector{Float64}, threshold::Float64, ratio::Float64, attack_coeff=0.01, release_coeff=0.001)
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
            # Atak:
            gain = (1 - attack_coeff) * gain + attack_coeff * target_gain
        else
            # Release:
            gain = (1 - release_coeff) * gain + release_coeff * target_gain
        end

        compressed_signal[i] = sample * gain
    end
    compressed_signal = clamp.(compressed_signal, -1.0, 1.0)
    return compressed_signal
end

function custom_distortion(filepath::String, tone:: Float64, level:: Float64, gain:: Float64)
    # 5.0 = neutral, 0.0 = minimum, 10.0 = maximum


    # Funkcje pomocnicze do zamiany parametrów Distortion na wartości dla przetwarzania audio
    function calculate_gain_multiplier(gain::Float64)
        min_gain = 1.1
        max_gain = 5.0

        if gain <= 0.0
            return min_gain
        elseif gain >= 10.0
            return max_gain
        else
            return min_gain + (max_gain - min_gain) * (gain / 10.0)
        end
    end

    function calculate_level_multiplier(level::Float64)
        min_level = 0.4
        neutral_level = 1.0
        max_level = 2.0

        if level <= 0.0
            return min_level
        elseif level >= 10.0
            return max_level
        elseif level <= 5.0
            return min_level + (neutral_level - min_level) * (level / 5.0)
        else
            return neutral_level + (max_level - neutral_level) * ((level - 5.0) / 5.0)
        end
    end

    function calculate_equalizer_values(tone::Float64)
        min_bass = 4.0
        min_low_mid = 2.6
        min_high_mid = 1.5
        min_treble = 0.5

        mid_bass = 1.2
        mid_low_mid = 4.0
        mid_high_mid = 4.0
        mid_treble = 1.2

        max_bass = 0.5
        max_low_mid = 1.5
        max_high_mid = 2.6
        max_treble = 4.0

        if tone <= 0.0
            return min_bass, min_low_mid, min_high_mid, min_treble
        elseif tone >= 10.0
            return max_bass, max_low_mid, max_high_mid, max_treble
        elseif tone <= 5.0
            bass = min_bass + (mid_bass - min_bass) * (tone / 5.0)
            low_mid = min_low_mid + (mid_low_mid - min_low_mid) * (tone / 5.0)
            high_mid = min_high_mid + (mid_high_mid - min_high_mid) * (tone / 5.0)
            treble = min_treble + (mid_treble - min_treble) * (tone / 5.0)
        else
            bass = mid_bass + (max_bass - mid_bass) * ((tone - 5.0) / 5.0)
            low_mid = mid_low_mid + (max_low_mid - mid_low_mid) * ((tone - 5.0) / 5.0)
            high_mid = mid_high_mid + (max_high_mid - mid_high_mid) * ((tone - 5.0) / 5.0)
            treble = mid_treble + (max_treble - mid_treble) * ((tone - 5.0) / 5.0)
        end

        return bass, low_mid, high_mid, treble
    end

    # Wczytywanie audio z WAV
    input_signal, fs = wavread(filepath)

    is_stereo = size(input_signal, 2) == 2
    if is_stereo
        println("Stereo input detected. Converting to mono...")
        input_signal = sum(input_signal, dims=2)[:, 1] ./ 2 
    else
        println("Mono input detected.")
        input_signal = input_signal[:, 1]
    end

    input_signal ./= maximum(abs, input_signal)  # zabezpieczenie przed clippingiem

    # Kalkulacja parametrów na podstawie potencjometrów
    println("Calculating gain multiplier...")
    gain_multiplier = calculate_gain_multiplier(gain)
    println("Gain multiplier: $gain_multiplier")

    println("Calculating level multiplier...")
    level_multiplier = calculate_level_multiplier(level)
    println("Level multiplier: $level_multiplier")

    println("Calculating equalizer values...")
    bass_multiplier, low_mid_multiplier, high_mid_multiplier, treble_multiplier = calculate_equalizer_values(tone)
    println("Equalizer values: Bass: $bass_multiplier, Low Mid: $low_mid_multiplier, High Mid: $high_mid_multiplier, Treble: $treble_multiplier")

    # Przetwarzanie sygnału
    println("Applying compression...")
    compressed = apply_compression_to_signal(input_signal, 0.5, 4.0)

    println("Applying equalization...")
    equalized = equalize_audio_signal(compressed, Int(fs);
        bass_multiplier=bass_multiplier,
        low_mid_multiplier=low_mid_multiplier,
        high_mid_multiplier=high_mid_multiplier,
        treble_multiplier=treble_multiplier)

    println("Applying distortion...")
    distorted = add_gain_to_signal(equalized; drive=gain_multiplier)

    println("Managing final level...")
    output_signal = distorted .* level_multiplier
    output_signal ./= maximum(abs, output_signal)  # końcowa normalizacja

    # Zapis pliku
    out_path = new_audio_path(filepath, "_custom-dist")
    println("Saving output to $out_path")
    wavwrite(output_signal, fs, out_path)
end




end # module AudioPlugin
