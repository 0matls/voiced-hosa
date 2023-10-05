% function [voiceQuefrencies, shift] = computeQuefrenciesForVoice(k, top_k_Quefrencies)
%     Construct the voice filename
%     filename = sprintf('voice%03d.txt', k);
% 
%     Load the voice data (assuming it's a single column txt file)
%     soundVector = load(filename);
% 
%     n = length(soundVector);
%     Fs = 8000;  % Sample rate (modify if different for other files)
%     Ts = 1/Fs;
%     periods = 5;
%     startingTime = 2;
%     estimatePitch = 0.006;
%     sFrame = soundVector(round(startingTime/Ts):round((startingTime + periods*estimatePitch)/Ts));
%     nFrame = length(sFrame);
%     shift = (-nFrame/2:nFrame/2 - 1) * (Ts) * 1000;
% 
%     sHamming = sFrame .* hamming(nFrame);
%     sRcepsHamming = fftshift(rceps(sHamming)); % Real cepstrum with hamming
%     voiceQuefrencies = sRcepsHamming;
%     if nargin > 1
%         shift = shift + 0.0625;
%         [peaks, locs] = findpeaks(sRcepsHamming);  % find peaks on the positive side
%         [~, idx] = maxk(peaks, 2*top_k_Quefrencies);
%         voiceQuefrencies = shift(locs(idx)); % Return the shifted positions of the top peaks
% 
%         Sort the voiceQuefrencies in descending order based on the peak heights
%         [~, order] = sort(peaks(idx), 'descend');
%         voiceQuefrencies = voiceQuefrencies(order);
%     end
% 
% end

function [voiceQuefrencies, shift] = computeQuefrenciesForVoice(k, top_k_Quefrencies)
    % Construct the voice filename
    filename = sprintf('voice%03d.txt', k);

    % Load the voice data (assuming it's a single column txt file)
    soundVector = load(filename);

    % Extract a frame from the soundVector for analysis
    Fs = 8000;  % Sample rate (modify if different for other files)
    Ts = 1/Fs;
    periods = 5;
    startingTime = 2;
    estimatePitch = 0.006;
    sFrame = soundVector(round(startingTime/Ts):round((startingTime + periods*estimatePitch)/Ts));

    % Compute quefrencies using the ACEP_Method
    h = ACEP_Method(sFrame, 100);
    enhancedVoiceQuefrencies = h;

    nFrame = length(enhancedVoiceQuefrencies);
    shift = (-nFrame/2:nFrame/2 - 1) * (Ts) * 1000;

    if nargin > 1
        [peaks, locs] = findpeaks(enhancedVoiceQuefrencies);  % find peaks on the positive side
        [~, idx] = maxk(peaks, 2*top_k_Quefrencies);
        voiceQuefrencies = shift(locs(idx)); % Return the shifted positions of the top peaks

        % Sort the voiceQuefrencies in descending order based on the peak heights
        [~, order] = sort(peaks(idx), 'descend');
        voiceQuefrencies = voiceQuefrencies(order);
    else
        voiceQuefrencies = enhancedVoiceQuefrencies;
    end
end

function enhancedQuefrencies = compute_quefrencies_from_impulse_response(h)
    % Placeholder for computing quefrencies from impulse response
    enhancedQuefrencies = abs(h);  
end
