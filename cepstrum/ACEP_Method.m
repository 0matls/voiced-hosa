function h = ACEP_Method(s, Nc)
    % Step 1: Real Cepstrum Calculation
    Cr = real_cepstrum(s);
    
    % Step 2: Long-Pass Liftering
    CLP = long_pass_lifter(Cr, Nc);
    
    % Step 3: Spectral Subtraction
    Cres = CLP;
    
    % Step 4: Peak Detection
    locs = peak_detection(Cres);
    
    % Step 5: Comb Liftering
    C_hat = comb_liftering(Cr, locs);
    
    % Step 6: Impulse Response Estimation
    h = impulse_response_estimation(C_hat);
end

function Cr = real_cepstrum(s)
    n = length(s);
    sHamming = s .* hamming(n);
    Cr = fftshift(rceps(sHamming));
end

function CLP = long_pass_lifter(Cr, Nc)
    N = length(Cr);
    len1 = floor((N - Nc)/2);  % Length of the ones on either side of zeros
    len2 = N - Nc - len1;  % Adjusting the length if N is odd
    
    lLP = [ones(1, len1), zeros(1, Nc), ones(1, len2)];
    
    if length(lLP) ~= N  % Safety check
        error('lLP is not of length N. Check the construction of lLP.');
    end
    
    CLP = Cr .* lLP;
end


function locs = peak_detection(Cres)
    thres = mean(Cres) + 3 * std(Cres);
    locs = find(Cres > thres);
end

function C_hat = comb_liftering(Cr, locs)
    lc = ones(size(Cr));
    for i = 1:length(locs)
        lc(max(1,locs(i)-2):min(length(Cr),locs(i)+2)) = 0;
    end
    C_hat = Cr .* lc;
end

function h = impulse_response_estimation(C_hat)
    H = exp(fft(C_hat));
    h = ifft(H);
end
