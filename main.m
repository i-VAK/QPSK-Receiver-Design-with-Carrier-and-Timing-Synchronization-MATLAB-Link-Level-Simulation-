clc
clear
close all

function [symOut, mu_hist, e_hist] = el_simple(x, sps, tau_hat)

x = x(:);
N = length(x);

mu   = tau_hat-1;         % fractional clock phase
step = 0.1;       % phase correction step
thr  = 0.02;      % error threshold
win  = 10;        % averaging window (symbols)

symOut = zeros(ceil(N/sps),1);
mu_hist = [];     % <-- clock phase history
e_hist = [];      % TED history (optional)

k = 1;
acc = 0;
cnt = 0;

while (mu + sps/8) < N-1
    
    % interpolate
    x_center = interp1(0:N-1, x, mu,        'linear', 0);
    x_early  = interp1(0:N-1, x, mu-sps/8,  'linear', 0);
    x_late   = interp1(0:N-1, x, mu+sps/8,  'linear', 0);

    % Early–Late error
    e = real(x_center)*(real(x_late)-real(x_early)) + ...
        imag(x_center)*(imag(x_late)-imag(x_early));

    e_hist(end+1) = e;
    mu_hist(end+1) = mod(mu,sps);     % normalize to 0…2 for plotting

    acc = acc + e;
    cnt = cnt + 1;

    % correction after every window
    if cnt == win
        e_avg = acc / win;

        if e_avg > thr
            % sampling late -> move earlier
            mu = mu + sps - step;
        elseif e_avg < -thr
            % sampling early -> move later
            mu = mu + sps + step;
        else
            mu = mu + sps;
        end

        acc = 0;
        cnt = 0;
    else
        mu = mu + sps;
    end

    symOut(k) = x_center;
    k = k + 1;
end

symOut = symOut(1:k-1);
end




Rs   = 500e3;                     % Symbol rate
ppm  = 50e-6;                     % 50 ppm tolerance
M    = 4;                         % QPSK
preamble_len = 200;               % Length of preamble

% -------- Frame structure --------
N_total = 20e4;                    % Total number of data symbols
N_frm   = 20e4;                   % Symbols per frame
num_frm = N_total / N_frm;        % Number of frames

% Generate all random bits once
bitStream = randi([0 1], 2*N_total, 1);

% Map bits → symbol indices
symbolIndex = bi2de(reshape(bitStream, 2, []).', 'left-msb');

% QPSK modulation
txSymbolStream = pskmod(symbolIndex, M, pi/4, 'gray');

% Preamble (fixed alternating sequence)
syncPreamble = (1/sqrt(2)) * (1+1j) * (-1).^(0:preamble_len-1).';
preamble_indices = mod((0:preamble_len-1), 2).' * 3;

Fs_nom = 2*Rs;                    % Nominal sampling rate
alpha_r = 1;                    % Random clock offset factor
Fs = Fs_nom * (1 + ppm*(2*alpha_r-1));   % Actual sampling rate

sps = Fs_nom / Rs;                % Samples per symbol

% RRC pulse shaping
rolloff = 0.5;
span = 10;
rrcPulse = rcosdesign(rolloff, span, sps, 'sqrt');
halfFilterLen = (length(rrcPulse)-1)/2;

EbNo_dB = 0:1:10;
SER_sim = zeros(size(EbNo_dB));
SNR_dB = EbNo_dB + 10*log10(2) - 10*log10(sps);

% Carrier offset range
cfo_max = 1e9*ppm;
% cfo_actual = cfo_max * (2*0.6 - 1);
cfo_actual = 0;

for p = 1:length(EbNo_dB)

    symErrCount = 0;
    disp(EbNo_dB(p))

    for frm = 1:num_frm

        % -------- TX build frame --------
        dataIdx = (frm-1)*N_frm + (1:N_frm);
        frameData = txSymbolStream(dataIdx);

        frameSymbols = [syncPreamble ; frameData];

        % Upsample + pulse shape
        txUpsampled = upsample(frameSymbols, sps);
        txBaseband = conv(txUpsampled, rrcPulse);
        txBaseband = txBaseband(halfFilterLen+1:end-halfFilterLen);

        % TX/RX clocks
        txTime = (0:length(txBaseband)-1).' / Fs_nom;
        rxTime = (0:ceil(length(txBaseband)*Fs/Fs_nom)-1).' / Fs;

        rxClockSignal = interp1(txTime, txBaseband, rxTime, 'linear', 0);

        % Apply carrier offset
        rxCarrierShifted = rxClockSignal .* exp(1j*2*pi*cfo_actual*rxTime);

        % Noise 
        % rxNoisy = awgn(rxCarrierShifted, SNR_dB(p), 'measured');
        % rxNoisy = awgn(rxCarrierShifted, 8, 'measured');
        rxNoisy = rxCarrierShifted;

        % Random sample drops
        maxDropSamples = 15;
        numDropped = randi([0 maxDropSamples]);
        % numDropped = 0;
        rxDropped = rxNoisy(numDropped+1:end);

        % -------- COARSE CFO ----------
        Ns_pre = preamble_len * sps;
        rxPre = rxDropped(1:Ns_pre);

        numSegments = 3;
        segLength = floor(length(rxPre) / numSegments);

        freqStep = Fs / (160 * sps);
        freqComb = (-2*cfo_max : freqStep : 2*cfo_max).';

        powerAvg = zeros(length(freqComb),1);

        for seg = 1:numSegments
            idx = (seg-1)*segLength + (1:segLength);

            x = rxPre(idx).^2;
            w = hamming(segLength);
            x = x .* w;

            n = (0:segLength-1).';

            for kk = 1:length(freqComb)
                testTone = exp(-1j*2*pi*freqComb(kk)*n/Fs);
                powerAvg(kk) = powerAvg(kk) + abs(sum(x .* testTone)).^2;
            end
        end

        powerAvg = powerAvg / numSegments;
        [~, kmax] = max(powerAvg);
        f2_est = freqComb(kmax);
        cfo_est_coarse = f2_est / 2;

        % -------- Apply coarse CFO correction --------
        totalSamples = (N_frm + preamble_len)*sps;
        droppedCount = totalSamples - length(rxDropped);

        n_all = (droppedCount : droppedCount+length(rxDropped)-1).';

        rxAfterCFO = rxDropped .* exp(-1j*2*pi*cfo_est_coarse*n_all/Fs);

        % -------- COARSE TIMING --------
        rxPreCFO = rxAfterCFO(1:Ns_pre);

        timingMetric = zeros(sps,1);
        for tau = 1:sps
            timingMetric(tau) = sum(abs(rxPreCFO(tau:sps:end)).^2);
        end
        [~, tau_hat] = max(timingMetric);

        % -------- MATCHED FILTER --------
        rxMF = conv(rxAfterCFO, rrcPulse);
        rxMF = rxMF(halfFilterLen+1:end-halfFilterLen);

        % symbsync = comm.SymbolSynchronizer(SamplesPerSymbol=sps,TimingErrorDetector='Early-Late (non-data-aided)');
        % [symbols,terr] = symbsync(rxMF);
        % [sym, terr] = el_sync(rxMF, sps);
        clckppm = ppm*(2*alpha_r-1)*1000000;
        fname = sprintf('clock_phase_%d.png', clckppm);
        [sym, terr, er] = el_simple(rxMF, sps,tau_hat);
        rxClockRecovered = sym;
        figure;
        subplot(2,1,1);
        plot(terr); grid on;
        title('Fractional Timing Offset');
        xlabel('Symbols');
        ylabel('Phase normalized by sps');
        legend(['' num2str(clckppm) ' ppm']);
        saveas(gcf,fname);

        % n = (0:length(rxMF)-1).';
        % mu = tau_hat-1;
        % alpha = 0.01;
        % 
        % maxSymbols = floor((length(rxMF) - sps) / sps);
        % rxClockRecovered = zeros(maxSymbols, 1);
        % 
        % mu_history = [];
        % e_history = [];
        % sps_e = sps * (1 + ppm*(2*alpha_r-1));
        % ted_acc=0;
        % k = 1;
        % while (mu + sps/8) < length(rxMF)-1
        %     x_center = interp1(n, rxMF, mu,        'linear', 0);
        %     x_early  = interp1(n, rxMF, mu-sps/8,  'linear', 0);
        %     x_late   = interp1(n, rxMF, mu+sps/8,  'linear', 0);
        % 
        % 
        %     % --- early-late timing detector ---
        %     ted = real(x_center)*(real(x_late)-real(x_early)) + ...
        %           imag(x_center)*(imag(x_late)-imag(x_early));
        % 
        %     ted_acc = ted_acc+ted;
        %     % --- advance pointer smoothly ---
        %     if mod(k,10)==0                   
        %         mu = mu + sps;
        %         ted_acc=0;
        %     else
        %         mu = mu +sps;
        %     end
        % 
        %         % --- store histories ---
        %     mu_history = [mu_history; mu];
        %     e_history  = [e_history; ted];
        % 
        %     rxClockRecovered(k) = x_center;
        %     k = k + 1;
        % end
        % e_history = movmean(e_history,500);
        % 
        % rxClockRecovered = rxClockRecovered(1:k-1);
        
        % 
        % figure('Name', 'Timing Recovery Diagnostics');
        % % 
        % subplot(2,1,1);
        % plot(mod(mu_history,sps)); grid on;
        % title('Fractional Timing Offset (should converge)');
        % ylabel('μ mod sps');
        % 
        % subplot(2,1,2);
        % plot(e_history); grid on;
        % title('unwrap');
        % ylabel('μ ');

        % -------- FINE CFO / PHASE TRACKING --------
        z = rxClockRecovered(:);
        Ns = length(z);

        symDropEst = N_frm + preamble_len - Ns;

        Lseg = 10;
        Np_use = preamble_len - Lseg;
        kvec = (0:Ns-1).';

        phi_blk = [];
        k_blk = [];

        for start = 1:Lseg:Np_use
            idx = start : start + Lseg - 1;
            phi = angle(z(idx) .* conj(syncPreamble(idx+symDropEst)));
            phi_blk(end+1,1) = mean(unwrap(phi));
            k_blk(end+1,1)   = mean(kvec(idx));
        end

        fpoly = polyfit(k_blk, phi_blk, 1);

        dphi_hat = fpoly(1);
        phi0_hat = fpoly(2);

        residual_freq = dphi_hat*Rs/2*pi;


        z_corr = z .* exp(-1j * (phi0_hat + dphi_hat * kvec));

        phi_track = 0;

        for start = preamble_len+1 : Lseg : Ns-Lseg+1
            idx = start : start + Lseg - 1;

            s_hat = pskmod(pskdemod(z_corr(idx), M, pi/4,'gray'), ...
                           M, pi/4,'gray');

            phi_err = angle(z_corr(idx) .* conj(s_hat));
            dphi = mean(unwrap(phi_err));

            phi_track = phi_track + dphi;

            z_corr(idx(1):end) = z_corr(idx(1):end) .* exp(-1j*dphi);
        end

        rxFineSymbols = z_corr;

        % -------- DECISION --------
        detectedData = rxFineSymbols(end-N_frm+1:end);

        rx_idx = pskdemod(detectedData, M, pi/4, 'gray');
        tx_idx = symbolIndex(dataIdx);

        symErrCount = symErrCount + sum(rx_idx ~= tx_idx);

    end

    SER_sim(p) = symErrCount / N_total;

end

EbNo_lin = 10.^(EbNo_dB/10);
SER_theory = 2*qfunc(sqrt(2*EbNo_lin)) - qfunc(sqrt(2*EbNo_lin)).^2;

figure;
semilogy(EbNo_dB, SER_sim, 'o-', 'LineWidth', 1.5); hold on;
semilogy(EbNo_dB, SER_theory, '--', 'LineWidth', 1.5);
grid on;
xlabel('E_b / N_0 (dB)');
ylabel('SER');
legend('40ppm clk, 30ppm cfo', 'Ideal');
title('QPSK SER vs Eb/N0');
saveas(gcf,'clk_carr_offset_ser.png');
