clc
clear all
close all
%%
fc = 77e9;
c = 3e8;
lambda = c/fc;

%%
% The sweep time can be computed based on the time needed for the signal to
% travel the unambiguous maximum range. In general, for an FMCW radar
% system, the sweep time should be at least 5 to 6 times the round trip
% time. This example uses a factor of 5.5.

range_max = 200;
tau=2*range_max/c;
tm = 5.5*tau;

%%
% The sweep bandwidth can be determined according to the range resolution
% and the sweep slope is calculated using both sweep bandwidth and sweep
% time.

range_res = 0.5;
bw = c/(2*range_res);
sweep_slope = bw/tm;

% In this example, the beat frequency corresponding to the maximum range is
% given by

fr_max=2*range_max*sweep_slope/c;

%%
% In addition, the maximum speed of a traveling car is about 230 km/h.
% Hence the maximum Doppler shift and the maximum beat frequency can be
% computed as

v_max = 230*1000/3600;
fd_max = 2*v_max/lambda;

fb_max = fr_max+fd_max;


%%
% This example adopts a sample rate of the larger of twice the maximum beat
% frequency and the bandwidth.

fs = 2*max(fb_max,bw);

%%

% System parameters                     Value
% --------------------------------------------
% Operating frequency (GHz)  [fc]         77
% Maximum target range (m)   [rang_max]   200
% Range resolution (m)       [range_res]  0.5
% Maximum target speed (km/h)[v_max]      230
% Sweep time (microseconds)  [tm]         7.33
% Sweep bandwidth (MHz)      [bw]         300
% Maximum beat frequency (MHz) [fb_max]   27.30
% Sample rate (MHz)          [fs]         300


%% *****   Target  ******* %
Range = [0 80 160 60 110];
Vel = [0 10 70 30 15]./3.6;
%Target 1
InitialRange_1 = Range(1);
TarVel_1 = Vel(1);
freD_1=2.*TarVel_1/lambda;
%Target 2
InitialRange_2 = Range(2);
TarVel_2 = Vel(2);
freD_2=2.*TarVel_2/lambda;
%Target 3
InitialRange_3 = Range(3);
TarVel_3 = Vel(3);
freD_3=2.*TarVel_3/lambda;
%Target 4
InitialRange_4 = Range(4);
TarVel_4 = Vel(4);
freD_4=2.*TarVel_4/lambda;
%Target 5
InitialRange_5 = Range(5);
TarVel_5 = Vel(5);
freD_5=2.*TarVel_5/lambda;
%% FMCW waveforms
Nchirp = 100;
singleT = 0: 1/fs : (tm-1/fs);

NFFT = 2048;
%% *****   Up-chirp  ******* %
Matrix=[]; Matrix1=[]; Matrix2=[];
coef = 2*pi*sweep_slope*tm;

for n=1:Nchirp
    t = singleT + (n-1)*tm;
    % Target 1
    delRange_1 = TarVel_1*(n-1)*tm;
    delay_1 = 2*(InitialRange_1 +delRange_1)/c;
     % Target 2
    delRange_2 = TarVel_2*(n-1)*tm;
    delay_2 = 2*(InitialRange_2 +delRange_2)/c;
    % Target 3
    delRange_3 = TarVel_3*(n-1)*tm;
    delay_3 = 2*(InitialRange_3 +delRange_3)/c;
     % Target 4
    delRange_4 = TarVel_4*(n-1)*tm;
    delay_4 = 2*(InitialRange_4 +delRange_4)/c;
     % Target 5
    delRange_5 = TarVel_5*(n-1)*tm;
    delay_5 = 2*(InitialRange_5 +delRange_5)/c;
    % Transmitted signal
    phaseT = pi*sweep_slope.*t.^2 - (n-1)*coef.*t;
    Tx_up = exp (1i * phaseT);
    
    % Returned signal
    % Target 1
    phaseR_1 = pi*sweep_slope.*t.^2 + 2*pi*freD_1.*t - 2*pi*sweep_slope*delay_1.*t...
        -(n-1)*coef.*t;
    Rx_up_1 = exp (1i * phaseR_1);
     % Target 2
    phaseR_2 = pi*sweep_slope.*t.^2 + 2*pi*freD_2.*t - 2*pi*sweep_slope*delay_2.*t...
        -(n-1)*coef.*t;
    Rx_up_2 = exp (1i * phaseR_2);
      % Target 3
    phaseR_3 = pi*sweep_slope.*t.^2 + 2*pi*freD_3.*t - 2*pi*sweep_slope*delay_3.*t...
        -(n-1)*coef.*t;
    Rx_up_3 = exp (1i * phaseR_3);
     % Target 4
    phaseR_4 = pi*sweep_slope.*t.^2 + 2*pi*freD_4.*t - 2*pi*sweep_slope*delay_4.*t...
        -(n-1)*coef.*t;
    Rx_up_4 = exp (1i * phaseR_4);
     % Target 5
    phaseR_5 = pi*sweep_slope.*t.^2 + 2*pi*freD_5.*t - 2*pi*sweep_slope*delay_5.*t...
        -(n-1)*coef.*t;
    Rx_up_5 = exp (1i * phaseR_5);
    
    Rx_up = Rx_up_1 + Rx_up_2 + Rx_up_3 + Rx_up_4 + Rx_up_5;
   % IF signal
    IF = Tx_up.*conj(Rx_up);
    Matrix(:,n) = IF;
    Matrix1(:,n) = fft(IF,NFFT);
end

%% Second FFT
for k=1:NFFT
     Matrix2(k,:)=fftshift( fft( Matrix1(k,:)) );
 end


%% Waveform plot   
 figure(1)  % Transmitted signal
 subplot(2,1,1);
 plot(t,real(Tx_up),'b'); % transmitted signal - time
 xlabel("Time");
 ylabel("Amplitude");
 title("Transmitted signal -- up-chirp"); 
 subplot(2,1,2); % spectrogram
 spectrogram(Tx_up,256,250,256,fs,'yaxis');

 figure(2) % Returned signal
 subplot(2,1,1);
 plot(t,real(Rx_up_5),'b'); % received signal - time
 xlabel("Time");
 ylabel("Amplitude");
 title("Reveive signal -- up-chirp"); 
 subplot(2,1,2); % spectrogram
 spectrogram(Rx_up_5,256,250,256,fs,'yaxis');

 figure(3) % IF signal
 subplot(2,1,1);
 plot(t,real(IF),'b'); % time
 xlabel("Time");
 ylabel("Amplitude");
 title("IF signal"); 
 subplot(2,1,2); % spectrogram
 spectrogram(IF,256,250,256,fs,'yaxis');
 
 %% *********   Range-Doppler plot   ***********%%
 % range vector
 freqR=(1:NFFT).* fs/NFFT;
 rangeAxis = freqR*c/(2*sweep_slope);
 
 % velocity vector
 fm=1/tm;
 dfm=fm/Nchirp;
 freqV=(0:dfm : (fm-dfm))-(fm/2);
 velAxis= freqV * (lambda/2);
 
 figure
 s = surf(-velAxis ,rangeAxis , 10.*log10(squeeze(abs(Matrix2))) )
 ylim([0 range_max])
 xlabel("Velocity of Targets (m/s)");
 ylabel("Range of Targets (m)");
 zlabel("Amplitude")
 s.EdgeColor = 'none'
 colormap default
 