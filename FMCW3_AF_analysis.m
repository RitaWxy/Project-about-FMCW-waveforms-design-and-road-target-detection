clc;
clear all; 
close all;

c=3e8; %(speed light)

%%Trget information

Target_range_resolution=0.5;    %(in meters)
max_target_speed=230;  %(m/s)
max_unambiguos_range=150;  %(in meters)
%%Signal information


B=c/(2*Target_range_resolution);   %(bamdwidth in Hz)




Tr=2*max_unambiguos_range/c;       %Time repetition interval
T=6e-4;         %duration of the whole signal
Np = T/Tr; %number of pulses of the whole signal

%% frequency and time axis definition

%Generating a single chirp
c=3e8; 
f_s=2*B;     %%sample frequency
dt=1/f_s;
t=0:dt:Tr-dt;
df=1/Tr;
f=0:df:(f_s-df);
Ns=length(t);     %Length of single chirp


mu=2*pi*B/Tr;               %%Ramp - with 2pi factor
ramp_frequency=(mu/2*pi).*t;

figure(1);
plot(t,ramp_frequency);
grid on;
title('Single chirp')
xlabel('Time,s')
ylabel('Freq, Hz')

train_ramp_ferquency=repmat(ramp_frequency,1,round(T/Tr));
t1=0:dt:T-dt;

figure(2);
plot(t1,train_ramp_ferquency);
grid on;
title('Single chirp')
xlabel('Time,s')
ylabel('Freq, Hz')

s=exp(1i*(mu/2)*t.^2);    %%complex transmit signal


figure(3);
plot(t,real(s));
grid on;
title('Up-Chirp')
xlabel('Time,s')
ylabel('Amplitude, norm')


S=fft(s);

figure(4); 
plot(f-f_s/2,abs(fftshift(S))); 
grid on;
title('FFT of Chirp')
xlabel('Freq,Hz')
ylabel('Amplitude, norm')


R0=0+0*t;           %%change of distance of scatterer point target

tau0=2*R0/c; 
Sr = S.*exp(-1j*2*pi.*f.*tau0);
Sro= Sr.*conj(S);    %signal after matched filter


r=t*c/2;

figure(6);
plot(r,abs(ifft((Sro))));
%plot(r,fftshift(abs(ifft((Sro)))));
grid on
title('conversion from beat freq - Received signal after mixer');
xlabel('Range [m]'); ylabel('Single Recieved Signal');


Sc=repmat(Sr,1,Np);   %%vector with returnof Np pulses
Scmatrix=reshape(Sc,Ns,Np);   %matrix with phase-shift due to scatterer point-slike target at distace R0

Smf=conj(S);   
HRR=zeros(Ns,Np);

for i=1:Np
  Y=Scmatrix(:,i).*Smf.';

HRR(:,i)=fftshift(ifft(Y));


end

% 
% figure(7);
% mesh(abs(HRR));
Im = HRR;

for i= 1:Ns
    %Im(i,:) = fftshift(fft( HRR(i,:)));
    
    Im(i,:) = fft( HRR(i,:));
end

 
 PRI=size(Im,1)*dt;
PRF=1/PRI;
 N=size(Im,2);
fk=(-1/2:1/N:(1/2-1/N))*PRF;

figure(8);
mesh(abs(Im));

NFFT=2^10;


IM1=fft2(Im);

Im1=ifft2(IM1,NFFT,NFFT);


rzp=linspace(0,150,1024);
fkzp=linspace(-PRF/2,PRF/2,NFFT);

AF=fftshift(abs(Im1)/max(abs(Im1(:))));  %Ambiguity function

figure(10); 
surf(fkzp,rzp,10*log10(fftshift(abs(Im1)/max(abs(Im1(:))))));
shading interp
colormap jet;
xlabel('Doppler [Hz]'); ylabel('Range [m]');
caxis([-40 0]); 
colorbar;
ylim([0 100]);
title('High Resolution Map');

PositionTargetRange=max(AF);

[Value1,IndexRange]=max(PositionTargetRange);
PositionTargetVelocity=max(AF');
[Value1,IndexVel]=max(PositionTargetVelocity);
% Place=find(abs(Matrix2)==PositionTarget;
[IndexRange IndexVel]
figure(11)

plot(fkzp,10.*log10((abs(AF(IndexVel,:)))))
grid on;
title('zeroRange cut')
xlabel('Doppler, Hz'); ylabel('Norm Amplitude');

figure(12)

plot(rzp,10.*log10((abs(AF(:,IndexRange)))))
grid on;

title('zeroDoppler cut')
xlabel('Range,m'); ylabel('Norm Amplitude');

 xlim([0 100])


