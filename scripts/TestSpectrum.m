T = 80;
Fs = 2*pi;

y = sin([1:2000]*2*pi/T);
yfilt = conv_band_pass(y, [3 100]);
%yfilt2 = conv_band_pass(y, [6 10]);

[b,a] = butter(2, sort([2*pi/100 2*pi/3]/(Fs/2)), 'bandpass');
yfilt = filter(b,a,y);
%Hd = designfilt('lowpassfir','FilterOrder',20,'CutoffFrequency',2*pi/150, ...
%                'DesignMethod','window','Window',{@butter,3},'SampleRate',Fs);

subplot(211);
PlotSpectrum(y);
PlotSpectrum(yfilt);
%PlotSpectrum(yfilt2);
linex(1/T);

subplot(212);
plot(y); hold on;
plot(yfilt);
%plot(yfilt2);
xlim([100 200]);

%%

rng default

Fs = 1;
n = 1:365;
x = cos(2*pi*(1/7)*n)+cos(2*pi*(1/30)*n-pi/4);
trend = 3*sin(2*pi*(1/1480)*n);
y = x+trend+0.5*randn(size(n));

[pxx,f] = periodogram(y,[],length(y),Fs);

subplot(2,1,1)
plot(n,y)
xlim([1 365])
xlabel('Days')
grid

subplot(2,1,2)
plot(f,10*log10(pxx))
xlabel('Cycles/day')
ylabel('dB')
grid
