close all
clear all
clc

%%% OFDM in GNU Octave %%%
%%% by Kukuh Nomikusain %%%
%%% 11/319359/TK/38488 %%%

%%%%% INISIALISASI

fo = 24000000;
to = 1/fo;

t = 0.01*to:0.01*to:8*to;
carr_real = cos(2*pi*fo*t);
carr_imaj = sin(2*pi*fo*t);

%%%%% INPUT
% input data biner
% inp = [0010 0110 1110 1010 0011 0111 1111 1011];

inp = [0 0 1 0 0 1 1 0 1 1 1 0 1 0 1 0 0 0 1 1 0 1 1 1 1 1 1 1 1 0 1 1]

figure(1)
stem(1:32,inp)

%%%%% MAPPING DAN S/P
% mapping 16-QAM
% map(hasil_map, pos)
% pos=1 => real , pos=2 => imajiner

l = 0;

for k=1:4:length(inp)
	temp_real = [inp(k) inp(k+1)];
	temp_imaj = [inp(k+2) inp(k+3)];
	
	h = k+l-l*4;
	
	if temp_real == [0 0]
		map(h,1) = -3;
	elseif temp_real == [0 1]
		map(h,1) = -1;
	elseif temp_real == [1 1]
		map(h,1) = 1;
	else
		map(h,1) = 3;
	endif

	if temp_imaj == [0 0]
		map(h,2) = -3;
	elseif temp_imaj == [0 1]
		map(h,2) = -1;
	elseif temp_imaj == [1 1]
		map(h,2) = 1;
	else
		map(h,2) = 3;
	endif

	l = l+1;
endfor

mapping = map(:,1) + i*map(:,2);

figure(2)
stem(1:8,map(:,1))
figure(3)
stem(1:8,map(:,2))

%%%%% IFFT
ifftmap = ifft(mapping);

figure(4)
stem(1:8,abs(ifftmap))

%%%%% INTERPOLASI
inter_real = real(ifftmap)';
inter_imaj = imag(ifftmap)';

k = ones(1,100);
interpol_real = zeros(1,800);
interpol_imaj = zeros(1,800);

% real
for m = 1:8
	for n = 1:100 
		o = (m*100+n)-100;
		interpol_real(o) = inter_real(m) * k(n);
	endfor
endfor

% imajiner
for m = 1:8
	for n = 1:100 
		o = (m*100+n)-100;
		interpol_imaj(o) = inter_imaj(m) * k(n);
	endfor
endfor

figure(5)
stem(1:800,interpol_real)
figure(6)
stem(1:800,interpol_imaj)

%%%%% PERKALIAN SEBELUM TRANSMIT
% perkalian dengan carrier
channel_real = carr_real .* interpol_real;
channel_imaj = carr_imaj .* interpol_imaj;

%%%%% TRANSMIT
% channel
channel = channel_real + channel_imaj;

figure(7)
plot(1:800,channel_real)
figure(8)
plot(1:800,channel_imaj)
figure(9)
plot(1:800,channel)


%%%%% RECEIVE
% receiver real
rec_real = channel .* carr_real;
figure(10)
plot(1:800,rec_real)

% receiver imajiner
rec_imaj = channel .* carr_imaj;
figure(11)
plot(1:800,rec_real)


% integral
% integral real
for m = 1:8
	o = m*100-99;		
	integral_real(m) = sum(rec_real(o:o+99))*2/100;
endfor
figure(12)
stem(1:8,integral_real)

% integral imajiner
for m = 1:8
	o = m*100-99;		
	integral_imaj(m) = sum(rec_imaj(o:o+99))*2/100;
endfor
figure(13)
stem(1:8,integral_imaj)

integral = (integral_real' + i*integral_imaj');

%%%%% FFT
rec_fft = round(fft(integral));

figure(14)
stem(1:8,abs(rec_fft))


%%%%% DEMAPPING DAN P/S
% demapping 16-QAM

for k=1:length(rec_fft)
	temp_real = real(rec_fft(k));
	temp_imaj = imag(rec_fft(k));
	
	h = k*4-3;
	
	if temp_real == -3
		out(h:h+1) = [0 0];
	elseif temp_real == -1
		out(h:h+1) = [0 1];
	elseif temp_real == 1
		out(h:h+1) = [1 1];
	else
		out(h:h+1) = [1 0];
	endif

	if temp_imaj == -3
		out(h+2:h+3) = [0 0];
	elseif temp_imaj == -1
		out(h+2:h+3) = [0 1];
	elseif temp_imaj == 1
		out(h+2:h+3) = [1 1];
	else
		out(h+2:h+3) = [1 0];
	endif
endfor

%%%%% OUTPUT
out

figure(15)
stem(1:32,out)
