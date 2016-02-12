close all
clear all
clc

% input data biner
% inp = [0010 0110 1110 1010; 0011 0111 1111 1011; 0001 0101 1101 1001; 0000 0100 1100 1000];

inp = [0 0 1 0 0 1 1 0 1 1 1 0 1 0 1 0 0 0 1 1 0 1 1 1 1 1 1 1 1 0 1 1 0 0 0 1 0 1 0 1 1 1 0 1 1 0 0 1 0 0 0 0 0 1 0 0 1 1 0 0 1 0 0 0];

%konstelasi = (reshape(inp,4,16))';

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


% ifft
ifft_real = ifft(map(:,1));
ifft_imaj = ifft(map(:,2));



% interpolasi
inter_real = ifft_real';
inter_imaj = ifft_imaj';

figure, plot(1:16,inter_real)


