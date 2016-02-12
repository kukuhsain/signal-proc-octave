%clean everything
clc all; clear all; close all;

%inisialisasi_1
fo = 1000000;
to = 1/fo;

t = 0.01*to:0.01*to:10*to;
pembawa = cos(2*pi*fo*t);
data = [1 0 1 1 0 0 1 0 1 0];

%mengubah nilai [1 0] menjadi [1 -1]
for m = 1:10
	informasi(m) = data(m);
	if (informasi(m)==0)
		informasi(m)=-1;
	endif
endfor

%inisialisasi_2
k = ones(1,100);
masuk = zeros(1,1000);

%inisialisasi sinyal masukan, agar dimensi vektornya setara
for m = 1:10
	for n = 1:100 
		o = (m*100+n)-100;
		masuk(o) = informasi(m) * k(n);
	endfor
endfor

%perkalian untuk mendapatkan sinyal bpsk
for m = 1:1000
	bpsk(m) = masuk(m)*pembawa(m);
endfor 

%penampilan gambar keluaran
figure(1)
plot(t,pembawa)

figure(2)
plot(t,masuk)

figure(3)
plot(t,bpsk)

%penerima menerima sinyal bpsk, kemudian dikalikan dengan sinyal pembawa
for m = 1:1000
	terima(m) = bpsk(m)*pembawa(m);
endfor 

%hasil perkalian selanjutnya diintegralkan
%setelah diintegralkan, hasilnya di thresholding, untuk mendapat data biner
for m = 1:10
	for n = 1:100 
		o = (m*100+n)-100;		
		if (n != 1)		
			itemp(n) = itemp(n-1) + terima(o);		
		else itemp(n) = terima(o);
		endif
		integral(o) = itemp(n)/100;
		if (integral(o)>=0)
			output(o) = 1;
		else output(o) = 0;
		endif
	endfor
endfor

figure(4)
plot(t,terima)

figure(5)
plot(t,integral)

figure(6)
plot(t,output)
