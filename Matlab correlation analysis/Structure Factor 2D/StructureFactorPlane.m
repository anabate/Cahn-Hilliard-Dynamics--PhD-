
% This scrip is used to get the structure factor of a plane pattern
% Write the input file and the structure factor is given as Sk
% A=h;%load('Out1.dat');
clear all

arch = 'C:\Users\Anabella\Desktop\EXP MARTIN\BUSCAR CENTROS Y g6\all .dat\U3.dat';
load(arch);

U = real(U3);

FTA=abs(fftshift(fft2(U)));
Sk=rscan(FTA);

Sk2 = transpose(Sk);

 archivo_destino = 'C:\Users\Anabella\Desktop\EXP MARTIN\BUSCAR CENTROS Y g6\all .dat\Sk1.dat';  
 save(archivo_destino, 'Sk2', '-ascii');