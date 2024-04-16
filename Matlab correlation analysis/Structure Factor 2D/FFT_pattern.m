%programa para calcular el espectro de Fourier

clear all

archivos_a_leer = {
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\1\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\2\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\3\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\4\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\5\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\6\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\7\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\8\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\9\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\10\U.dat',
	'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\11\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\12\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\13\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\14\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\15\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\16\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\17\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\18\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\19\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\20\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\21\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\22\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\23\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\24\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\25\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\26\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\27\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\28\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\29\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\30\U.dat',
	'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\31\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\32\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\33\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\34\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\35\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\36\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\37\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\38\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\39\U.dat',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\40\U.dat'
    };
	
archivos_a_salvar2 = {
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\1\FFT_01',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\2\FFT_02',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\3\FFT_03',     
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\4\FFT_04',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\5\FFT_05',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\6\FFT_06',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\7\FFT_07',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\8\FFT_08',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\9\FFT_09',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\10FFT_10',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\11\FFT_11',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\12\FFT_12',     
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\13\FFT_13',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\14\FFT_14',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\15\FFT_15',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\16\FFT_16',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\17\FFT_17',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\18\FFT_18',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\19\FFT_19',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\20\FFT_20'  
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\21\FFT_21',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\22\FFT_22',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\23\FFT_23',     
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\24\FFT_24',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\25\FFT_25',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\26\FFT_26',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\27\FFT_27',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\28\FFT_28',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\29\FFT_29',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\30\FFT_30',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\31\FFT_31',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\32\FFT_32',     
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\33\FFT_33',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\34\FFT_34',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\35\FFT_35',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\36\FFT_36',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\37\FFT_37',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\38\FFT_38',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\39\FFT_39',
     'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\40\FFT_40'  
    };  

archivos_figura_1 = {%Orientational map
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT1.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT3.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT4.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT5.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT6.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT7.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT8.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT9.tif',
	'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT10.tif',
	'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT11.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT12.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT13.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT14.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT15.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT16.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT17.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT18.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT19.tif',
	'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT20.tif',        
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT21.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT22.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT23.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT24.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT25.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT26.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT27.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT28.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT29.tif',
	'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT30.tif',
	'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT31.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT32.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT33.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT34.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT35.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT36.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT37.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT38.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT39.tif',
	'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT40.tif'       
    };

archivos_figura_2 = {%defectos positivos y negativos
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_01.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_02.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_03.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_04.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_05.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_06.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_07.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_08.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_09.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_10.tif',
	'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_11.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_12.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_13.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_14.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_15.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_16.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_17.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_18.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_19.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_20.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_21.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_22.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_23.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_24.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_25.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_26.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_27.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_28.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_29.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_30.tif',
	'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_31.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_32.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_33.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_34.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_35.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_36.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_37.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_38.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_39.tif',
    'C:\Users\Anabella\Documents\Visual Studio 2012\Projects\CDS\CDS\ANALISIS DEFECTOS MATLAB2\FFT2_40.tif'
    };


	
%Get the lenght of the list
cantidad_archivos_a_leer   = length(archivos_a_leer);
cantidad_archivos_a_salvar2 = length(archivos_a_salvar2);

if cantidad_archivos_a_leer == cantidad_archivos_a_salvar2
    cantidad_archivos = cantidad_archivos_a_leer;
else
    cantidad_archivos = 0;
end


for archivo_actual = 1:cantidad_archivos;
    
    %Get the current file path
    arch1 = archivos_a_leer{archivo_actual};
    load(arch1);
    
UF=fft2(U);
UFS = fftshift(UF);
RUFS = real(UFS);

R=120;
P = zeros(R+1,1); 
j=1;

for r=0:R
    for theta=0:2*pi
        
        x = round ((128 + r * cos(theta)));
        y = round ((128 + r * sin(theta)));
        P(j) = P(j) + abs(UFS(x, y));
    end
    
    j=j+1;  
end

 FFT = 100*log(1+abs(UFS));
   
      arch = archivos_a_salvar2{archivo_actual};
      save(arch, 'P');

    handle = figure;
	plot(P); colormap(gray); 
    title('FFT');
    hold on
	axis('off'), axis('equal')
    set(gcf, 'Visible', 'off');
		
	archivo_figura_1 = archivos_figura_1{archivo_actual};
	saveas(handle, archivo_figura_1)

    %Plot defectos positivos y negativos
	handle2 = figure;
	imagesc(real(FFT)); colormap(gray); 
    title('magnitude spectrum'); 
	axis('off'), axis('equal')
    axis tight
    axis 'auto y'
    axis 'auto x'
	hold on
    set(gcf, 'Visible', 'off');
	
	archivo_figura_2 = archivos_figura_2{archivo_actual};
	saveas(handle2, archivo_figura_2)
 

%     figure;
%     imagesc(angle(UFS));  colormap(gray);
%     title('phase spectrum');


 clearvars -except ColorMap mtotal archivo_actual cantidad_archivos archivos_a_leer archivos_a_salvar archivos_figura_1 archivos_a_salvar2 archivos_a_leer2 archivos_figura_2 

end %End for each loop

