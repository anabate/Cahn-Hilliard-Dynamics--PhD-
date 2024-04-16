%programa para calcular <h(x)h(x+r)> 20/8/2011

clear all

%Files to travel
archivos_a_leer = {
      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\01\U.dat',
      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\02\U.dat',
      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\07\U.dat',
      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\15\U.dat',
      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\30\U.dat',
      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\40\U.dat'
    };

archivos_a_salvar2 = {
    'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\01\Rang_01',
    'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\02\Rang_02',
    'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\07\Rang_07',
    'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\15\Rang_15',
    'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\30\Rang_30',
    'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\40\Rang_40'
    };

archivos_a_salvar = {
     'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\01\correlation_01',
     'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\02\correlation_02',
     'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\07\correlation_07',
     'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\15\correlation_15',     
     'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\30\correlation_30',
     'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\2 CAPAS\Matlab\40\correlation_40'
    };


%Get the lenght of the list
cantidad_archivos_a_leer   = length(archivos_a_leer);
cantidad_archivos_a_salvar = length(archivos_a_salvar);

if cantidad_archivos_a_leer == cantidad_archivos_a_salvar
    cantidad_archivos = cantidad_archivos_a_leer;
else
    cantidad_archivos = 0;
end
  
for archivo_actual = 1:cantidad_archivos;
    
    %Get the current file path
    arch1 = archivos_a_leer{archivo_actual};
    load(arch1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%CALCULO ORIENTACION%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

del1=1
del2=1
mx=1024
my=1024

%%%Derivada respecto a X
for n=1:mx;
    for m=1:my;
        A1(n,m)=(i/del1)*(sin(2*pi*(n-1)/mx));
    end
end

%%%Derivada respecto a y
for n=1:mx;
    for m=1:my;
        A2(n,m)=(i/del2)*(sin(2*pi*(m-1)/my));
    end
end

%%%%matriz de promedio
for n=1:mx;
    for m=1:my;
        Aprom(n,m)=(1/4)*(cos(2*pi*(n-1)/mx)+cos(2*pi*(m-1)/my))+(1/2);
    end
end

%%%%%%%%%%%%%%%%%%%% extraccion de nx,ny
UF=fft2(U);
nx=ifft2(UF.*A1);
ny=ifft2(UF.*A2);
modul=sqrt(nx.^2+ny.^2);
nx=nx./modul;
ny=ny./modul;
cosX=2*nx.^2-1;
senX=2*nx.*ny;

%%%promedio
for h=1:10
senp=ifft2(fft2(senX).*Aprom);
senX=real(senp);
end

for h=1:10
cosp=ifft2(fft2(cosX).*Aprom);
cosX=real(cosp);
end

%%%determinacion del angulo
ang=atan2(senX,cosX);
Rang = real(ang);

 arch = archivos_a_salvar2{archivo_actual};
 save(arch,'Rang');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%CALCULO LONGITUD DE CORRELACIÓN%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    R=150;
    H=fft2(Rang);
    tam=1024;

    for n=1:1024;
        for m=1:1024;
            A(n,m)=exp(i*2*pi*(n-1)/1024);
        end
    end
    for n=1:1024;
        for m=1:1024;
            B(n,m)=exp(i*2*pi*(m-1)/1024);
        end
    end

      
    for r = 0:R
        C = zeros(1024);
        for x=0:r
            y=round(sqrt(r.^2-x.^2));
            C=(A.^x).*(B.^y) + (A.^x).*(B.^-y) + (A.^-x).*(B.^y) + (A.^-x).*(B.^-y) + C;
        end
        RR=ifft2((1/(4*(r+1)))*C.*H);
        
        HH(r+1)=sum(sum(RR.*Rang));
        
        disp(r);
        
       arch = archivos_a_salvar{archivo_actual};
       save(arch, 'HH');
    end

end %End foreach loop
        
    