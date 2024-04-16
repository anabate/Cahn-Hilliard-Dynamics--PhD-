%programa para calcular <h(x)h(x+r)> 20/8/2011

clear all

%Files to travel
archivos_a_leer = {
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\1\U.dat',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\2\U.dat',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\3\U.dat',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\5\U.dat',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\8\U.dat',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\10\U.dat',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\13\U.dat',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\16\U.dat',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\20\U.dat',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\25\U.dat'
   
    };

archivos_a_salvar2 = {
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\1\Rang_01',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\2\Rang_02',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\3\Rang_03',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\5\Rang_04',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\8\Rang_08',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\10\Rang_10',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\Rang_13',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\Rang_16',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\20\Rang_20',
    'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\20\Rang_25'
   
   
    };

archivos_a_salvar = {
     'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\1\correlation_01',
     'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\2\correlation_02',
     'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\3\correlation_03',
     'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\5\correlation_05',     
     'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\8\correlation_08',
     'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\10\correlation_10',
     'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\13\correlation_13',
     'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\16\correlation_16',
     'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\20\correlation_20',
     'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\Orden + desorden\Matlab6 poco anneling\20\correlation_25'

    
    };

% archivos_figura_1 = {%defectos positivos
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\01\Figura_1_01.tif',
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\02\Figura_1_02.tif',
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\07\Figura_1_07.tif',
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\15\Figura_1_15.tif',     
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\30\Figura_1_30.tif',
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\40\Figura_1_40.tif'
%     };
	
% archivos_figura_2 = {%defectos positivos y negativos
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\01\Figura_2_01.tif',
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\02\Figura_2_02.tif',
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\07\Figura_2_07.tif',
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\15\Figura_2_15.tif',     
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\30\Figura_2_30.tif',
%      'd:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\TAU=0.19 f=0.47 -1024-\1 CAPA\Matlab\40\Figura_2_40.tif'
%     };


	
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

	del1=1;
	del2=1;
	mx=512;
	my=512;

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
%%%%%%%%%%%%%%%%%%%%ENCUENTRO LOS DEFECTOS%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 	for i=2:mx-1
% 		for j=2:my-1
% 		  
% 		   dif=ang(i+1,j)-ang(i-1,j);
% 		   difvec=[dif,dif+2*pi,dif-2*pi];
% 		   [Y,I]=min(abs(difvec));   
% 		   PHAX(i-1,j-1)=difvec(I)/2;      
% 		   
% 		end
% 	end
% 
% 	%%%%%derivada y
% 	for i=2:mx-1
% 		for j=2:my-1
% 		   
% 		   dif=ang(i,j+1)-ang(i,j-1);
% 		   difvec=[dif,dif+2*pi,dif-2*pi];
% 		   [Y,I]=min(abs(difvec));   
% 		   PHAY(i-1,j-1)=difvec(I)/2;  
% 		   
% 		end
% 	end
% 
% 	PHA=PHAX.^2+PHAY.^2;
% 
% 	%%%%%%%%%%%%%%%%%%ahora encuentro los puntos que conforman el defecto
% 	%para eso determino que  > max(max(PHA))*0.38
% 	t=1;
% 	maximo=1;%;max(max(PHA))*0.3;
% 	for i=11:mx-13
% 		for j=11:my-13
% 		   
% 			if PHA(i,j)>=maximo
% 			defecto(t,1)=i;
% 			defecto(t,2)=j;
% 			t=t+1;   
% 			end
% 			
% 		end
% 	end
% 	%%%%determino el centro de cada defecto
% 	etiq=1
% 	t=1
% 	r=12
% 	clear A
% 	A=defecto;
% 
% 	while etiq ~=0
% 		j=1;
% 		k=1;
% 		b=1;
% 		dimA=size(A,1);
% 		
% 		%%hasta largo de A
% 		
% 		for n=1:dimA
% 			
% 			if sqrt((A(1,1)-A(n,1)).^2+(A(1,2)-A(n,2)).^2)<r
% 			 
% 				puntox(j)=A(n,1);
% 				puntoy(j)=A(n,2);
% 				j=j+1;
% 				
% 			else
% 				
% 				AA(k,1)=A(b,1);
% 				AA(k,2)=A(b,2);
% 				k=k+1; 
% 			
% 			end
% 	   
% 		b=b+1;
% 		end
% 
% 
% 		if k~=1
% 			clear A
% 			A=AA;
% 			clear AA
% 			etiq=1;
% 
% 		else
% 			etiq=0;
% 		end
% 
% 		%encuentro centros de defectos
% 		defectox(t)=sum(puntox)/length(puntox);
% 		defectoy(t)=sum(puntoy)/length(puntoy);
% 		t=t+1;
% 
% 		clear puntox
% 		clear puntoy
% 
% 	end
% 
% 	f=1
% 	g=1
% 
% 	%%%calculo de la carga topologica
% 	for m=1:length(defectox)
% 
% 		pos1=round(defectox(m));
% 		pos2=round(defectoy(m));
% 		
% 		for i=-5:5
% 			l1(i+6)=PHAX(pos1+i,pos2+5);
% 		end
% 
% 		for i=-5:5
% 			l2(i+6)=PHAX(pos1+i,pos2-5);
% 		end
% 
% 		for i=-5:5
% 			l3(i+6)=PHAX(pos1+5,pos2+i);
% 		end
% 
% 		for i=-5:5
% 			l4(i+6)=PHAX(pos1-5,pos2+i);
% 		end
% 
% 		carga=sum(l1)-sum(l2)-sum(l3)+sum(l4)
% 		
% 		%%%armo archivos con defectos positivos y negativos
% 		if carga >=0.2
% 			defectopos(f,1)=pos1;
% 			defectopos(f,2)=pos2;
% 		f=f+1;
% 		end
% 
% 		if carga <=-0.2
% 			defectoneg(g,1)=pos1;
% 			defectoneg(g,2)=pos2;
% 		g=g+1;
% 		end
% 
% 	end
% 
% 	%Plot defectos positivos y negativos
% 	handle = figure;
% 	pcolor(real(PHA)), shading interp, ...
% 	axis('off'), axis('equal')
% 	hold on
%     plot(defectopos(:,2),defectopos(:,1),'* w')
% 	plot(defectoneg(:,2),defectoneg(:,1),'* r')
% 	
% 	archivo_figura_2 = archivos_figura_2{archivo_actual};
% 	saveas(handle, archivo_figura_2)
% 	
% 	handle = figure;
% 	pcolor(real(PHA)), shading interp, ...
%     hold on
% 	axis('off'), axis('equal')
% 	plot(defectoy,defectox,'* w')
% 	
% 	archivo_figura_1 = archivos_figura_1{archivo_actual};
% 	saveas(handle, archivo_figura_1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%CALCULO LONGITUD DE CORRELACIÓN%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    R=120;
    H=fft2(Rang);
    tam=512;

    for n=1:mx;
        for m=1:my;
            A(n,m)=exp(i*2*pi*(n-1)/tam);
        end
    end
    for n=1:mx;
        for m=1:my;
            B(n,m)=exp(i*2*pi*(m-1)/tam);
        end
    end

      
    for r = 0:R
        C = zeros(tam);
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
    clear Rang
    clear HH
    clear A
    clear B
    clear C
    clear A1
    clear A2
    clear Aprom
end %End for each loop
        
    