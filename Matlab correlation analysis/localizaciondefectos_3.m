% function [pentax,pentay,heptax,heptay]=localizaciondefectos(U)
%%%prgrama para calcular la posicion de las esferas y calcular la posicion
%%%de los defectos 5 y 7. Los parametros iniciales son la matriz del parametro de orden y los datos son la
%%%posicion del los hepta (heptax, heptay),y los penta (pentax,
%%%pentay).1/10/2009.

clear all

%Files to travel
archivos_a_leer = {%U
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\1\U.dat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\2\U.dat',       
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\3\U.dat',
	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\4\U.dat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\5\U.dat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\6\U.dat', 
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\7\U.dat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\8\U.dat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\9\U.dat',
	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\10\U.dat'
       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\11\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\12\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\13\U.dat',
% 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\14\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\15\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\16\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\17\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\18\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\19\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\20\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\21\U.dat', 
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\22\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\23\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\24\U.dat',
%   	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\25\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\26\U.dat'
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\27\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\28\U.dat',
% 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\29\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\30\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\31\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\32\U.dat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\33\U.dat'
 }

 archivos_a_salvar = {%RANGUP
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\1\x1.mat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\2\x2.mat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\3\x3.mat',
	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\4\x4.mat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\5\x5.mat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\6\x6.mat', 
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\7\x7.mat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\8\x8.mat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\9\x9.mat',
	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\10\x10.mat',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\11\x11.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\12\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\13\x1.mat',
% 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\14\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\15\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\16\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\17\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\18\x1.mat',
% 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\19\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\20\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\21\x1.mat', 
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\22\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\23\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\24\x1.mat',
% 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\25\x1.mat'
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\26\x1.mat'
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\27\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\28\x1.mat',
% 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\29\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\30\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 15\31\x.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 15\32\x1.mat',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 15\33\x1.mat'
}

 archivos_figura_1 = {
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig1.tif', 
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig2.tif',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig3.tif',
 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig4.tif',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig5.tif',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig6.tif', 
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig7.tif',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig8.tif',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig9.tif',
	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig10.tif'       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig12.tif',
      'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 14\fig11.tif',
% 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig14.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig15.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig16.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig17.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig18.tif',
% 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig19.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig20.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig21.tif', 
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig22.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig23.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig24.tif',
% 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig25.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig26.tif'
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig27.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig28.tif',
% 	  'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig29.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 20\fig30.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 15\31\fig.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 15\32\fig.tif',
%       'D:\Mis documentos\Ana\Proyecto\proyect 3D\CDS RUN\Alineados\orden + esferas\MATLAB\PRUEBA 15\33\fig.tif'       
}

%Get the lenght of the list
cantidad_archivos_a_leer   = length(archivos_a_leer);
cantidad_archivos_a_salvar = length(archivos_a_salvar);

%if cantidad_archivos_a_leer == cantidad_archivos_a_salvar
    cantidad_archivos = cantidad_archivos_a_leer;
%else
%    cantidad_archivos = 0;
%end
  
for archivo_actual = 1:cantidad_archivos;
    
    %Get the current file path
    arch1 = archivos_a_leer{archivo_actual}; %cargo a U
    load(arch1);

    U=real(U);
    %dim=256;%size(U);
    m=512%dim(1,1);
    n=512%dim(1,2);
    s=1;
    dp=-0.4%-0.1

    for i=2:m-1
        for j=2:10
            if U(i,j)<dp
            if U(i-1,j)>U(i,j)
                   if U(i-1,j+1)>U(i,j)
                          if U(i-1,j-1)>U(i,j)
                if U(i+1,j)>U(i,j)
                       if U(i+1,j+1)>U(i,j)
                           if U(i+1,j-1)>U(i,j)
                    if U(i,j-1)>U(i,j)
                        if U(i,j+1)>U(i,j) 
                            X(s,1)=j;
                        X(s,2)=i;
                    s=s+1;    
                    end
                    end
                        end
                    end
                end
            end
        end
    end
    end
    end
    end
    %%%%%%%%%%%%
    for i=2:m-1
        for j=n-10:n-1
            if U(i,j)<dp
            if U(i-1,j)>U(i,j)
                   if U(i-1,j+1)>U(i,j)
                          if U(i-1,j-1)>U(i,j)
                if U(i+1,j)>U(i,j)
                       if U(i+1,j+1)>U(i,j)
                           if U(i+1,j-1)>U(i,j)
                    if U(i,j-1)>U(i,j)
                        if U(i,j+1)>U(i,j) 
                            X(s,1)=j;
                        X(s,2)=i;
                    s=s+1;    
                    end
                    end
                        end
                    end
                end
            end
        end
    end
    end
    end
    end
    %%%%%%%%%%%%%%%%%
    for i=2:10
        for j=11:m-11
            if U(i,j)<dp
            if U(i-1,j)>U(i,j)
                   if U(i-1,j+1)>U(i,j)
                          if U(i-1,j-1)>U(i,j)
                if U(i+1,j)>U(i,j)
                       if U(i+1,j+1)>U(i,j)
                           if U(i+1,j-1)>U(i,j)
                    if U(i,j-1)>U(i,j)
                        if U(i,j+1)>U(i,j) 
                         X(s,1)=j;
                        X(s,2)=i;
                    s=s+1;    
                    end
                    end
                        end
                    end
                end
            end
        end
    end
    end
    end
    end
    %%%%%%%%%%%%%%%%%%%
    for i=n-10:n-1
        for j=11:m-11
            if U(i,j)<dp
            if U(i-1,j)>U(i,j)
                   if U(i-1,j+1)>U(i,j)
                          if U(i-1,j-1)>U(i,j)
                if U(i+1,j)>U(i,j)
                       if U(i+1,j+1)>U(i,j)
                           if U(i+1,j-1)>U(i,j)
                    if U(i,j-1)>U(i,j)
                        if U(i,j+1)>U(i,j) 
                         X(s,1)=j;
                        X(s,2)=i;
                    s=s+1;    
                    end
                    end
                        end
                    end
                end
            end
        end
    end
    end
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    primero=s;
    %%%%%%%%%datos del cubo central
    for i=11:n-11
        for j=11:m-11
            if U(i,j)<dp
            if U(i-1,j)>U(i,j)
                   if U(i-1,j+1)>U(i,j)
                          if U(i-1,j-1)>U(i,j)
                if U(i+1,j)>U(i,j)
                       if U(i+1,j+1)>U(i,j)
                           if U(i+1,j-1)>U(i,j)
                    if U(i,j-1)>U(i,j)
                        if U(i,j+1)>U(i,j) 
                          X(s,1)=j;
                        X(s,2)=i;
                    s=s+1;    
                    end
                    end
                        end
                    end
                end
            end
        end
    end
    end
    end
    end
    ultimo=s-1;
    [V, C] = voronoin(X);
     pent=1;
        hept=1;
    for r=primero:ultimo
        A=size(C{r});
        if A(1,2)==5
            pentax(pent)=X(r,1);
            pentay(pent)=X(r,2);
             pent=pent+1;
        end
            if A(1,2)==7
           heptax(hept)=X(r,1);
           heptay(hept)=X(r,2);
           hept=hept+1;  
           end
    end

     arch = archivos_a_salvar{archivo_actual};
     save(arch,'X');

    handle = figure;
    pcolor(U), shading interp, colormap(gray), ...
    axis('off'), axis('equal')
    hold on

    plot(X(:,1),X(:,2),'* w')
    plot(pentax,pentay,'* r')
    plot(heptax,heptay,'* g')
    
    
     archivo_figura_1 = archivos_figura_1{archivo_actual};
     saveas(handle, archivo_figura_1)

    %print(handle, archivo_figura_1, '-dpng','-r600', '-opengl');
   
  %  saveas(fig,archivo_actual.png)     
    clear X
    clear A
    clear pentax
    clear pentay
    clear heptax
    clear heptay
end
