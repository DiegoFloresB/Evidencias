%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%++++ Equipo 6: Mariachis del Amazonas| Modelo del volcán El Reventador ++++
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%El siguiente programa modela las trayectorias de un vestigio cuando son
%eruptados del volcán considerando la resistencia al aire.

%CONSIDERACIONES-----------------------------------------------------------
%No se descrbibirán funciones o líneas de código que ya se hayan descrito
%anteriormente. Es decir, la programación de los tres proyectiles tienen
%elementos similares por lo que solo se agregará una pequeña nota que
%mencione el funcionamiento de dicha línea de código.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%++++++++++++++++++++++++++++ C Ó D I G O +++++++++++++++++++++++++++++++++
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Borrar figuras, valores y Command Window----------------------------------
clear all
close all
clc

%Formulas Esfera-----------------------------------------------------------

m=67000; %Masa del vestigio [kg]
r=10; %Radio del vestigio [m]
area_sph=(4*pi*r)/2; %[m^2]


%**************************************************************************
%----------------Primer proyectil: Trayectoria pre-establecida--------------
%**************************************************************************


v0=500;         %Velocidad Inicial [m/s]
angle=75.77;    %Ángulo de disparo [°]


%-----------------------------------------------NOTA:DATOS GENERALES A USAR


y0=3560;                    %Altura del volcán [m]
cd_sphere=0.47;             %Constante
Temperature1=0;             %Temperatura en condiciones normales[C°]
Temperature2=600;           %Temperatura en erupción [°C]
g=9.81;                     %Gravedad [m/s^2]

%Densidad del aire a 3500 msnm [kg/m^3]
air_d_Reventador_nc=(348.42*(1-y0*1.05*10^-4))/(Temperature1+273);

%Densidad del aire considerando la temperatura del aire [kg/m^3] %No se usa
air_d_Reventadort_temp=(348.42*(1-y0*1.05*10^-4))/(273+Temperature2);


%Fórmulas para calcular la fuerza de arrastre [FD]-------------------------

%Fuerza de Arrastre [N]
Fd=((1/2)*cd_sphere*air_d_Reventador_nc*area_sph*v0^2);

ax=(-Fd*cosd(angle))/m;     %Aceleración en "x" considerando FD [m/s^2]
ay=-g-(Fd*sind(angle)/m);   %Aceleración en "y" considerando FD y g [m/s^2]
time_v=((-2*v0*sind(angle))/ay)+y0;     %Tiempo de vuelo [s]


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Gráficar 1ra Parábola-----------------------------------------------------
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt=input('Elige el paso temporal (dt) '); %Tamaño del paso en el tiempo
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t=0:dt:time_v;                   %Intervalo de tiempo                         
n=length(t);                     %Numero de posiciones dentro del intervalo    

y=[];                            %Lista Y
y(1)=3560;                       %Posición Inicial Y(1)
a=3560-v0*sind(angle)*dt;        %Posición anterior a Y(1)
y(2)=2*y(1)-a+ay*dt^2;           %Calcular Y(2)
 
vxf=v0*cosd(angle)+ax*t;
vyf=v0*sind(angle)+ay*t;
magnitud= sqrt(vxf.^2+vyf.^2);
ang=atand(vyf/vxf);

%Metodo de Verlet
for c=3:n                          
    y(c)=2*y(c-1)-y(c-2)+ay*dt^2;   %Calcular Y(c)
    if y(c)<0                     
        y(c)=0;
        magnitud(c)=0;

    end
end
    
hmax=(v0.^2)*(sind(angle).^2)/(2*-ay);
thmax=(v0)*sind(angle)/-ay;
R= (v0.^2)*(sind(2*angle))/-ay;

x=[];
x(1)=0;
b=0-v0*cosd(angle)*dt;
x(2)=2*x(1)-b+ax*dt^2;
x(2)=2*x(1)+ax*dt^2;

for c=3:n
    x(c)=2*x(c-1)-x(c-2)-+ax*dt^2;
    if y(c)==0
        x(c)=x(c-1);
    end
end
    

%**************************************************************************
%-------------------Segundo proyectil: Elegir medidas ---------------------
%**************************************************************************

v01=input("¿Qué velocidad inicial desea? ");
angle1=input("¿Qué ángulo de disparo desea? ");

Fd1=((1/2)*cd_sphere*air_d_Reventador_nc*area_sph*v01^2);

ax1=(-Fd1*cosd(angle1))/m;    
ay1=-g-(Fd1*sind(angle1)/m);  
time_v1=((-2*v01*sind(angle1))/ay1)+y0;     

t1=0:dt:time_v1;                  
n1=length(t1);                     

y1=[];                             
y1(1)=3560;                        
a1=3560-v01*sind(angle1)*dt;       
y1(2)=2*y1(1)-a1+ay1*dt^2;       

vxf1=v01*cosd(angle1)+ax1*t1;
vyf1=v01*sind(angle1)+ay1*t1;
magnitud1= sqrt(vxf1.^2+vyf1.^2);
ang1=atand(vyf1/vxf1);
        
for c1=3:n1
    y1(c1)=2*y1(c1-1)-y1(c1-2)+ay1*dt^2;
    if y1(c1)<0
        y1(c1)=0;
        magnitud1(c1)=0;

    end
end

hmax1=(v01.^2)*(sind(angle1).^2)/(2*-ay1);
thmax1=(v01)*sind(angle1)/-ay1;
R1=(v01.^2)*(sind(2*angle1))/-ay1;

x1=[];
x1(1)=0;
b1=0-v01*cosd(angle1)*dt;
x1(2)=2*x1(1)-b1+ax1*dt^2;
%x(2)=2*x(1)+ax*dt^2;

for c1=3:n1
    x1(c1)=2*x1(c1-1)-x1(c1-2)-+ax1*dt^2;
    if y1(c1)==0
        x1(c1)=x1(c1-1);
    end
end
    

   

%**************************************************************************
%-----------------Tercer proyectil: Medidas Aleatorias---------------------
%**************************************************************************

v02=(500-100).*rand(1,1)+100;
angle2=(90-10).*rand(1,1)+10;
Fd2=((1/2)*cd_sphere*air_d_Reventador_nc*area_sph*v02^2);

ax2=(-Fd2*cosd(angle2))/m;    
ay2=-g-(Fd2*sind(angle2)/m);  
time_v2=((-2*v02*sind(angle2))/ay2)+y0;     

t2=0:dt:time_v2;                  
n2=length(t2);                     

y2=[];                             
y2(1)=3560;                        
a2=3560-v02*sind(angle2)*dt;       
y2(2)=2*y2(1)-a2+ay2*dt^2;         

vxf2=v02*cosd(angle2)+ax2*t2;
vyf2=v02*sind(angle2)+ay2*t2;
magnitud2= sqrt(vxf2.^2+vyf2.^2);
ang2=atand(vyf2/vxf2);
        
for c2=3:n2
    y2(c2)=2*y2(c2-1)-y2(c2-2)+ay2*dt^2;
    if y2(c2)<0
        y2(c2)=0;
        magnitud2(c2)=0;
    end
end

hmax2=(v02.^2)*(sind(angle2).^2)/(2*-ay2);
thmax2=(v02)*sind(angle2)/-ay2;
R2=(v02.^2)*(sind(2*angle2))/-ay2;

x2=[];
x2(1)=0;
b2=0-v02*cosd(angle2)*dt;
x2(2)=2*x2(1)-b2+ax2*dt^2;
%x(2)=2*x(1)+ax*dt^2;

for c2=3:n2
    x2(c2)=2*x2(c2-1)-x2(c2-2)-+ax2*dt^2;
    if y2(c2)==0
        x2(c2)=x2(c2-1);
end
end
    


%**************************************************************************
%-------------------------------IMAGENES-----------------------------------
%**************************************************************************

%Llamamos las imágenes png. 
[reventador,~,alpha] = imread('volcan2.png');          %Inserta imagen del volcán
[roca_0,~,alpha1] = imread("eu pedrini manito.png");   %Inserta imagen de la piedra
[roca_1,~,alpha2] = imread("eu pedrini manito.png");
[roca_2,~,alpha3] = imread("eu pedrini manito.png");


%Posiciones----------------------------------------------------------------

%%Rotamos los vectores a columnas para posteriormente añadir una tabla(x,y)
%con información

%1ra Trayectoria
posicion_x=x';
posicion_y=y';

%2da Trayectoria
posicion_x1=x1';
posicion_y1=y1';

%3ra Trayectoria
posicion_x2=x2';
posicion_y2=y2';

%Graficar 1er proyectil----------------------------------------------------
%Creamos un ciclo en el que se reemplazarán los datos en cada iteración

for c=1:length(x)
    plot(x,y);
    grid on
    hold on
    plot(x1,y1)
    plot(x2,y2)

    %Graficamos la parábola completa, para que cada vez que se reinicie el
    %ciclo ésta se vuelva a graficar
    plot(x,y);
    grid on
    hold on

    %Usamos más espacio del necesario para los límites de la grafica
    xlim([-20000,20000]);
    ylim([0,10000]);


    %Imprimimos las imágenes png en cada posición de la parábola
    image(reventador,'alphadata',im2double(alpha),'xdata',[8500 -8500],'ydata',[4500 -5])
    hold on

    %Escribimos el texto que se irá actualizando en cada iteración
    formatSpec = "Trayectoria Establecida (Morada) \nTiempo = %.2f\nX = %.2f\nY = %.2f\nVelocidad = %.2f.";
    formatSpec1 = "Trayectoria Dada (Naranja) \nTiempo = %.2f\nX = %.2f\nY = %.2f.\nVelocidad = %.2f";
    formatSpec2 = "Trayectoria Aleatoria (Amarilla) \nTiempo = %.2f\nX = %.2f\nY = %.2f\nVelocidad = %.2f.";

    %Definimos la cadena de carácteres
    txt=sprintf(formatSpec,t(c),x(c),y(c),magnitud(c));
    txt1=sprintf(formatSpec1,t1(c),x1(c),y1(c),magnitud1(c));
    txt2=sprintf(formatSpec2,t2(c),x2(c),y2(c),magnitud2(c));

    %Mandamos a imprimir nuestro texto
    text(6000,9000,txt)
    text(6000,7000,txt1)
    text(6000,5000,txt2)

    image(roca_0,'alphadata',im2double(alpha1),'xdata',[x(c)+500 x(c)-500],'ydata',[y(c) y(c)-300]);
    image(roca_1,'alphadata',im2double(alpha2),'xdata',[x1(c)+500 x1(c)-500],'ydata',[y1(c) y1(c)-300]);
    image(roca_2,'alphadata',im2double(alpha3),'xdata',[x2(c)+500 x2(c)-500],'ydata',[y2(c) y2(c)-300]);


    %Stop-motion:Cada cuanto se actualizan los elementos (imágenes y texto)
    pause(0.001)
    hold off
    
end



%Definimos tabla (x,y) hasta este punto para que se añada después de
%terminar el ciclo "for"
T= table(posicion_x,posicion_y);
T1= table(posicion_x1,posicion_y1);
T2= table(posicion_x2,posicion_y2);


%Nombramos las columnas
T.Properties.VariableNames = ["Posición en X","Posición en Y"];
T1.Properties.VariableNames = ["Posición en X","Posición en Y"];
T2.Properties.VariableNames = ["Posición en X","Posición en Y"];

%Formato de la tabla

uitable("Data",T{:,:},"ColumnName",T.Properties.VariableNames,"RowName","");
uitable("Data",T1{:,:},"ColumnName",T1.Properties.VariableNames,"RowName","");
uitable("Data",T2{:,:},"ColumnName",T2.Properties.VariableNames,"RowName","");

disp(" ")
disp("Proyectil 1")
disp("La altura máxima es "+hmax);
disp("El tiempo en alcanzar la altura máxima es "+thmax);
disp("El alcance horizontal es "+R);
disp(" ")
disp("Proyectil 2")
disp("La altura máxima es "+hmax1);
disp("El tiempo en alcanzar la altura máxima es "+thmax1);
disp("El alcance horizontal es "+R1);
disp(" ")
disp("Proyectil 3")
disp("La altura máxima es "+hmax2);
disp("El tiempo en alcanzar la altura máxima es "+thmax2);
disp("El alcance horizontal es "+R2);