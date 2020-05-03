% PART_1

% This script allows you to plot NACA-4-digit airfoil by entering the value
% of the 4 digits.
clear
close all


% 1) Define NACA airfoil

NACA=input('insert NACA foil to analyse _ NACA');
AoA=input('enter the value of AoA in degrees _ AoA=');
x=linspace(0,1,500);

if NACA<=1000 % SYM NACA airfoil
    
    y_upper=5*(NACA/100)*(0.2969*x.^(0.5)-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4-0.0021*x.^5)';     %% COMMENT: I've changed added 0.0021*x.^5 in order to close the airfoil at the TE
    y_lower=-5*(NACA/100)*(0.2969*x.^(0.5)-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4-0.0021*x.^5)';
    y_c=0;
    foil_name="00"+num2str(NACA);
else % NON-SYM NACA airfoil
    
    foil=num2str(NACA);
        out=sscanf(foil, '%1f%1f%2f',3);
        foil_name=num2str(NACA);
        
        a=out(1,1)/100;
        b=out(2,1)/10;
        
        y_t=5*(out(3,1)/100)*(0.2969*x.^(0.5)-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4-0.0021*x.^5)';
       
        y_c=zeros(length(x),1);
        theta=zeros(length(x),1);
        
    for n=1:length(x)
        
    if x(n)<=b  
        
     theta(n)=atan((2*a/(b^2))*(b-x(n)));
     y_c(n)=(a/(b^2))*(2*b*x(n)-x(n)^2);
     
    elseif x(n)>b 
        
         theta(n)=atan((2*a/((1-b)^2))*((1-2*b)+2*b.*x(n)-x(n)^2));
         y_c(n)=(a/((1-b)^2))*((1-2*b)+2*b*x(n)-x(n)^2);
         
    end
     
  
     
    end
    
    y_upper=y_c+y_t.*cos(theta);
     y_lower=y_c-y_t.*cos(theta);
    
     
     
end

% rotate NACA airfoil plot

x_r_up=(x'-0.25).*cos(AoA*pi/180)+y_upper.*sin(AoA*pi/180)+0.25;         %
x_r_low=(x'-0.25).*cos(AoA*pi/180)+y_lower.*sin(AoA*pi/180)+0.25;        % Rotating the airfoil around (Xc,Yc)=(0.25,0)
y_r_up=-(x'-0.25).*sin(AoA*pi/180)+y_upper.*cos(AoA*pi/180);             %
y_r_low=-(x'-0.25).*sin(AoA*pi/180)+y_lower.*cos(AoA*pi/180);            %
x_r_c=(x'-0.25).*cos(AoA*pi/180)+y_c.*sin(AoA*pi/180)+0.25;              % Rotating the MCL
y_r_c=-(x'-0.25).*sin(AoA*pi/180)+y_c.*cos(AoA*pi/180);                  %


hold on
axis equal
plot(x_r_up,y_r_up,'b-')        % plotting upper surface
plot(x_r_low,y_r_low,'b-')      % plotting lower surface
plot(x_r_c,y_r_c,'g-')          % plotting MCL
grid on
title(['AoA=',num2str(AoA),'º'])
legend(strcat('NACA',foil_name))

% inBetween = [y_r_low, fliplr(y_r_up)];            % add those rows to
% fill(x, inBetween, 'c');                          % shade the airfoil

hold off

% ---------------------------------------------------------------- %

% building panels

aux_y_r_low=flipud(y_r_low);                                      % using aux_y_r_low as an auxiliary vector to build Airfoil_coord
aux_x_r_low=flipud(x_r_low);


Airfoil_coord=zeros(length(y_r_up)+length(y_r_low)-1,2);          % Airfoil_coord=(x,y)
                                                                  % coordinates of airfoil surface "counted" clockwise starting from (x,y)=(0,0)
Airfoil_coord(1:length(y_r_up),2)=y_r_up;
Airfoil_coord(1:length(y_r_up),1)=x_r_up;
Airfoil_coord(length(y_r_up)+1:length(y_r_up)+length(y_r_low)-1,2)=aux_y_r_low(2:length(y_r_low));
Airfoil_coord(length(y_r_up)+1:length(y_r_up)+length(y_r_low)-1,1)=aux_x_r_low(2:length(x_r_low));

% panel inclination of each airfoil panel
%
% panel inclination [rad] counted counter-clockwise from x-direction
panel_inclination=zeros(length(y_r_up)+length(y_r_low)-2,1);
n_y=zeros(length(y_r_up)+length(y_r_low)-2,1);
n_x=zeros(length(y_r_up)+length(y_r_low)-2,1);
tau_y=zeros(length(y_r_up)+length(y_r_low)-2,1);
tau_x=zeros(length(y_r_up)+length(y_r_low)-2,1);
mid_panel_x=zeros(length(y_r_up)+length(y_r_low)-2,1);
mid_panel_y=zeros(length(y_r_up)+length(y_r_low)-2,1);
for i=1:length(y_r_up)+length(y_r_low)-2
    
    panel_inclination(i)=atan((Airfoil_coord(i+1,2)-Airfoil_coord(i,2))/(Airfoil_coord(i+1,1)-Airfoil_coord(i,1)));
    % define normal and tangential direction of each panel 
    n_y(i)=(Airfoil_coord(i+1,1)-Airfoil_coord(i,1))/(((Airfoil_coord(i,1)-Airfoil_coord(i+1,1))^2+(Airfoil_coord(i+1,2)-Airfoil_coord(i,2))^2)^(0.5));
    n_x(i)=(Airfoil_coord(i,2)-Airfoil_coord(i+1,2))/(((Airfoil_coord(i,1)-Airfoil_coord(i+1,1))^2+(Airfoil_coord(i+1,2)-Airfoil_coord(i,2))^2)^(0.5));
    tau_x(i)=-(Airfoil_coord(i,1)-Airfoil_coord(i+1,1))/(((Airfoil_coord(i,1)-Airfoil_coord(i+1,1))^2+(Airfoil_coord(i+1,2)-Airfoil_coord(i,2))^2)^(0.5));
    tau_y(i)=(Airfoil_coord(i+1,2)-Airfoil_coord(i,2))/(((Airfoil_coord(i,1)-Airfoil_coord(i+1,1))^2+(Airfoil_coord(i+1,2)-Airfoil_coord(i,2))^2)^(0.5));
    mid_panel_x(i)=(Airfoil_coord(i,1)+Airfoil_coord(i+1,1))*0.5;
    mid_panel_y(i)=(Airfoil_coord(i,2)+Airfoil_coord(i+1,2))*0.5;
end
panel_inclination(length(y_r_up):end,1)=panel_inclination(length(y_r_up):end,1)+pi;
normal_tan=[n_x,n_y,tau_x,tau_y];
mid_panel_coord=[mid_panel_x,mid_panel_y];
