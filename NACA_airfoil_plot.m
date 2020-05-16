% NACA airfoil plot

% This script allows you to plot NACA-4-digit airfoil by entering the value
% of the 4 digits.

clear
clc
close all
format long


prompt = {'NACA airfoil','Angle of attack in degrees'};         % 
dlg_title = 'NACA airfoil design';                              % Creating an input-dialog window
num_lines = 1;                                                  %           
defaultans = {'0000','0'};                                      %   
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);       %

NACA = str2double(answer{1,1}) ;                                % Coverting user's answer
AoA = str2double(answer{2,1}) ;                                 % in double precision number
    
x = linspace(0,1,500);

if NACA<=1000 % SYM NACA airfoil
    
    y_upper=5*(NACA/100)*(0.2969*x.^(0.5)-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4-0.0021*x.^5)';     %% COMMENT: I've added 0.0021*x.^5 in order to close the airfoil at the TE
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
plot(x_r_c,y_r_c,'g-')          % plotting MCL
plot(x_r_low,y_r_low,'b-')      % plotting lower surface

grid on
title(['AoA=',num2str(AoA),'º'])
legend(strcat('NACA',foil_name),'MCL')

%  inBetween = [y_r_low, fliplr(y_r_up)];            % add those rows to
%  fill(x, inBetween, 'y');                          % shade the airfoil

hold off

% ---------------------------------------------------------------- %

% clear unuseful variables
clear a b defaultans dlg_title num_lines prompt foil ...
        n answer out