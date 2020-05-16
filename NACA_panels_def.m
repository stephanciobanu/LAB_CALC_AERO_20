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
%     n_y(i)=(Airfoil_coord(i+1,1)-Airfoil_coord(i,1))/(((Airfoil_coord(i,1)-Airfoil_coord(i+1,1))^2+(Airfoil_coord(i+1,2)-Airfoil_coord(i,2))^2)^(0.5));
%     n_x(i)=(Airfoil_coord(i,2)-Airfoil_coord(i+1,2))/(((Airfoil_coord(i,1)-Airfoil_coord(i+1,1))^2+(Airfoil_coord(i+1,2)-Airfoil_coord(i,2))^2)^(0.5));
    tau_x(i)=-(Airfoil_coord(i,1)-Airfoil_coord(i+1,1))/(((Airfoil_coord(i,1)-Airfoil_coord(i+1,1))^2+(Airfoil_coord(i+1,2)-Airfoil_coord(i,2))^2)^(0.5));
    tau_y(i)=(Airfoil_coord(i+1,2)-Airfoil_coord(i,2))/(((Airfoil_coord(i,1)-Airfoil_coord(i+1,1))^2+(Airfoil_coord(i+1,2)-Airfoil_coord(i,2))^2)^(0.5));
    
    mid_panel_x(i)=(Airfoil_coord(i,1)+Airfoil_coord(i+1,1))*0.5;
    mid_panel_y(i)=(Airfoil_coord(i,2)+Airfoil_coord(i+1,2))*0.5;
end
n_x(:)=-tau_y(:);
n_y(:)=tau_x(:);
panel_inclination(length(y_r_up):end,1)=panel_inclination(length(y_r_up):end,1)+pi;
normal_tan=[n_x,n_y,tau_x,tau_y];
mid_panel_coord=[mid_panel_x,mid_panel_y];

% clear useless variables
clear i