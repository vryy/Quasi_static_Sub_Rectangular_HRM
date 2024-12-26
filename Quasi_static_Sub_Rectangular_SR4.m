clc
clear
close
%FEMSTSF code applied for subrectangular tunnel under quasi-static loading
%in full slip condition - The case of reference tunnel SR4
%========================= INPUT PARAMETERS ===============================
% ------------------------Seismic parameters-------------------------------
a_max = 0.5; % the peak horizontal acceleration of earthquake event
shear_strain = 0.3837; % shear_strain calculated by a_max value, in percent)
% -------------------------------------------------------------------------
num_steps=50; % number iterative steps of the calculation
num_el = 360; % Total number of element in a lining ring
% ---------------------Tunnel Lining parameters----------------------------
Size_tunnel = 1.0; % (times) size of tunnel compared to reference case (SR4)
Tunnel_height = 7.2*Size_tunnel; % (m) height of the tunnel (7.2 is height of the reference tunnel SR4)
Tunnel_width = 9.7*Size_tunnel; % (m) width of tunnel (9.7 is width of the reference tunnel SR4)
thick=0.5; % thickness of lining (m)
E=35000;   % elastic modulus of concrete (MPa)
B=1;       % width of the lining ring in longitudinal direction (m)
A=B*thick; % Cross-section area of the lining in longitudinal direction
J=B*thick^3/12; % inertia moment of lining section
%======Establishing the Geometry of Sub-Rectangular tunnel case (SR4)======
    % Excavation radius of the tunnel lining (m)
    Rex1=9.95; % Radius of top and bottom parts, SR4 case
    Rex2=1;    % Radius of shoulder parts, SR4 case
    Rex3=5.35; % Radius of side walls, SR4 case
    % Center O:
    XO=0.0;
    YO=0.0;
    % Center O1:
    XO1=0.0;
    YO1=6.35*Size_tunnel;
    % Center O2:
    XO2=0.0;
    YO2=-6.35*Size_tunnel;
    % Center O3:
    XO3=0.5*Size_tunnel;
    YO3=0.0;
    % Center O4:
    XO4=-0.5*Size_tunnel;
    YO4=0.0;
    % Center O5:
    XO5=3.4*Size_tunnel;
    YO5=1.93*Size_tunnel;
    % Center O6:
    XO6=3.4*Size_tunnel;
    YO6=-1.93*Size_tunnel;
    % Center O7:
    XO7=-3.4*Size_tunnel;
    YO7=-1.93*Size_tunnel;
    % Center O8:
    XO8=-3.4*Size_tunnel;
    YO8=1.93*Size_tunnel;

    % Excavation radius of the tunnel lining (m)
    Rex1=Rex1*Size_tunnel-0.5 + thick; % Radius of top and bottom parts
    Rex2=Rex2*Size_tunnel-0.5 + thick;    % Radius of shoulder parts
    Rex3=Rex3*Size_tunnel-0.5 + thick; % Radius of side walls
    D1=2*Rex1;
    D2=2*Rex2;
    D3=2*Rex3;
    % radius of middle line of the tunnel lining (m)
    R1=Rex1 - thick/2;
    R2=Rex2 - thick/2;
    R3=Rex3 - thick/2;
    
    % Angle of tunnel sections corresponding to each radius Rex:
    Angle1 = round(22)      % Half of Angle range at top and bottom
    Angle2 = round(41)      % Angle range at shoulders
    Angle3 = round(54/2)    % Half of Angle range at side walls
      
    A0 = 1;              % First point at the tunnel bottom 
    A1 = A0 + Angle1 - 1 % R1_Bottom of tunnel
    A2 = A1 + Angle2     % R2_shoulder of tunnel
    A3 = A2 + Angle3     % R3_Sidewall of tunnel
    A4 = A3 + Angle3     % R3_Sidewall of tunnel
    A5 = A4 + Angle2     % R2_shoulder of tunnel
    A6 = A5 + Angle1     % R1_Top of tunnel
    A7 = A6 + Angle1     % R1_Top of tunnel
    A8 = A7 + Angle2     % R2_shoulder of tunnel
    A9 = A8 + Angle3     % R3_Sidewall of tunnel
    A10 = A9 + Angle3    % R3_Sidewall of tunnel
    A11 = A10 + Angle2   % R2_shoulder of tunnel
    A12 = A11 + Angle1   % R1_Bottom of tunnel 
    
    % Coordinators of point along tunnel boundary measured counter-clockwise from tunnel bottom
    for i=A0:A1
        x(i)=XO1    + R1*cos((i-90)/180*pi());
        y(i)=YO1    + R1*sin((i-90)/180*pi());
        num_ele1=    A1-A0;                    % number of element in section 1
        delta(i)=Angle1/(180/pi())/num_ele1;   % calculation of the angle between two elements
        Le(i)=2*R1*sin(delta(i)/2);            % calculation of the length of each element
        Rex(i)=R1;
    end
    for i=A1:A2
        x(i)=XO6    + R2*cos((i-90)/180*pi());
        y(i)=YO6    + R2*sin((i-90)/180*pi());
        num_ele2=    A2-A1;                     % number of element in section 2
        delta(i)=Angle2/(180/pi())/num_ele2;    % calculation of the angle between two elements
        Le(i)=2*R2*sin(delta(i)/2);             % calculation of the length of each element
        Rex(i)=R2;
    end
    for i=A2:A3
        x(i)=XO4    + R3*cos((i-90)/180*pi());
        y(i)=YO4    + R3*sin((i-90)/180*pi());
        num_ele3=    A3-A2;                     % number of element in section 3
        delta(i)=Angle3/(180/pi())/num_ele3;    % calculation of the angle between two elements
        Le(i)=2*R3*sin(delta(i)/2);             % calculation of the length of each element
        Rex(i)=R3;
    end
    for i=A3:A4
        x(i)=XO4    + R3*cos((i-90)/180*pi());
        y(i)=YO4    + R3*sin((i-90)/180*pi());
        num_ele4=    A4-A3;                    % number of element in section 4
        delta(i)=Angle3/(180/pi())/num_ele4;   % calculation of the angle between two elements
        Le(i)=2*R3*sin(delta(i)/2);            % calculation of the length of each element
        Rex(i)=R3;
    end
    for i=A4:A5
        x(i)=XO5    + R2*cos((i-90)/180*pi());
        y(i)=YO5    + R2*sin((i-90)/180*pi());
        num_ele5=    A5-A4;                    % number of element in section 5
        delta(i)=Angle2/(180/pi())/num_ele5;   % calculation of the angle between two elements
        Le(i)=2*R2*sin(delta(i)/2);            % calculation of the length of each element
        Rex(i)=R2;
    end
    for i=A5:A6
        x(i)=XO2    + R1*cos((i-90)/180*pi());
        y(i)=YO2    + R1*sin((i-90)/180*pi());
        num_ele6=    A6-A5;                    % number of element in section 6
        delta(i)=Angle1/(180/pi())/num_ele6;   % calculation of the angle between two elements
        Le(i)=2*R1*sin(delta(i)/2);            % calculation of the length of each element
        Rex(i)=R1;
    end
    for i=A6:A7
        x(i)=XO2    + R1*cos((i-90)/180*pi());
        y(i)=YO2    + R1*sin((i-90)/180*pi());
        num_ele7=    A7-A6;                    % number of element in section 7
        delta(i)=Angle1/(180/pi())/num_ele7;   % calculation of the angle between two elements
        Le(i)=2*R1*sin(delta(i)/2);            % calculation of the length of each element
        Rex(i)=R1;
    end
    for i=A7:A8
        x(i)=XO8    + R2*cos((i-90)/180*pi());
        y(i)=YO8    + R2*sin((i-90)/180*pi());
        num_ele8=    A8-A7;                    % number of element in section 8
        delta(i)=Angle2/(180/pi())/num_ele8;   % calculation of the angle between two elements
        Le(i)=2*R2*sin(delta(i)/2);            % calculation of the length of each element
        Rex(i)=R2;
    end
    for i=A8:A9
        x(i)=XO3    + R3*cos((i-90)/180*pi());
        y(i)=YO3    + R3*sin((i-90)/180*pi());
        num_ele9=    A9-A8;                    % number of element in section 9
        delta(i)=Angle3/(180/pi())/num_ele9;   % calculation of the angle between two elements
        Le(i)=2*R3*sin(delta(i)/2);            % calculation of the length of each element
        Rex(i)=R3;
    end
    for i=A9:A10
        x(i)=XO3    + R3*cos((i-90)/180*pi());
        y(i)=YO3    + R3*sin((i-90)/180*pi());
        num_ele10=    A10-A9;                   % number of element in section 10
        delta(i)=Angle3/(180/pi())/num_ele10;   % calculation of the angle between two elements
        Le(i)=2*R3*sin(delta(i)/2);             % calculation of the length of each element
        Rex(i)=R3;
    end
    for i=A10:A11
        x(i)=XO7    + R2*cos((i-90)/180*pi());
        y(i)=YO7    + R2*sin((i-90)/180*pi());
        num_ele11=    A11-A10;                  % number of element in section 11
        delta(i)=Angle2/(180/pi())/num_ele11;   % calculation of the angle between two elements
        Le(i)=2*R2*sin(delta(i)/2);             % calculation of the length of each element
        Rex(i)=R2;
    end
    for i=A11:A12
        x(i)=XO1    + R1*cos((i-90)/180*pi());
        y(i)=YO1    + R1*sin((i-90)/180*pi());
        num_ele12=    A12-A11;                  % number of element in section 12
        delta(i)=Angle1/(180/pi())/num_ele12;   % calculation of the angle between two elements
        Le(i)=2*R1*sin(delta(i)/2);             % calculation of the length of each element
        Rex(i)=R1;
    end
  
    for i=1:359
        if i==1
            beta(i)=delta(i)/2; % beta is the mean angle of each element (with reference to the vertical towards down)
            eps(1)=-(pi()/2); % eps is the angle of each node (not element) evaluated starting from x axis (in counter-clockwise)
            teta(1)=0; % teta is the angle of each node with reference to the vertical (towards down) which will be used for the reaction force at nodes
        else
            beta(i)=beta(i-1)+delta(i);
            eps(i)=(beta(i)+beta(i-1))/2-(pi()/2);
            teta(i)=teta(i-1)+delta(i);
        end
    end
    beta(360)=360*pi()/180-beta(1);
    eps(360)=(beta(360)+beta(360-1))/2-(pi()/2); 
%========================Properties of soil ===============================
for i=1:360
        % Properties of the Primary Ground ;       
        c(i) = 0/1000; % cohesion of ground (MPa);      
        fi(i) = 33/(180/pi()); % friction angle of ground (degrees)        
        Eground(i)= 100; % elastic modulus of ground (MPa)       
        muy(i) = 0.34; % Unit weight of ground (MN/m3)       
        gama_ground(i) = 0.018; % input parameter       
        G_ground(i)= Eground(i)/2/(1+muy(i)); % shear modulus of ground (MPa)
end
%======================FINISH OF INPUT PARAMETERS==========================
c_ = c(1);
fi_ = fi(1);
gama_ground_ = gama_ground(1);
%--------------------sign used for the application of loads----------------
sign_vertical_1_360=1.0;
sign_vertical_90_270=1.0;
sign_horizontal_1_360=-1.;
%----------------------------active loads ---------------------------------
    % vertical Loads at TOP half of tunnel
    Pv_l_1 = ((-250.1*thick^4+458.4*thick^3-281.86*thick^2+72.864*thick-7.6155)+(59.325*a_max^4-96.445*a_max^3+60.362*a_max^2-17.465*a_max+1.9939) + (29.241*Eground(i)^-0.41))*G_ground(i)*shear_strain/100;  
    Pv_l_2 = ((-103.85*thick^4+191.41*thick^3-119.06*thick^2+30.571*thick-2.9479)+(4.3*50*a_max/(50*a_max)^2.5-0.0344) + 2.3736*Eground(i)^-0.124)*G_ground(i)*shear_strain/100;
    Pv_l_3 = ((2.0302*thick^4-5.7045*thick^3+5.9116*thick^2-2.6744*thick+0.4451)+(-1*50*a_max/(50*a_max)^2.7+0.0042)+(7.65*Eground(i)/Eground(i)^1.5-1.54116))*G_ground(i)*shear_strain/100;
    Pv_r_1 = ((122.61*thick^4-202.9*thick^3+95.628*thick^2-11.077*thick-0.6517)+(4.7*50*a_max/(50*a_max)^2.5-0.0376)+ 29.241*Eground(i)^-0.41)*G_ground(i)*shear_strain/100;
    Pv_r_2 = ((-6.0999*thick^3+4.7769*thick^2+0.6707*thick-0.767)+(1.85*50*a_max/(50*a_max)^2.7-0.00777)+ 2.451*Eground(i)^-0.145)*G_ground(i)*shear_strain/100;
    Pv_r_3 = ((0.5458*thick^4-2.9313*thick^3+4.303*thick^2-2.4475*thick+0.4803)+(-1*50*a_max/(50*a_max)^2.7+0.0042)+(7.55*Eground(i)/Eground(i)^1.5-1.46749))*G_ground(i)*shear_strain/100;
    Pv_top_left1 = Pv_l_1; 
    Pv_top_left2 = -Pv_l_2;
    Pv_top_left3 = Pv_l_3;
    Pv_top_right1 = -Pv_r_1; 
    Pv_top_right2 = Pv_r_2;
    Pv_top_right3 = -Pv_r_3;
    Gradient_ver_top_left1 = Pv_top_left1/(0.7127*Size_tunnel);
    Gradient_ver_top_left2 = (Pv_top_left2 - Pv_top_left3)/(1.3440*Size_tunnel);
    Gradient_ver_top_left3 = Pv_top_left3/(2.5105*Size_tunnel);
    Gradient_ver_top_right1 = Pv_top_right1/(0.7127*Size_tunnel);
    Gradient_ver_top_right2 = (Pv_top_right2 - Pv_top_right3)/(1.3440*Size_tunnel);
    Gradient_ver_top_right3 = Pv_top_right3/(2.5105*Size_tunnel);
    for i=A3+1:A4
        pv(i) = (Gradient_ver_top_right1*(R3*(1-cos((i-90)/180*pi()))));
    end
    for i=A4+1:A4+20
        pv(i) = (Gradient_ver_top_right1*(0.5559*Size_tunnel + R2*(cos(27/180*pi())-cos((i-90)/180*pi()))));
    end
    for i=A4+21:A5
        pv(i) = (Pv_top_right2 - Gradient_ver_top_right2*(R2*(cos(48/180*pi())-cos((i-90)/180*pi()))));
    end
    for i=A5+1:A5+7
        pv(i) = (Pv_top_right2 - Gradient_ver_top_right2*(0.2209*Size_tunnel + R1*(cos(68/180*pi())-cos((i-90)/180*pi()))));
    end
    for i=A5+8:A6
        pv(i) = (Pv_top_right3 - Gradient_ver_top_right3*(R1*(cos(75/180*pi())-cos((i-90)/180*pi()))));
    end
    for i=A6:A6+15
        pv(i) = Gradient_ver_top_left3*(R1*(cos(90/180*pi())-cos((i-90)/180*pi())));
    end
    for i=A6+16:A7
        pv(i) = Pv_top_left3 + Gradient_ver_top_left2*(R1*(cos(105/180*pi())-cos((i-90)/180*pi())));
    end
    for i=A7+1:A7+20
        pv(i) = Pv_top_left3 + Gradient_ver_top_left2*(1.1231*Size_tunnel + R2*(cos(112/180*pi())-cos((i-90)/180*pi())));
    end
    for i=A7+21:A8
        pv(i) = (Pv_top_left1 - Gradient_ver_top_left1*(R2*(cos(133/180*pi())-cos((i-90)/180*pi()))));
    end
    for i=A8+1:A9
        pv(i) = (Pv_top_left1 - Gradient_ver_top_left1*(0.1568*Size_tunnel + R3*(cos(153/180*pi())-cos((i-90)/180*pi()))));
    end     
    % ---------------------------------------
    % vertical Loads at LOWER half of tunnel
    Pv_lower_left1 = -Pv_r_1; 
    Pv_lower_left2 = Pv_r_2;
    Pv_lower_left3 = -Pv_r_3;
    Pv_lower_right1 = Pv_l_1; 
    Pv_lower_right2 = -Pv_l_2;
    Pv_lower_right3 = Pv_l_3;
    Gradient_ver_lower_left1 = Pv_lower_left1/(0.7127*Size_tunnel);
    Gradient_ver_lower_left2 = (Pv_lower_left2 - Pv_lower_left3)/(1.3440*Size_tunnel);
    Gradient_ver_lower_left3 = Pv_lower_left3/(2.5105*Size_tunnel);
    Gradient_ver_lower_right1 = Pv_lower_right1/(0.7127*Size_tunnel);
    Gradient_ver_lower_right2 = (Pv_lower_right2 - Pv_lower_right3)/(1.3440*Size_tunnel);
    Gradient_ver_lower_right3 = Pv_lower_right3/(2.5105*Size_tunnel);
    for i=A0:A0+15
        pv(i) = Gradient_ver_lower_right3*(R1*cos((90-i)/180*pi()));
    end
    for i=A0+16:A1
        pv(i) = Pv_lower_right3 + Gradient_ver_lower_right2*(R1*(cos((90-i)/180*pi())-cos((90-15)/180*pi())));
    end
    for i=A1+1:A1+20
        pv(i) = Pv_lower_right3 + Gradient_ver_lower_right2*(1.1231*Size_tunnel + R2*(cos((90-i)/180*pi())-cos((90-22)/180*pi())));
    end
    for i=A1+21:A2
        pv(i) = (Pv_lower_right1 - Gradient_ver_lower_right1*(R2*(cos((90-i)/180*pi())-cos((90-43)/180*pi()))));
    end
        for i=A2+1:A3
        pv(i) = (Pv_lower_right1 - Gradient_ver_lower_right1*(0.1568*Size_tunnel + R3*(cos((90-i)/180*pi())-cos((90-63)/180*pi()))));
    end
    for i=A9+1:A10
        pv(i) = (Gradient_ver_lower_left1*(R3*((cos((i-90)/180*pi()))-(cos((270-90)/180*pi())))));
    end  
    for i=A10+1:A10+20
        pv(i) = (Gradient_ver_lower_left1*(0.5559*Size_tunnel + R2*((cos((i-90)/180*pi()))-(cos((297-90)/180*pi())))));
    end
    for i=A10+21:A11
        pv(i) = (Pv_lower_left2 - Gradient_ver_lower_left2*(R2*((cos((i-90)/180*pi()))-(cos((318-90)/180*pi())))));
    end
    for i=A11+1:A11+7
        pv(i) = (Pv_lower_left2 - Gradient_ver_lower_left2*(0.2209*Size_tunnel + R1*((cos((i-90)/180*pi()))-(cos((338-90)/180*pi())))));
    end 
    for i=A11+8:A12
        pv(i) = Pv_lower_left3 - Gradient_ver_lower_left3*(R1*((cos((i-90)/180*pi()))-(cos((345-90)/180*pi()))));
    end
    % ---------------------------------------
    % horizontal Loads at RIGHT side of tunnel
    Ph_rt_1 = ((-5.7523*thick^2+5.9451*thick-1.5353)+15.315*Eground(i)^-0.34)*G_ground(i)*shear_strain/100;
    Ph_rt_2 = ((1.424*thick^4+0.8914*thick^3-5.294*thick^2+3.8368*thick-0.795)+0.3632*Eground(i)^0.1688)*G_ground(i)*shear_strain/100;
    Ph_rt_3 = ((5.4297*thick^4-13.643*thick^3+12.774*thick^2-5.2909*thick+0.8195)-0.3632*Eground(i)^0.1688)*G_ground(i)*shear_strain/100;
    Ph_rb_1 = ((2.1236*thick - 1.0632)+15.315*Eground(i)^-0.34)*G_ground(i)*shear_strain/100;
    Ph_rb_2 = ((0.6164*thick - 0.3103)+0.3632*Eground(i)^0.1688)*G_ground(i)*shear_strain/100;
    Ph_rb_3 = ((5.4297*thick^4-13.643*thick^3+12.774*thick^2-5.2909*thick+0.8195)-0.3632*Eground(i)^0.1688)*G_ground(i)*shear_strain/100;
    Ph_right_top1 = -Ph_rt_1; 
    Ph_right_top2 = Ph_rt_2;
    Ph_right_top3 = -Ph_rt_3; 
    Ph_right_bot1 = Ph_rb_1; 
    Ph_right_bot2 = -Ph_rb_2;
    Ph_right_bot3 = Ph_rb_3;
    Gradient_hor_right_top1 = Ph_right_top1/(2.5234*Size_tunnel);
    Gradient_hor_right_top2 = (Ph_right_top2 - Ph_right_top3)/(0.5138*Size_tunnel);
    Gradient_hor_right_top3 = (Ph_right_top3)/(0.3305*Size_tunnel);
    Gradient_hor_right_bot1 = Ph_right_bot1/(2.5234*Size_tunnel);
    Gradient_hor_right_bot2 = (Ph_right_bot2 - Ph_right_bot3)/(0.5138*Size_tunnel);
    Gradient_hor_right_bot3 = (Ph_right_bot3)/(0.3305*Size_tunnel);
    for i=A0:A0+15
        ph(i) = (Gradient_hor_right_bot3*(R1*(1-cos(i/180*pi()))));
    end
    for i=A0+16:A1
        ph(i) = (Ph_right_bot3 + Gradient_hor_right_bot2*(R1*(cos(15/180*pi())-cos(i/180*pi()))));
    end
    for i=A1+1:A1+20
        ph(i) = (Ph_right_bot3 + Gradient_hor_right_bot2*(0.3758*Size_tunnel + R2*(cos(22/180*pi())-cos(i/180*pi()))));
    end
    for i=A1+21:A2
        ph(i) = (Ph_right_bot1 - Gradient_hor_right_bot1*(R2*(cos(43/180*pi())-cos(i/180*pi()))));
    end
    for i=A2+1:A3
        ph(i) = (Ph_right_bot1 - Gradient_hor_right_bot1*(0.2080*Size_tunnel + R3*(cos(63/180*pi())-cos(i/180*pi()))));
    end
    for i=A3+1:A4
        ph(i) = Gradient_hor_right_top1*(R3*(cos(90/180*pi())-cos(i/180*pi())));
    end
    for i=A4+1:A4+20
        ph(i) = Gradient_hor_right_top1*(2.3154*Size_tunnel + R2*(cos(110/180*pi())-cos(i/180*pi())));
    end
    for i=A4+21:A5
        ph(i) = (Ph_right_top2 - Gradient_hor_right_top2*(R2*(cos(138/180*pi())-cos(i/180*pi()))));
    end
    for i=A5+1:A5+7
        ph(i) = (Ph_right_top2 - Gradient_hor_right_top2*(0.1380*Size_tunnel + R1*(cos(158/180*pi())-cos(i/180*pi()))));
    end
    for i=A5+8:A6
        ph(i) = (Ph_right_top3 - Gradient_hor_right_top3*(R1*(cos(165/180*pi())-cos(i/180*pi()))));
    end
    % ---------------------------------------
    % horizontal Loads at LEFT side of tunnel  
    Ph_left_top1 = Ph_rb_1; 
    Ph_left_top2 = -Ph_rb_2;
    Ph_left_top3 = Ph_rb_3; 
    Ph_left_bot1 = -Ph_rt_1; 
    Ph_left_bot2 = Ph_rt_2;
    Ph_left_bot3 = -Ph_rt_3; 
    Gradient_hor_left_top1 = Ph_left_top1/(2.5234*Size_tunnel);
    Gradient_hor_left_top2 = (Ph_left_top2 - Ph_left_top3)/(0.5138*Size_tunnel);
    Gradient_hor_left_top3 = (Ph_left_top3)/(0.3305*Size_tunnel);
    Gradient_hor_left_bot1 = Ph_left_bot1/(2.5234*Size_tunnel);
    Gradient_hor_left_bot2 = (Ph_left_bot2 - Ph_left_bot3)/(0.5138*Size_tunnel);
    Gradient_hor_left_bot3 = (Ph_left_bot3)/(0.3305*Size_tunnel);
    for i=A6+1:A6+15
        ph(i) = (Gradient_hor_left_top3*(R1*(cos(i/180*pi()) - cos(pi()))));
    end
    for i=A6+16:A7
        ph(i) = (Ph_left_top3 + Gradient_hor_left_top2*(R1*(cos(i/180*pi()) - cos(195/180*pi()))));
    end
    for i=A7+1:A7+20
        ph(i) = (Ph_left_top3 + Gradient_hor_left_top2*(0.3758*Size_tunnel + R2*(cos(i/180*pi()) - cos(202/180*pi()))));
    end
    for i=A7+21:A8
        ph(i) = (Ph_left_top1 - Gradient_hor_left_top1*(R2*(cos(i/180*pi()) - cos(223/180*pi()))));
    end
    for i=A8+1:A9
        ph(i) = (Ph_left_top1 - Gradient_hor_left_top1*(0.2080*Size_tunnel + R3*(cos(i/180*pi()) - cos(243/180*pi()))));
    end
    for i=A9+1:A10
        ph(i) = Gradient_hor_left_bot1*(R3*(cos(i/180*pi())-cos(270/180*pi())));
    end
    for i=A10+1:A10+20
        ph(i) = Gradient_hor_left_bot1*(2.3154*Size_tunnel + R2*(cos(i/180*pi())-cos(297/180*pi())));
    end    
    for i=A10+21:A11
        ph(i) = (Ph_left_bot2 - Gradient_hor_left_bot2*(R2*(cos(i/180*pi())-cos(318/180*pi()))));
    end
    for i=A11+1:A11+7
        ph(i) = (Ph_left_bot2 - Gradient_hor_left_bot2*(0.1380*Size_tunnel + R1*(cos(i/180*pi())-cos(338/180*pi()))));
    end
    for i=A11+8:A12
        ph(i) = (Ph_left_bot3 - Gradient_hor_left_bot3*(R1*(cos(i/180*pi())-cos(345/180*pi()))));
    end  
%--------------------------------------------------------------------------
for i=1:360
% estimation of initial stiffness of normal springs
KNin(i)= (0.1002*log(Eground(i))+1.1136+(0.0812*Tunnel_width - 0.7874))*Eground(i)/(1+muy(i))/Rex(i); 
% estimation of initial stiffness of shear springs
KSin(i)=0.00000001; %1/3*KNin(i)
%--------------------------------------------------------------------------
% limit of normal pressure of the ground
pclim(i)=2*c(i)*cos(fi(i))/(1-sin(fi(i)))+(1+sin(fi(i)))/(1-sin(fi(i)))*(pv(i)-ph(i))/2*(muy(i)/(1-muy(i)));
end
% sign used for the application of loads
sign_vertical(1:360,1)=sign_vertical_1_360;
sign_vertical(90:270,1)=sign_vertical_90_270;
sign_horizontal(1:360,1)=sign_horizontal_1_360;
%--------------------------------------------------------------------------
for i=1:360
% parameter for the stress-strain response of the ground (asymptotic curve)
bN(i)=1/pclim(i);
aN(i)=KNin(i)/pclim(i)^2;
% limit of shear stress of the ground-lining interface
taulim(i)=abs((ph(i))+pv(i))/2*tan(fi(i));
% parameter for the stress-strain response of the ground (asymptotic curve)
bS(i)=1/taulim(i);
aS(i)=KSin(i)/taulim(i)^2;
% estimation of initial stiffness of rotational springs (according to
% Janssen's method)
KRin(i)=B*thick^2*E/12;
end
%--------------------------------------------------------------------------
%% global stiffness matrix (initially all zero)
Kctot=zeros(3*360,3*360);
% z4 vector is a vector of normal displacements for each node
z4=zeros(360,1);
% z5 vector is a vector of shear displacements for each node
z5=zeros(360,1);
% KN vector of normal spring
KN=zeros(360,1);
% KS vector of tangential spring
KS=zeros(360,1);
% p vector of reaction force
p=zeros(360,1);
% z6 vector of ratio of M/(N*a) at each node
z6=zeros(360,1);
% Rj vector of fixity factor at each node
Rj=zeros(360,1);
% C correction matrix to take into account the effect of joints
C=zeros(6,6);
% Kc_ele stiffness matrix of element in local system
Kc_ele=zeros(6,6);
% T transformation matrix between the local system and the global system
T=zeros(6,6);
% Kel matrix of element in the global system
Kel=zeros(6,6);
for i=1:num_el
    Jnode(i)=J; % Assumed inertia moment of lining section at each node
    Jele(i)=J; % Assumed average inertia moment of lining section at each element
    Anode(i)=A; % Assumed inertia moment of lining section at each node
    Aele(i)=A; % Assumed average inertia moment of lining section at each element
    Enode(i)=E; % Assumed Young modulus of lining section at each node
    Eele(i)=E; % Assumed average Young modulus of lining section at each element
end
%--------------------------------------------------------------------------
for i=1:num_el
    bc(i)=B; 
    bj(i)=0;  
    Sj(i)=thick; 
end
% Assuming a critical value of the different value in maximum bending 
% moment between 2 sucessively iterative process 
M_critical = 0.00001; % MN.m/m
M_diference = 1; % MN.m/m (arbitrarily assumed)
M_previous = 1; % MN.m/m (arbitrarily assumed) 
%--------------------------------------------------------------------------
%=================iterative procedcure to reach convergence================
while M_diference > M_critical
% iterative procedcure to reach convergence
for j=1:num_steps
    % evaluation of spring stiffness for each node
    for i=1:num_el     
        L=Le(i);       
        if z4(i)>0
            KN(i)=KNin(i)*(1-j/num_steps);
            p(i)=0;
        elseif z4(i)==0
            KN(i)=KNin(i)*j/num_steps;
            p(i)=0;
        else
            p(i)=1*(pclim(i)-1/(aN(i)*(-z4(i))+bN(i)));
            KN(i)=1*p(i)/(-z4(i));
        end
        if z5(i)==0
            KS(i)=KSin(i)*j/num_steps;
        else
            tau=taulim(i)-1/(aS(i)*abs(z5(i))+bS(i));
            KS(i)=tau/abs(z5(i))*j/num_steps;
        end
        if abs(z6(i))<=1/6
            KR(i)=KRin(i)*j/num_steps;
        else
            alfa=9*thick*B*E/8;
            KR(i)=(alfa*abs(M(i))*(2*z6(i)-1)^2/abs(Ax(i)))*j/num_steps;
        end
        KN(i) = KNin(i)/3*j/num_steps;
        KS(i) = KSin(i)/3*j/num_steps;
        Rj(i) = 1;
        if i<=(num_el-1);
            C(1,1)=1;
            C(2,2)=(4*Rj(i+1)-2*Rj(i)+Rj(i)*Rj(i+1))/(4-Rj(i)*Rj(i+1));
            C(2,3)=-2*L*Rj(i)*(1-Rj(i+1))/(4-Rj(i)*Rj(i+1));
            C(3,2)=6*(Rj(i)-Rj(i+1))/(L*(4-Rj(i)*Rj(i+1)));
            C(3,3)=3*Rj(i)*(2-Rj(i+1))/(4-Rj(i)*Rj(i+1));
            C(4,4)=1;
            C(5,5)=(4*Rj(i)-2*Rj(i+1)+Rj(i)*Rj(i+1))/(4-Rj(i)*Rj(i+1));
            C(5,6)=2*L*Rj(i+1)*(1-Rj(i))/(4-Rj(i)*Rj(i+1));
            C(6,5)=6*(Rj(i)-Rj(i+1))/(L*(4-Rj(i)*Rj(i+1)));
            C(6,6)=3*Rj(i+1)*(2-Rj(i))/(4-Rj(i)*Rj(i+1));
        elseif i==num_el
            C(1,1)=1;
            C(2,2)=(4*Rj(1)-2*Rj(num_el)+Rj(num_el)*Rj(1))/(4-Rj(num_el)*Rj(1));
            C(2,3)=-2*L*Rj(num_el)*(1-Rj(1))/(4-Rj(num_el)*Rj(1));
            C(3,2)=6*(Rj(num_el)-Rj(1))/(L*(4-Rj(num_el)*Rj(1)));
            C(3,3)=3*Rj(num_el)*(2-Rj(1))/(4-Rj(num_el)*Rj(1));
            C(4,4)=1;
            C(5,5)=(4*Rj(num_el)-2*Rj(1)+Rj(num_el)*Rj(1))/(4-Rj(num_el)*Rj(1));
            C(5,6)=2*L*Rj(1)*(1-Rj(num_el))/(4-Rj(num_el)*Rj(1));
            C(6,5)=6*(Rj(num_el)-Rj(1))/(L*(4-Rj(num_el)*Rj(1)));
            C(6,6)=3*Rj(1)*(2-Rj(num_el))/(4-Rj(num_el)*Rj(1));
        end
        %------------------------------------------------------------------
        if i==1 
            Jele(i) = 1/2*(Jnode(i)+Jnode(i+1));
            Aele(i) = 1/2*(Anode(i)+Anode(i+1));
            Eele(i) = 1/2*(Enode(i)+Enode(i+1));
            alfa=beta(1);
            T(1,1)=cos(alfa);
            T(2,1)=-sin(alfa);
            T(1,2)=sin(alfa);
            T(2,2)=cos(alfa);
            T(3,3)=1;
            T(4,4)=cos(alfa);
            T(5,4)=-sin(alfa);
            T(4,5)=sin(alfa);
            T(5,5)=cos(alfa);
            T(6,6)=1;
            
            Kc_ele(1,1)=Eele(i)*Aele(i)/L;
            Kc_ele(4,1)=-Eele(i)*Aele(i)/L;
            Kc_ele(2,2)=12*Eele(i)*Jele(i)/L^3;
            Kc_ele(3,2)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(5,2)=-12*Eele(i)*Jele(i)/L^3;
            Kc_ele(6,2)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(2,3)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(3,3)=4*Eele(i)*Jele(i)/L;
            Kc_ele(5,3)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(6,3)=2*Eele(i)*Jele(i)/L;
            Kc_ele(1,4)=-Eele(i)*Aele(i)/L;
            Kc_ele(4,4)=Eele(i)*Aele(i)/L;
            Kc_ele(2,5)=-12*Eele(i)*Jele(i)/L^3;
            Kc_ele(3,5)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(5,5)=12*Eele(i)*Jele(i)/L^3;
            Kc_ele(6,5)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(2,6)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(3,6)=2*Eele(i)*Jele(i)/L;
            Kc_ele(5,6)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(6,6)=4*Eele(i)*Jele(i)/L;
            
            Kc_temp=Kc_ele*C;
            
            Kc_temp(1,1)=Kc_temp(1,1);
            Kc_temp(2,2)=Kc_temp(2,2);
            Kc_temp(4,4)=Kc_temp(4,4);
            Kc_temp(5,5)=Kc_temp(5,5);
            Kc=T'*Kc_temp*T;
            Kctot(1:6,1:6)=Kc;
        elseif (i>=2)&(i<=(num_el-1))
            M3=Kc;
            alfa=beta(i);
            Jele(i) = 1/2*(Jnode(i)+Jnode(i+1));
            Aele(i) = 1/2*(Anode(i)+Anode(i+1));
            Eele(i) = 1/2*(Enode(i)+Enode(i+1)); 
            T(1,1)=cos(alfa);
            T(2,1)=-sin(alfa);
            T(1,2)=sin(alfa);
            T(2,2)=cos(alfa);
            T(3,3)=1;
            T(4,4)=cos(alfa);
            T(5,4)=-sin(alfa);
            T(4,5)=sin(alfa);
            T(5,5)=cos(alfa);
            T(6,6)=1;
            Kc_ele(1,1)=Eele(i)*Aele(i)/L;
            Kc_ele(4,1)=-Eele(i)*Aele(i)/L;
            Kc_ele(2,2)=12*Eele(i)*Jele(i)/L^3;
            Kc_ele(3,2)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(5,2)=-12*Eele(i)*Jele(i)/L^3;
            Kc_ele(6,2)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(2,3)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(3,3)=4*Eele(i)*Jele(i)/L;
            Kc_ele(5,3)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(6,3)=2*Eele(i)*Jele(i)/L;
            Kc_ele(1,4)=-Eele(i)*Aele(i)/L;
            Kc_ele(4,4)=Eele(i)*Aele(i)/L;
            Kc_ele(2,5)=-12*Eele(i)*Jele(i)/L^3;
            Kc_ele(3,5)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(5,5)=12*Eele(i)*Jele(i)/L^3;
            Kc_ele(6,5)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(2,6)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(3,6)=2*Eele(i)*Jele(i)/L;
            Kc_ele(5,6)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(6,6)=4*Eele(i)*Jele(i)/L;
            Kc_temp=Kc_ele*C;
            Kc_temp(1,1)=Kc_temp(1,1);
            Kc_temp(2,2)=Kc_temp(2,2);
            Kc_temp(4,4)=Kc_temp(4,4);
            Kc_temp(5,5)=Kc_temp(5,5);
            Kc=T'*Kc_temp*T;     
            Kctot((3*i-2):(3*i+3),(3*i-2):(3*i+3))=Kc;
            Kctot((3*i-2):(3*i),(3*i-2):(3*i))=Kctot((3*i-2):(3*i),(3*i-2):(3*i))+M3(4:6,4:6);
        elseif i==num_el
            alfa=beta(num_el);
            Jele(num_el) = 1/2*(Jnode(num_el)+Jnode(1));
            Aele(num_el) = 1/2*(Anode(num_el)+Anode(1));
            Eele(num_el) = 1/2*(Enode(num_el)+Enode(1));
            T(1,1)=cos(alfa);
            T(2,1)=-sin(alfa);
            T(1,2)=sin(alfa);
            T(2,2)=cos(alfa);
            T(3,3)=1;
            T(4,4)=cos(alfa);
            T(5,4)=-sin(alfa);
            T(4,5)=sin(alfa);
            T(5,5)=cos(alfa);
            T(6,6)=1;
            Kc_ele(1,1)=Eele(i)*Aele(i)/L;
            Kc_ele(4,1)=-Eele(i)*Aele(i)/L;
            Kc_ele(2,2)=12*Eele(i)*Jele(i)/L^3;
            Kc_ele(3,2)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(5,2)=-12*Eele(i)*Jele(i)/L^3;
            Kc_ele(6,2)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(2,3)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(3,3)=4*Eele(i)*Jele(i)/L;
            Kc_ele(5,3)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(6,3)=2*Eele(i)*Jele(i)/L;
            Kc_ele(1,4)=-Eele(i)*Aele(i)/L;
            Kc_ele(4,4)=Eele(i)*Aele(i)/L;
            Kc_ele(2,5)=-12*Eele(i)*Jele(i)/L^3;
            Kc_ele(3,5)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(5,5)=12*Eele(i)*Jele(i)/L^3;
            Kc_ele(6,5)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(2,6)=6*Eele(i)*Jele(i)/L^2;
            Kc_ele(3,6)=2*Eele(i)*Jele(i)/L;
            Kc_ele(5,6)=-6*Eele(i)*Jele(i)/L^2;
            Kc_ele(6,6)=4*Eele(i)*Jele(i)/L;
            Kc_temp=Kc_ele*C;
            Kc_temp(1,1)=Kc_temp(1,1);
            Kc_temp(2,2)=Kc_temp(2,2);
            Kc_temp(4,4)=Kc_temp(4,4);
            Kc_temp(5,5)=Kc_temp(5,5);
            Kc=T'*Kc_temp*T;     
            Kctot(3*num_el-2:3*num_el,3*num_el-2:3*num_el)=Kc(1:3,1:3)+Kctot(3*num_el-2:3*num_el,3*num_el-2:3*num_el);
            Kctot(3*num_el-2:3*num_el,1:3)=Kc(1:3,4:6);
            Kctot(1:3,3*num_el-2:3*num_el)=Kc(4:6,1:3);
            Kctot(1:3,1:3)=Kc(4:6,4:6)+Kctot(1:3,1:3);
        end
      % modification of the global matrix in order to consider the presence of normal and shear springs
      Kctot(3*i-2,3*i-2)=Kctot(3*i-2,3*i-2)+cos(eps(i))*KN(i)*L*B*cos(delta(i)/2)*cos(eps(i))+sin(eps(i))*KS(i)*L*B*cos(delta(i)/2)*sin(eps(i));
      Kctot(3*i-1,3*i-1)=Kctot(3*i-1,3*i-1)+sin(eps(i))*KN(i)*L*B*cos(delta(i)/2)*sin(eps(i))+cos(eps(i))*KS(i)*L*B*cos(delta(i)/2)*cos(eps(i));
      Kctot(3*i-2,3*i-1)=Kctot(3*i-2,3*i-1)+sin(eps(i))*KN(i)*L*B*cos(delta(i)/2)*cos(eps(i))-cos(eps(i))*KS(i)*L*B*cos(delta(i)/2)*sin(eps(i));
      Kctot(3*i-1,3*i-2)=Kctot(3*i-1,3*i-2)+cos(eps(i))*KN(i)*L*B*cos(delta(i)/2)*sin(eps(i))-sin(eps(i))*KS(i)*L*B*cos(delta(i)/2)*cos(eps(i));
    end
    %----------------------------------------------------------------------
    ext_force=zeros(3*num_el,1);
    % begin of determining the external forces
    for i=1:num_el      
        L=Le(i);
        if i<=(num_el-1)
            Fy_eq(i) = pv(i)*sign_vertical(i)*L*cos(beta(i))*B/2;
            ext_force(3*i-1)=ext_force(3*i-1)+Fy_eq(i);
            ext_force(3*i+2)=ext_force(3*i+2)+Fy_eq(i);
        elseif i==num_el
            Fy_eq(num_el) = pv(i)*sign_vertical(num_el)*L*cos(beta(num_el))*B/2;
            ext_force(3*num_el-1)=ext_force(3*num_el-1)+Fy_eq(num_el);
            ext_force(2)=ext_force(2)+Fy_eq(num_el);
        end
        %------------------------------------------------------------------
        if i<=(num_el-1)
            Fx_eq(i) = ph(i)*sign_horizontal(i)*L*sin(beta(i))*B/2;
            ext_force(3*i-2)=ext_force(3*i-2)+Fx_eq(i);
            ext_force(3*i+1)=ext_force(3*i+1)+Fx_eq(i);
        elseif i==num_el
            Fx_eq(num_el)=ph(i)*sign_horizontal(num_el)*L*sin(beta(num_el))*B/2;
            ext_force(3*num_el-2)=ext_force(3*num_el-2)+Fx_eq(num_el);
            ext_force(1)=ext_force(1)+Fx_eq(num_el);
        end
    end
    %end of determining the external forces
    %----------------------------------------------------------------------
    % determination of global displacement vector
    q=Kctot\ext_force;
    Jtotal = 0;
    for i=1:(num_el-1)        
        L=Le(i);
        displ=q((3*(i-1)+1):(3*(i-1)+6),1);
        alfa=beta(i);
        if i<=(num_el-1);    
            C(1,1)=1;
            C(2,2)=(4*Rj(i+1)-2*Rj(i)+Rj(i)*Rj(i+1))/(4-Rj(i)*Rj(i+1));
            C(2,3)=-2*L*Rj(i)*(1-Rj(i+1))/(4-Rj(i)*Rj(i+1));
            C(3,2)=6*(Rj(i)-Rj(i+1))/(L*(4-Rj(i)*Rj(i+1)));
            C(3,3)=3*Rj(i)*(2-Rj(i+1))/(4-Rj(i)*Rj(i+1));
            C(4,4)=1;
            C(5,5)=(4*Rj(i)-2*Rj(i+1)+Rj(i)*Rj(i+1))/(4-Rj(i)*Rj(i+1));
            C(5,6)=2*L*Rj(i+1)*(1-Rj(i))/(4-Rj(i)*Rj(i+1));
            C(6,5)=6*(Rj(i)-Rj(i+1))/(L*(4-Rj(i)*Rj(i+1)));
            C(6,6)=3*Rj(i+1)*(2-Rj(i))/(4-Rj(i)*Rj(i+1));
        elseif i==num_el
            C(1,1)=1;
            C(2,2)=(4*Rj(1)-2*Rj(num_el)+Rj(num_el)*Rj(1))/(4-Rj(num_el)*Rj(1));
            C(2,3)=-2*L*Rj(num_el)*(1-Rj(1))/(4-Rj(num_el)*Rj(1));
            C(3,2)=6*(Rj(num_el)-Rj(1))/(L*(4-Rj(num_el)*Rj(1)));
            C(3,3)=3*Rj(num_el)*(2-Rj(1))/(4-Rj(num_el)*Rj(1));
            C(4,4)=1;
            C(5,5)=(4*Rj(num_el)-2*Rj(1)+Rj(num_el)*Rj(1))/(4-Rj(num_el)*Rj(1));
            C(5,6)=2*L*Rj(1)*(1-Rj(num_el))/(4-Rj(num_el)*Rj(1));
            C(6,5)=6*(Rj(num_el)-Rj(1))/(L*(4-Rj(num_el)*Rj(1)));
            C(6,6)=3*Rj(1)*(2-Rj(num_el))/(4-Rj(num_el)*Rj(1));
        end
        T(1,1)=cos(alfa);
        T(2,1)=-sin(alfa);
        T(1,2)=sin(alfa);
        T(2,2)=cos(alfa);
        T(3,3)=1;
        T(4,4)=cos(alfa);
        T(5,4)=-sin(alfa);
        T(4,5)=sin(alfa);
        T(5,5)=cos(alfa);
        T(6,6)=1;
        Kc_ele(1,1)=Eele(i)*Aele(i)/L;
        Kc_ele(4,1)=-Eele(i)*Aele(i)/L;
        Kc_ele(2,2)=12*Eele(i)*Jele(i)/L^3;
        Kc_ele(3,2)=6*Eele(i)*Jele(i)/L^2;
        Kc_ele(5,2)=-12*Eele(i)*Jele(i)/L^3;
        Kc_ele(6,2)=6*Eele(i)*Jele(i)/L^2;
        Kc_ele(2,3)=6*Eele(i)*Jele(i)/L^2;
        Kc_ele(3,3)=4*Eele(i)*Jele(i)/L;
        Kc_ele(5,3)=-6*Eele(i)*Jele(i)/L^2;
        Kc_ele(6,3)=2*Eele(i)*Jele(i)/L;
        Kc_ele(1,4)=-Eele(i)*Aele(i)/L;
        Kc_ele(4,4)=Eele(i)*Aele(i)/L;
        Kc_ele(2,5)=-12*Eele(i)*Jele(i)/L^3;
        Kc_ele(3,5)=-6*Eele(i)*Jele(i)/L^2;
        Kc_ele(5,5)=12*Eele(i)*Jele(i)/L^3;
        Kc_ele(6,5)=-6*Eele(i)*Jele(i)/L^2;
        Kc_ele(2,6)=6*Eele(i)*Jele(i)/L^2;
        Kc_ele(3,6)=2*Eele(i)*Jele(i)/L;
        Kc_ele(5,6)=-6*Eele(i)*Jele(i)/L^2;
        Kc_ele(6,6)=4*Eele(i)*Jele(i)/L;
        Kel=Kc_ele*C*T;
        %------------------------------------------------------------------
        % Internal forces at nodes in the global system  
        F_nod=Kel*displ;
        if i==1
            M(i)=F_nod(3,1);
            Ax(i)=F_nod(1,1);
            z6(i)=M(i)/Ax(i)/thick;
        end
        M(i+1)=-F_nod(6,1);
        Ax(i+1)=-F_nod(4,1);
        z6(i)=M(i)/Ax(i)/thick;       
        % determining height of compressive stress at join section (Sj)
        if bj(i)>0
            a_temp=Ax(i)*bj(i)/3;
            b_temp=M(i)*bj(i)-Ax(i)*bj(i)*thick/2;
            c_temp=M(i)*bc(i)*thick;
            d_temp=Ax(i)*bc(i)*(thick^3)/3;
            e_temp=-Ax(i)*bc(i)*(thick^4)/2;
            equation_temp = [a_temp,b_temp,c_temp,d_temp,e_temp];
            x_temp = roots(equation_temp);
            x_temp=sort(x_temp);
            variable_number=numel(x_temp);
            for m=1:variable_number; 
                phan_ao=imag(x_temp(m));
                if phan_ao==0 ; 
                    if x_temp(m)>0
                        if x_temp(m)<thick
                            Sj(i)=x_temp(m); 
                        end
                    end
                end
            end
            Nj(i)= Ax(i)*bj(i)*Sj(i)/(bj(i)*Sj(i)+bc(i)*thick);
            Mj(i)= Nj(i)*(thick/2-Sj(i)/3);
            Sj(i)=thick;
            Jnode(i)=bc(i)*(thick^3)/12+bj(i)*(Sj(i)^3)/12;
            Anode(i)=bc(i)*thick+bj(i)*Sj(i);
        end
        Jtotal= Jtotal + Jnode(i);           
     end     
    for i=1:(num_el-1)
        displ=q((3*(i-1)+1):(3*(i-1)+6),1);
        z4(i)=-(displ(1,1)*cos(eps(i))+displ(2,1)*sin(eps(i))); 
        z5(i)=-displ(1,1)*sin(eps(i))+displ(2,1)*cos(eps(i));
    end
   displ(1:3)=q(3*num_el-2:3*num_el);
   displ(4:6)=q(1:3);
   z4(num_el)=-(displ(1,1)*cos(eps(num_el))+displ(2,1)*sin(eps(num_el)));
   z5(num_el)=-displ(1,1)*sin(eps(num_el))+displ(2,1)*cos(eps(num_el));
end
    M_max = max(M)
    M_diference = abs(M_max - M_previous)
    M_previous = max(M);
    Jmean=Jtotal/num_el;
end
%==========================================================================
lungh=0;
% evaluation of internal forces (axial force, shear force, bending moments) in all nodes
for i=1:(num_el-1)   
   L=Le(i);
   displ=q((3*(i-1)+1):(3*(i-1)+6),1);
   alfa=beta(i);
   if i<=(num_el-1); 
       C(1,1)=1;
       C(2,2)=(4*Rj(i+1)-2*Rj(i)+Rj(i)*Rj(i+1))/(4-Rj(i)*Rj(i+1));
       C(2,3)=-2*L*Rj(i)*(1-Rj(i+1))/(4-Rj(i)*Rj(i+1));
       C(3,2)=6*(Rj(i)-Rj(i+1))/(L*(4-Rj(i)*Rj(i+1)));
       C(3,3)=3*Rj(i)*(2-Rj(i+1))/(4-Rj(i)*Rj(i+1));
       C(4,4)=1;
       C(5,5)=(4*Rj(i)-2*Rj(i+1)+Rj(i)*Rj(i+1))/(4-Rj(i)*Rj(i+1));
       C(5,6)=2*L*Rj(i+1)*(1-Rj(i))/(4-Rj(i)*Rj(i+1));
       C(6,5)=6*(Rj(i)-Rj(i+1))/(L*(4-Rj(i)*Rj(i+1)));
       C(6,6)=3*Rj(i+1)*(2-Rj(i))/(4-Rj(i)*Rj(i+1));
   elseif i==num_el
       C(1,1)=1;
       C(2,2)=(4*Rj(1)-2*Rj(num_el)+Rj(num_el)*Rj(1))/(4-Rj(num_el)*Rj(1));
       C(2,3)=-2*L*Rj(num_el)*(1-Rj(1))/(4-Rj(num_el)*Rj(1));
       C(3,2)=6*(Rj(num_el)-Rj(1))/(L*(4-Rj(num_el)*Rj(1)));
       C(3,3)=3*Rj(num_el)*(2-Rj(1))/(4-Rj(num_el)*Rj(1));
       C(4,4)=1;
       C(5,5)=(4*Rj(num_el)-2*Rj(1)+Rj(num_el)*Rj(1))/(4-Rj(num_el)*Rj(1));
       C(5,6)=2*L*Rj(1)*(1-Rj(num_el))/(4-Rj(num_el)*Rj(1));
       C(6,5)=6*(Rj(num_el)-Rj(1))/(L*(4-Rj(num_el)*Rj(1)));
       C(6,6)=3*Rj(1)*(2-Rj(num_el))/(4-Rj(num_el)*Rj(1));
   end
   % local matrix of stiffness of the element in the global system
   T(1,1)=cos(alfa);
   T(2,1)=-sin(alfa);
   T(1,2)=sin(alfa);
   T(2,2)=cos(alfa);
   T(3,3)=1;
   T(4,4)=cos(alfa);
   T(5,4)=-sin(alfa);
   T(4,5)=sin(alfa);
   T(5,5)=cos(alfa);
   T(6,6)=1;
   Kc_ele(1,1)=Eele(i)*Aele(i)/L;
   Kc_ele(4,1)=-Eele(i)*Aele(i)/L;
   Kc_ele(2,2)=12*Eele(i)*Jele(i)/L^3;
   Kc_ele(3,2)=6*Eele(i)*Jele(i)/L^2;
   Kc_ele(5,2)=-12*Eele(i)*Jele(i)/L^3;
   Kc_ele(6,2)=6*Eele(i)*Jele(i)/L^2;
   Kc_ele(2,3)=6*Eele(i)*Jele(i)/L^2;
   Kc_ele(3,3)=4*Eele(i)*Jele(i)/L;
   Kc_ele(5,3)=-6*Eele(i)*Jele(i)/L^2;
   Kc_ele(6,3)=2*Eele(i)*Jele(i)/L;
   Kc_ele(1,4)=-Eele(i)*Aele(i)/L;
   Kc_ele(4,4)=Eele(i)*Aele(i)/L;
   Kc_ele(2,5)=-12*Eele(i)*Jele(i)/L^3;
   Kc_ele(3,5)=-6*Eele(i)*Jele(i)/L^2;
   Kc_ele(5,5)=12*Eele(i)*Jele(i)/L^3;
   Kc_ele(6,5)=-6*Eele(i)*Jele(i)/L^2;
   Kc_ele(2,6)=6*Eele(i)*Jele(i)/L^2;
   Kc_ele(3,6)=2*Eele(i)*Jele(i)/L;
   Kc_ele(5,6)=-6*Eele(i)*Jele(i)/L^2;
   Kc_ele(6,6)=4*Eele(i)*Jele(i)/L;
   % local matrix of stiffness of the element in the global system
   Kel=Kc_ele*C*T;
   %-----------------------------------------------------------------------
   % Internal forces at nodes in the global system  
   F_nod=Kel*displ;
   if i==1
       M(i)=F_nod(3,1);
       Ax(i)=F_nod(1,1);
       Sh(i)=F_nod(2,1);
       Epsilon_in(i)=(-Ax(i)/Eele(i)/A+M(i)/Eele(i)/Jele(i))*1000000; 
       Epsilon_ex(i)=(-Ax(i)/Eele(i)/A-M(i)/Eele(i)/Jele(i))*1000000; 
   end
   M(i+1)=-F_nod(6,1);      
   Ax(i+1)=-F_nod(4,1);
   Sh(i+1)=F_nod(5,1);
   Epsilon_in(i+1)=(-Ax(i+1)/Eele(i)/A+M(i+1)/Eele(i)/Jele(i))*1000000;
   Epsilon_ex(i+1)=(-Ax(i+1)/Eele(i)/A-M(i+1)/Eele(i)/Jele(i))*1000000;
   w1(i)=lungh;
   lungh=lungh+360/num_el;
   z1(i)=q((3*i),1)*(180/pi());
   z2(i)=q((3*i-1),1);
   z3(i)=q((3*i-2),1);
   % normal displacement of the node
   z4(i)=-(displ(1,1)*cos(eps(i))+displ(2,1)*sin(eps(i)));
   % shear displacement of the node
   z5(i)=-displ(1,1)*sin(eps(i))+displ(2,1)*cos(eps(i));
   % determination of the reaction pressure and shear stress of the ground
   if z4(i)>=0
       preaz(i)=0;
   else
       preaz(i)=pclim(i)-1/(aN(i)*(-z4(i))+bN(i));
   end
   if z5(i)>=0
       taustress(i)=taulim(i)-1/(aS(i)*z5(i)+bS(i));
   else
       taustress(i)=-(taulim(i)-1/(aS(i)*(-z5(i))+bS(i)));
   end      
end
%--------------------------------------------------------------------------
% distance from the origin of the last node
w1(num_el)=lungh;
z1(num_el)=q((3*num_el),1)*(180/pi());
z2(num_el)=q((3*num_el-1),1);
z3(num_el)=q((3*num_el-2),1);
displ(1:3)=q(3*num_el-2:3*num_el);
displ(4:6)=q(1:3);
% normal displacement of the last node
z4(num_el)=-(displ(1,1)*cos(eps(num_el))+displ(2,1)*sin(eps(num_el)));
% shear displacement of the last node
z5(num_el)=-displ(1,1)*sin(eps(num_el))+displ(2,1)*cos(eps(num_el));
% determination of the reaction pressure and shear stress of the ground
if z4(num_el)>=0
    preaz(num_el)=0;
else
    preaz(num_el)=pclim(num_el)-1/(aN(num_el)*(-z4(num_el))+bN(num_el));
end
if z5(num_el)>=0
    taustress(num_el)=taulim(num_el)-1/(aS(num_el)*z5(num_el)+bS(num_el));
else
    taustress(num_el)=-(taulim(num_el)-1/(aS(num_el)*(-z5(num_el))+bS(num_el)));
end 
%--------------------------------------------------------------------------
M(num_el+1)=M(1);
Ax(num_el+1)=Ax(1);
Sh(num_el+1)=Sh(1);
Epsilon_in(num_el+1)=(-Ax(num_el+1)/Eele(i)/A+M(num_el+1)/Eele(i)/Jele(i))*1000000;
Epsilon_ex(num_el+1)=(-Ax(num_el+1)/Eele(i)/A-M(num_el+1)/Eele(i)/Jele(i))*1000000;
w1(num_el+1)=w1(num_el)+360/num_el;
z1(num_el+1)=z1(1);
z2(num_el+1)=z2(1);
z3(num_el+1)=z3(1);
z4(num_el+1)=z4(1);
z5(num_el+1)=z5(1);
preaz(num_el+1)=preaz(1);
taustress(num_el+1)=taustress(1);
%=======================Export the results=================================
M=M/B;
Ax=Ax/B;
Sh=Sh/B;
%======================Total Static Forces=================================
M_s = M;
Ax_s = Ax;
Sh_s = Sh;
Epsilon_in_s = Epsilon_in;
Epsilon_ex_s = Epsilon_ex;
w1_s = w1;
z1_s = z1;
z2_s = z2;
z3_s = z3;
z4_s = z4;
z5_s = z5;
preaz_s = preaz;
taustress_s = taustress;
maximum_Bending_moment_MNm = max(M);
minimum_Bending_moment_MNm = min(M);
maximum_NormalForce_MN = max(Ax);
minimum_NormalForce_MN = min(Ax);
maximum_ShearForce_MN = max(Sh);
minimum_ShearForce_MN = min(Sh);
maximum_Radial_displacement_m = max(z4);
%---------------------- Display--------------------------------------------
disp(['maximum_Bending_moment_MNm = ',num2str(max(M))]);
disp(' ');
disp(['minimum_Bending_moment_MNm = ', num2str(min(M))]);
disp(' ');
disp(['maximum_NormalForce_MN = ', num2str(max(Ax))]);
disp(' ');
disp(['minimum_NormalForce_MN = ', num2str(min(Ax))]);
disp(' ');
disp(['maximum_ShearForce_MN = ', num2str(max(Sh))]);
disp(' ');
disp(['minimum_ShearForce_MN = ', num2str(min(Sh))]);
disp(' ');
disp(['maximum_Radial_displacement_m = ', num2str(max(z4))]);
%-------------------------------------------------------------------------- 
 a1=M.';
 a2=Ax.';
 a3=Sh.';
 a4=preaz.';
 a5=taustress.';
 a6=z1.';
 a7=z4;
 a8=-z5;
 a9=Epsilon_in.';
 a10=Epsilon_ex.';
 a0=w1';
%==========================================================================
 % Plotting
subplot(2,2,1)
   plot(w1,M), ylabel('Incremental Bending Moment M(MN.m)'),...
   xlabel('Angle (degrees) measured counter-clockwise from tunnel bottom'),... 
   grid on, axis([0 max(w1) min(M) max(M)])
   set(gca,'XTick',0:60:360)
subplot(2,2,2)
   plot(w1,Ax), ylabel('Incremental Normal Forces N(MN)'),...
   xlabel('Angle (degrees) measured counter-clockwise from tunnel bottom'),...
   grid on, axis([0 max(w1) min(Ax) max(Ax)])
   set(gca,'XTick',0:60:360)
%----------------------------------*****-----------------------------------
