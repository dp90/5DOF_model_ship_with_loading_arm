function Max_vector = determine_vessel_motion(zeta,period)

test1 = 0; % Hydrostatic forces
test2 = 0; % displacements
test3 = 0; % visual

Tend  = 100;% Seconds final time
Tstep = 0.05; % Timestep

H.zeta = zeta;
H.period = period;

%% Hydrostatic Parameters
H.rho = 1025;       % [kg/m^3]  Water density
H.g = 9.81;         % [kg*m/s^2]Gravitational constant
H.zeta = 3;         % [m]       Wave amplitude
H.period = 16;      % [s]       Wave period
H.omega = 1/H.period; % [rad]   Wave frequency
H.lambda = 200;     % [m]       Wavelength
H.k = ...           % [-]       Wave number
    2*pi/H.lambda;
H.v = H.omega/H.k;  % [m/s]     Wave velocity


%% Vessel Parameters
Lpp = 133;          % [m]       Length
Bpp = 50;           % [m]       Width
Hpp = 12;           % [m]       Depth
Mpp = 35E6;         % [kg]      Apparent mass
Tpp = ...           % [m]       Draft
    8.5; %Mpp/(H.rho*Bpp*Lpp);
GM  = ...           % [m]       Metacentric height of vessel
    0.5*(Tpp-Hpp)+(Bpp^2)/(12*Tpp); %given in assignment
a   = 10;           % [m]       Distance CoG platform to boom
theta1_0 = ...      % [rad]     Initial angle
    atan((Hpp-Tpp)/a);
Jpp = ...           % [kg*m^2]  Mass  moment of inertia
    Mpp/12*(Hpp^2 + Bpp^2);

H.Bpp=Bpp; H.Tpp=Tpp; H.Lpp=Lpp; H.Mpp=Mpp; H.Hpp=Hpp; % store in H-structure
P.Mpp=Mpp; P.Jpp=Jpp; % store in P-structure

%% Crane Parameters
H_0 = 30;           % [m]       Initial height of the crane
Lb  = 40;           % [m]       Length of the boom
Mb  = 5e6;          % [kg]      Mass of the boom
Jb  = ...           % [kg*m^2]  Mass  moment of inertia of the boom
    Mb*Lb^2/3;
d   = ...           % [m]       Distance CoG boom to CoG platform
    sqrt((a + sqrt(Lb^2-H_0^2)/2)^2 + (Hpp-Tpp+H_0/2)^2);
J3  = ...           % [kg*m^2]  Mass  moment of inertia of boom around platform
    1/12*Mb*Lb^2+Mb*d^2;
theta2_0 = ...      % [rad]     Initial angle
    acos(H_0/Lb);
B_0 = ...           % [m]       Initial horizontal distance of the crane tip
    Lb*sin(theta2_0);

H.Mb=Mb; H.B_0=B_0; H.a=a; H.H_0=H_0; % store in H-structure
P.Jb=Jb; P.J3=J3;  P.Mb=Mb; % store in P-structure

%% Cable Parameters
Lc  = 25;           % [m]       Length of cable supporting the load
Ec  = 210e9;        % [Pa]      Young's modulus
Dc  = 0.15;         % [m]       Cable diameter
Ml  = 400e3;        % [kg]      Mass of the load
theta3_0 = 0;       % [rad]     Initial angle
b   = theta3_0*Lc;  % [m]       Initial horizontal displacement of the load

H.Ml=Ml; H.b=b; H.Lc=Lc; % store in H-structure
P.Ml=Ml; % store in P-structure

%% Mooring Line Parameters
Lmh1= 300;          % [m]       Horizontal length of the left mooring line
Lmh2= 300;          % [m]       Horizontal length of the right mooring line
Lmv = 300;           % [m]       Vertical length of the mooring line

%% Restoring Forces

k1  = ...           % [N/m]     Equivalent stiffness in horizontal direction
    1/4*H.g*(2*Bpp*Lpp*Tpp*H.rho-(Mpp+Mb+Ml))*(Lmh1+2-Lmh2)/Lmv; %self-defined
k4  = k1;           % [N/m]     Equivalent stiffness in horizontal direction
kr  = ...           % [N/theta] Equivalent stiffness in roll rotation
    (Mpp+Mb+Ml)*H.g*GM; %given in assignment
k2  = ...           % [N/m]     Equivalent stiffness in vertical direction
    Bpp*Lpp*H.rho*H.g/2; %given and self-defined
k3  = k2;           % [N/m]     Equivalent stiffness in vertical direction
k5  = ...           % [N/m]     Equivalent stiffness of the cable
    Ec/4*pi*Dc^2/(a+B_0+sqrt(H_0^2+a^2));        
c2  = 1e7;          % [N/m/s]   Equivalent damping coefficient in vertical direction
c3  = c2;           % [N/m/s]   Equivalent damping coefficient in vertical direction

H.Fmasses = ...
    Tpp*(k2+k3);

P.k1=k1; P.k2=k2; P.k3=k3; P.k4=k4; P.k5=k5; P.c2=c2; P.c3=c3; % store in P-structure

%% Pressure test

if test1==1

for jj = 1:Tend %seconds
    
    t_n = jj;
    xlocjj=jj;

    z=Tpp;
    
for ii = 1:Bpp+1
    x=ii-16;
    p = compute_pressure(x, z, t_n, H);
    pressure_0(ii) = p;
    xlocii(ii)=x;
end

[Fvec] = compute_loads(t_n, H);

Fh=Fvec(1);
Fv=Fvec(2);
Mtot=Fvec(3);


figure(1)
subplot(4,1,1)
scatter(t_n,Mtot)
hold on
xlim([1 100])
title('Moment')

subplot(4,1,2)
scatter(t_n,Fh)
xlim([1 100])
hold on
title('Fh')

subplot(4,1,3)
scatter(t_n,Fv)
xlim([1 100])
hold on
title('Fv')

subplot(4,1,4)
plot(xlocii,pressure_0)
% ylim([1e3 7e4])
title('Pressure')

pause(0.1)

end

end

%% ODE running

q0 = zeros(1,10).';
tvec=[0:Tstep:Tend];

[T,Q] = ode45(@(t_n, q_n) solve_statespace_vector(t_n, q_n, P, H), tvec, q0);

%% Collecting maximums

max1=max(Q(:,1));
max2=max(Q(:,1));
max3=max(Q(:,3))*360/(2*pi);
max4=max(Q(:,4))*360/(2*pi);
max5=max(Q(:,5))*360/(2*pi);

Max_vector = [max1; max2; max3; max4; max5];



%% Check loads

if test2==1    

figure(2)
plot(T,Q(:,1),'Linestyle','-.','Linewidth',2)
hold on
plot(T,Q(:,2),'Linestyle','--','Linewidth',2)
hold on
plot(T,Q(:,3)*360/(2*pi),'Linestyle','-')
hold on 
plot(T,Q(:,4)*360/(2*pi),'Linestyle','-')
hold on
plot(T,Q(:,5)*360/(2*pi),'Linestyle','-')
ylabel('Rotation [Degrees] / Displacement [m]')
xlabel('Time [s]')
legend('X1','Y1','\theta1','\theta2','\theta3')
end

%% Check visual

if test3==1     
    for jj=1:numel(T);
        
      
    theta1=Q(jj,3);
    theta2=Q(jj,4);
    theta3=Q(jj,5);
    x1=Q(jj,1);
    y1=Q(jj,2);
    
    phi1 = atan((Hpp-Tpp)/(Bpp/2));
    phi2 = atan(Tpp/(Bpp/2));
    phi3 = atan((Hpp-Tpp)/(a));
    phi4 = atan(H_0/B_0);
    
    y5 = sin(phi1-theta1)*sqrt((Hpp-Tpp)^2+((1/2)*Bpp)^2)+y1;
    x5 = cos(phi1-theta1)*sqrt((Hpp-Tpp)^2+((1/2)*Bpp)^2)+x1;
    y6 = -sin(phi2+theta1)*sqrt(Tpp^2+((1/2)*Bpp)^2)+y1;
    x6 = cos(phi2+theta1)*sqrt(Tpp^2+((1/2)*Bpp)^2)+x1;
    y7 = -sin(phi2-theta1)*sqrt(Tpp^2+((1/2)*Bpp)^2)+y1;
    x7 = -cos(phi2-theta1)*sqrt(Tpp^2+((1/2)*Bpp)^2)+x1;
    y8 = sin(phi1+theta1)*sqrt((Hpp-Tpp)^2+((1/2)*Bpp)^2)+y1;
    x8 = -cos(phi1+theta1)*sqrt((Hpp-Tpp)^2+((1/2)*Bpp)^2)+x1;
    y9 = sin(phi3+theta1)*sqrt((Hpp-Tpp)^2+a^2)+y1;
    x9 = -cos(phi3+theta1)*sqrt((Hpp-Tpp)^2+a^2)+x1;
    y2 = sin(phi3-theta1)*sqrt((Hpp-Tpp)^2+a^2)+y1;
    x2 = cos(phi3-theta1)*sqrt((Hpp-Tpp)^2+a^2)+x1;
    x10 = sin(theta1)*(Hpp-Tpp+H_0)+x1;
    y10 = cos(theta1)*(Hpp-Tpp+H_0)+y1;
    x3 = cos(phi4-theta2)*sqrt(B_0^2+H_0^2)+x2;
    y3 = sin(phi4-theta2)*sqrt(B_0^2+H_0^2)+y2;
    x4 = -Lc*sin(theta3)+x3;
    y4 = -Lc*cos(theta3)+y3;

    Linevec = ...
    [x5 x8 y5 y8;
     x8 x7 y8 y7;
     x7 x6 y7 y6;
     x6 x5 y6 y5;
     x9 x10 y9 y10;
     x2 x3 y2 y3;
     x3 x10 y3 y10;
     x10 x2 y10 y2;
     x3 x4 y3 y4;
     -Bpp/2-10 Bpp/2+30 0 0;
     30 30 -10 5;
     30 60 5 5;
     60 60 5 -10];
 
     x11=linspace(-Bpp/2-10,Bpp/2+5,200);
     y11=H.zeta*sin(H.k*x11-H.omega*T(jj));
    
    figure(3)
    for ii=1:numel(Linevec(:,1))
        
        if ii == 10
            %line([Linevec(ii,1),Linevec(ii,2)],[Linevec(ii,3) Linevec(ii,4)],'Color','blue','Linewidth',1)
        elseif ii == 9 || ii == 7
            line([Linevec(ii,1),Linevec(ii,2)],[Linevec(ii,3) Linevec(ii,4)],'Color','black','Linewidth',1)
        else
            line([Linevec(ii,1),Linevec(ii,2)],[Linevec(ii,3) Linevec(ii,4)],'Color','red','Linewidth',2)
        end
        
        line(x11,y11,'Color','blue')     
        
        xlim([-Bpp/2-10 Bpp/2+30])
        ylim([-10 50])
        pause(0.01)
    end
    
    if mod(jj,3)==0        
        clf
        plot(x4,y4,'s','Color','black')
        
    end
    
    title('Vessel motion')
    legend(['T=' num2str(T(jj))])
    xlabel('[m]')
    ylabel('[m]')
    
    end    
end

end
