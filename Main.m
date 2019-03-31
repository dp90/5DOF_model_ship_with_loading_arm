clear all; 
clc; 
dbstop if error;
close('all');

zeta_vec = [0:0.4:10];
period_vec = [0:0.8:20];

limit_surf1 = zeros(numel(zeta_vec),numel(period_vec));
limit_surf2 = zeros(numel(zeta_vec),numel(period_vec));
limit_surf3 = zeros(numel(zeta_vec),numel(period_vec));
limit_surf4 = zeros(numel(zeta_vec),numel(period_vec));
limit_surf5 = zeros(numel(zeta_vec),numel(period_vec));

check=0;

totalcheck=numel(zeta_vec)*numel(period_vec);

for zz = 1:numel(zeta_vec);
    
    zeta = zeta_vec(zz);
    
    for pp = 1:numel(period_vec);
        
        check=check+1;
        
        clc
        procent = check/totalcheck*100
        
        period = period_vec(pp);
        
        [Max_vector] = determine_vessel_motion(zeta,period);
        
        limit_surf1(zz,pp) = Max_vector(1);
        limit_surf2(zz,pp) = Max_vector(2);
        limit_surf3(zz,pp) = Max_vector(3);
        limit_surf4(zz,pp) = Max_vector(4);
        limit_surf5(zz,pp) = Max_vector(5);
        
    end
    
end

limiter=ones(numel(zeta_vec),numel(period_vec));

save('limits.mat','limit_surf1','limit_surf2','limit_surf3','limit_surf4','limit_surf5')

figure()
surf(zeta_vec,period_vec,limit_surf3,'FaceColor','r','EdgeColor','none')
xlabel('Wave Amplitude [m]')
ylabel('Wave period [s]')
zlabel('Max Rotation \theta3 [Degrees]')
% zlim([0 45])
hold on
surf(zeta_vec,period_vec,limiter*9,'FaceColor','g','EdgeColor','none');

figure()
surf(zeta_vec,period_vec,limit_surf4,'FaceColor','r','EdgeColor','none')
xlabel('Wave Amplitude [m]')
ylabel('Wave period [s]')
zlabel('Max Rotation \theta2 [Degrees]')
% zlim([0 45])
hold on
surf(zeta_vec,period_vec,limiter*9,'FaceColor','g','EdgeColor','none');

figure()
surf(zeta_vec,period_vec,limit_surf5,'FaceColor','r','EdgeColor','none')
xlabel('Wave Amplitude [m]')
ylabel('Wave period [s]')
zlabel('Max Rotation \theta1 [Degrees]')
% zlim([0 45])
hold on
surf(zeta_vec,period_vec,limiter*9,'FaceColor','g','EdgeColor','none');

figure()
surf(zeta_vec,period_vec,limit_surf1,'FaceColor','r','EdgeColor','none')
xlabel('Wave Amplitude [m]')
ylabel('Wave period [s]')
zlabel('Max Horizontal displacement [m]')
% zlim([0 45])
hold on
surf(zeta_vec,period_vec,limiter*10,'FaceColor','g','EdgeColor','none');

figure()
surf(zeta_vec,period_vec,limit_surf2,'FaceColor','r','EdgeColor','none')
xlabel('Wave Amplitude [m]')
ylabel('Wave period [s]')
zlabel('Max Vertical displacement [m]')
% zlim([0 45])
hold on
surf(zeta_vec,period_vec,limiter*2,'FaceColor','g','EdgeColor','none');

