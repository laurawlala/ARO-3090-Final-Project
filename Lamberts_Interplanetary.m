function [V_inf_dep, V_inf_arr] = Lamberts_Interplanetary(Planet_1,Planet_2,yyyy,mm,dd,TOF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 132712440018; % Sun's gravitational parameter

%% Hohmann Transfer
[R1,~,mu1,T1,n1,V1,r_SOI1] = func_parameters_Planet(Planet_1);
[R2,~,mu2,T2,n2,V2,r_SOI2] = func_parameters_Planet(Planet_2);
T_syn = 2*pi/(abs(n1-n2));

at = (R2 + R1)/2;
et = abs(R2 - R1)/(R2 + R1);
ht = sqrt(mu*(1-et^2)*at);
Tt = 2*pi/sqrt(mu)*at^1.5;
phi_dep = mod(180 - (360/T2)*(Tt/2),360);
if phi_dep>180
    phi_dep = 360-phi_dep;
end

V_dep     = ht/R1;
V_inf_dep = abs(V_dep - V1);
V_arr     = ht/R2;
V_inf_arr = abs(V2 - V_arr);
if R1==R2
    fprintf('Please select two different planets.\n\n')
    return
end
%fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
%fprintf('             Hohmann Transfer from %s to %s \n\n',Planet_1,Planet_2)
%fprintf('     Observer (O): %s, Primary (P): Sun, Target (T): %s\n',Planet_1,Planet_2)
%fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
%fprintf('O-P-T Angle at Departure:      phi = %.1f deg   (0 < phi < 180 deg)\n',phi_dep)    
%fprintf('          Time of Flight:      TOF = %.0f days\n',Tt/86400/2)  
%fprintf('          Synodic Period:    T_Syn = %.0f days\n\n',T_syn/86400)  
%fprintf('V_inf (dep)  = %.3f km/s,  C3 (dep) =  %.1f km^2/s^2\n',V_inf_dep, V_inf_dep^2)
%fprintf('V_inf (arr)  = %.3f km/s,  C3 (arr) =  %.1f km^2/s^2\n\n',V_inf_arr, V_inf_arr^2)

tic
%% Julian Dates for Departure
JD_Dep = juliandate(yyyy,mm,dd);

Options = optimoptions('fsolve','Display','off','TolFun',1E-12,'TolX',1E-12);
% Using 432t as it is the latest Ephemeris available for the function
[r1_vec,VP1_dep] = planetEphemeris(JD_Dep','Sun',Planet_1,'432t');  
 % Creating column vectors
r1_vec=r1_vec';   VP1_dep=VP1_dep';   r1 = vecnorm(r1_vec);
JD_Arr = JD_Dep + TOF; % Julian Day

Arr = datetime(JD_Arr,'convertfrom','juliandate');
yyyyA=year(Arr); mmA=month(Arr); ddA=day(Arr);

% Using 432t as it is the latest Ephemeris available for the function
[r2_vec,VP2_arr] = planetEphemeris(JD_Arr','Sun',Planet_2,'432t');

% Creating column vectors
r2_vec=r2_vec';   VP2_arr=VP2_arr';   r2 = vecnorm(r2_vec);

% True Anomaly Difference
DTA_S = acos(dot(r1_vec,r2_vec)./(r1.*r2));
DTA_L = 2*pi-DTA_S;    

K = r1.*r2.*(1-cos(DTA_S));
L = r1 +r2;
M = r1.*r2.*(1+cos(DTA_S));

N_P=2000; % Has to be >1
  
p_min = K./(L+(2*M).^.5)*1.000001;
p_max = K./(L-(2*M).^.5)*0.999999;

dp = (p_max-p_min)/(N_P-1);
p = repmat(p_min,N_P,1);
p(2:N_P,1) = p_min + dp.*(1:N_P-1).';   

%% SHORT-way Orbit (_S)
a  = M.*K.*p./((2*M-L.^2).*p.^2+2*K.*L.*p-K.^2);
f  = 1-r2./p.*(1-cos(DTA_S));
dot_f = sqrt(mu./p).*((1-cos(DTA_S))./p-1./r1-1./r2).*tan(DTA_S/2);
cosDE = 1-r1./a.*(1-f);
sinDE = -r1.*r2.*dot_f./sqrt(mu.*a);
DE    = mod(atan2(sinDE,cosDE),2*pi);
tof   = ((a.^3./mu).^.5.*(DE-sinDE) + r1.*r2./(p.*mu).^.5.*sin(DTA_S))/86400;
d_tof = abs(tof - TOF);

[~,J_min] = min(d_tof);

p0 = p(J_min,1);

tof_min = tof(1,:);
tof_max = tof(end,:);

if (tof_min<TOF&&tof_max>TOF) || (tof_min>TOF&&tof_max<TOF)
    a0 = M.*K.*p0./((2*M-L.^2).*p0.^2+2*K.*L.*p0-K.^2);
    f0 = 1-r2./p0.*(1-cos(DTA_S));
    dot_f0 = sqrt(mu./p0).*((1-cos(DTA_S))./p0-1./r1-1./r2).*tan(DTA_S/2);
    cosDE0 = 1-r1./a0.*(1-f0);
    sinDE0 = -r1.*r2.*dot_f0./sqrt(mu.*a0);
    DE0 = mod(atan2(sinDE0,cosDE0),2*pi);          
    
    PAR = [r1,r2,DTA_S,TOF*86400,mu,1];
    X0 = [a0,DE0,p0]';
    
    X = fsolve(@(X)func_Lambert(X,PAR),X0,Options);

    a_S  = X(1);
    DE_S = X(2);
    p_S  = X(3);
       
    f    = 1-r2./p_S.*(1-cos(DTA_S));
    g    = r1.*r2./sqrt(p_S.*mu).*sin(DTA_S);
    dot_f = sqrt(mu./p_S).*((1-cos(DTA_S))./p_S-1./r1-1./r2).*tan(DTA_S/2);
    dot_g = 1-r1./p_S.*(1-cos(DTA_S));

    v1_vec_S = (r2_vec - f.*r1_vec)./g;
    v2_vec_S = dot_f.*r1_vec + dot_g./g.*(r2_vec - f.*r1_vec);

    TOF_check = ((a_S.^3./mu).^.5.*(DE_S-sin(DE_S)) + r1.*r2./(p_S.*mu).^.5.*sin(DTA_S))/86400;  
    delta_TOF = max(abs(TOF - TOF_check));

    if delta_TOF>1E-7
        fprintf('SHORT-Way Orbit may not reliable for TOF = %d days (Error: %.1f days)\n\n',TOF,delta_TOF)
    end
    % Departure
    vinf_dep_vec_S = v1_vec_S - VP1_dep;
    vinf_dep_S = vecnorm(vinf_dep_vec_S);
    C3_dep_S = vinf_dep_S^2;

    % Arrival
    vinf_arr_vec_S = VP2_arr - v2_vec_S;
    vinf_arr_S = vecnorm(vinf_arr_vec_S);
    C3_arr_S = vinf_arr_S^2;
else
    C3_dep_S = NaN;
    C3_arr_S = NaN;
end

%% LONG-way Orbit (_L)
dot_f = sqrt(mu./p).*((1-cos(DTA_L))./p-1./r1-1./r2).*tan(DTA_L/2);
sinDE = -r1.*r2.*dot_f./sqrt(mu.*a);
DE    = mod(atan2(sinDE,cosDE),2*pi);
tof   = ((a.^3./mu).^.5.*(DE-sinDE) + r1.*r2./(p.*mu).^.5.*sin(DTA_L))/86400;
d_tof = abs(tof - TOF);

[~,J_min] = min(d_tof);

p0 = p(J_min,1);

tof_min = tof(1,:);
tof_max = tof(end,:);

if (tof_min<TOF&&tof_max>TOF) || (tof_min>TOF&&tof_max<TOF)
    a0 = M.*K.*p0./((2*M-L.^2).*p0.^2+2*K.*L.*p0-K.^2);
    f0 = 1-r2./p0.*(1-cos(DTA_L));
    dot_f0 = sqrt(mu./p0).*((1-cos(DTA_L))./p0-1./r1-1./r2).*tan(DTA_L/2);
    cosDE0 = 1-r1./a0.*(1-f0);
    sinDE0 = -r1.*r2.*dot_f0./sqrt(mu.*a0);
    DE0 = mod(atan2(sinDE0,cosDE0),2*pi);          
    
    PAR = [r1,r2,DTA_L,TOF*86400,mu,1];
    X0 = [a0,DE0,p0]';
    
    X = fsolve(@(X)func_Lambert(X,PAR),X0,Options);

    a_L  = X(1);
    DE_L = X(2);
    p_L  = X(3);
       
    f    = 1-r2./p_L.*(1-cos(DTA_S));
    g    = r1.*r2./sqrt(p_L.*mu).*sin(DTA_L);
    dot_f = sqrt(mu./p_L).*((1-cos(DTA_S))./p_L-1./r1-1./r2).*tan(DTA_L/2);
    dot_g = 1-r1./p_L.*(1-cos(DTA_S));

    v1_vec_L = (r2_vec - f.*r1_vec)./g;
    v2_vec_L = dot_f.*r1_vec + dot_g./g.*(r2_vec - f.*r1_vec);

    TOF_check = ((a_L.^3./mu).^.5.*(DE_L-sin(DE_L)) + r1.*r2./(p_L.*mu).^.5.*sin(DTA_L))/86400;  
    delta_TOF = max(abs(TOF - TOF_check));

    if delta_TOF>1E-7
        fprintf('LONG-Way  Orbit may not reliable for TOF = %d days (Error: %.1f days)\n\n',TOF,delta_TOF)
    end
    % Departure
    vinf_dep_vec_L = v1_vec_L - VP1_dep;
    vinf_dep_L = vecnorm(vinf_dep_vec_L);
    C3_dep_L = vinf_dep_L^2;

    % Arrival
    vinf_arr_vec_L = VP2_arr - v2_vec_L;
    vinf_arr_L = vecnorm(vinf_arr_vec_L);
    C3_arr_L = vinf_arr_L^2;

else
    C3_dep_L = NaN;
    C3_arr_L = NaN;
end

if isnan(C3_dep_S) && isnan(C3_dep_L)
    fprintf('\nTransfer physically unrealizable for a closed orbit around the Sun. ')
    fprintf('Select a different time of flight (TOF).\n')
    return
end
if C3_dep_S<C3_dep_L || isnan(C3_dep_L)
    orbit_choice = sprintf('              The SHORT-Way Orbit is the viable solution.\n');
    p = p_S;
    C3_dep = C3_dep_S;
    C3_arr = C3_arr_S;
    color_traj = [0,0.6,1];
    v1_vec = v1_vec_S;
    v2_vec = v2_vec_S;
elseif C3_dep_S>=C3_dep_L || isnan(C3_dep_S)
    orbit_choice = sprintf('              The  LONG-Way Orbit is the viable solution.\n');
    p = p_L;
    C3_dep = C3_dep_L;
    C3_arr = C3_arr_L;
    color_traj = [.9,0,.9];
    v1_vec = v1_vec_L;
    v2_vec = v2_vec_L;
end


[e,a,i,w,Om,TA1,gam1,h] = func_RVtoOE(r1_vec,v1_vec,mu);
[~,~,~,~, ~,TA2,gam2,~] = func_RVtoOE(r2_vec,v2_vec,mu);

% Orbital Elements of Planet 1 at Departure
[e_P1,a_P1,i_P1,w_P1,Om_P1,~,~,h_P1] = func_RVtoOE(r1_vec,VP1_dep,mu);
% Orbital Elements of Planet 2 at Arrival
[e_P2,a_P2,i_P2,w_P2,Om_P2,~,~,h_P2] = func_RVtoOE(r2_vec,VP2_arr,mu);

%fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
%fprintf('        Interplanetary %s-to-%s Orbit Characteristics:\n',Planet_1,Planet_2)
%fprintf('\n%s',orbit_choice)
%fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
%fprintf('       Semi-Major Axis:    a  = %.1f km\n',a)
%fprintf('          Eccentricity:    e  = %.4f\n'   ,e)
%fprintf('      Periapsis Radius:   rp = %.1f km\n',a*(1-e))
%fprintf('       Apoapsis Radius:   ra = %.1f km\n',a*(1+e))
%fprintf('                Period:    T = %.1f days\n\n',2*pi/sqrt(mu)*a^1.5/86400)

%fprintf('R.A. of Ascending Node:   Om = %.1f  deg\n',Om)
%fprintf('           Inclination:    i = %.1f  deg\n',i)
%fprintf(' Argument of Periapsis:    w = %.1f  deg\n',w)

%fprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ \n')

%fprintf('True Anomaly of S/C at Departure:  TA_1 = %3.1f  deg\n'  ,TA1)
%fprintf('True Anomaly of S/C at   Arrival:  TA_2 = %3.1f  deg\n',TA2)
%fprintf('              Delta True Anomaly:         %3.1f  deg\n',mod(TA2-TA1,360))
%fprintf('                  Time of Flight:   TOF = %3.0f  days\n\n',TOF)

%fprintf('At Departure:  v_inf = %.3f  km/s, C3 = %.1f km^2/s^2\n',sqrt(C3_dep),C3_dep)
%fprintf('At   Arrival:  v_inf = %.3f  km/s, C3 = %.1f km^2/s^2\n',sqrt(C3_arr),C3_arr)

%% S/C Orbit
R_PI = func_RPI(Om,i,w);
R_IP = R_PI';

if TA1<TA2
    TA   = TA1:0.1:TA2;
    r_PF = zeros(3,size(TA,2));
    r_I  = zeros(3,size(TA,2));
    for j=1:size(TA,2)
        r_PF(:,j) = p/(1+e*cosd(TA(j)))*[cosd(TA(j)),sind(TA(j)),0]';
        r_I(:,j)  = R_IP*r_PF(:,j);
    end
else
    TAa   = 0:0.1:TA2;
    r_PFa = zeros(3,size(TAa,2));
    r_Ia  = zeros(3,size(TAa,2));
    for j=1:size(TAa,2)
        r_PFa(:,j) = p/(1+e*cosd(TAa(j)))*[cosd(TAa(j)),sind(TAa(j)),0]';
        r_Ia(:,j)  = R_IP*r_PFa(:,j);
    end
    TAb = TA1:0.1:360;
    r_PFb = zeros(3,size(TAb,2));
    r_Ib  = zeros(3,size(TAb,2));
    for j=1:size(TAb,2)
        r_PFb(:,j) = p/(1+e*cosd(TAb(j)))*[cosd(TAb(j)),sind(TAb(j)),0]';
        r_Ib(:,j)  = R_IP*r_PFb(:,j);
    end
end

%% Planet_1 Orbit
color_P1 = [0,1,0];
R_PI_P1 = func_RPI(Om_P1,i_P1,w_P1);
R_IP_P1 = R_PI_P1';

TA_P1   = 0:0.1:360;
r_PF_P1 = zeros(3,size(TA_P1,2));
r_I_P1  = zeros(3,size(TA_P1,2));
for j=1:size(TA_P1,2)
    r_PF_P1(:,j) = (h_P1^2/mu)/(1+e_P1*cosd(TA_P1(j)))*[cosd(TA_P1(j)),sind(TA_P1(j)),0]';
    r_I_P1(:,j)  = R_IP_P1*r_PF_P1(:,j);
end

%% Planet_2 Orbit
color_P2 = [1,.6,0];
R_PI_P2 = func_RPI(Om_P2,i_P2,w_P2);
R_IP_P2 = R_PI_P2';

TA_P2   = 0:0.1:360;
r_PF_P2 = zeros(3,size(TA_P2,2));
r_I_P2  = zeros(3,size(TA_P2,2));
for j=1:size(TA_P2,2)
    r_PF_P2(:,j) = (h_P2^2/mu)/(1+e_P2*cosd(TA_P2(j)))*[cosd(TA_P2(j)),sind(TA_P2(j)),0]';
    r_I_P2(:,j)  = R_IP_P2*r_PF_P2(:,j);
end

%% Plots
figure(1),hold on
set(gcf,'units','normalized','position',[0,0.06,0.5,0.83]) 
set(gca,'Color',[0,0,0])
set(gcf,'Color',[0,0,0])
view([Om,25])
strt1 = sprintf('Launch Date: %d-%d-%d,    Arrival Date: %d-%d-%d,    [TOF: %.0f days]\n',yyyy,mm,dd,yyyyA,mmA,ddA,TOF);
strt2 = sprintf('C3(dep.)=%.1f km^2/s^2,  v_{inf}(arr.)=%.2f km/s',C3_dep,sqrt(C3_arr));
tlt_t = title([strt1,strt2]);
tlt_t.Color = 'w'; 
tlt_t.FontSize = 13; 
l_S = plot3(0,0,0,'markeredgecolor',[1,1,0.4],'markerfacecolor',[1,1,0.4],'marker','o','markersize',20,'color','none');
if TA1<TA2
    l = plot3(r_I(1,:),r_I(2,:),r_I(3,:),'color',color_traj,'linewidth',1.5);
else
    l = plot3(r_Ia(1,:),r_Ia(2,:),r_Ia(3,:),'color',color_traj,'linewidth',1.5);
        plot3(r_Ib(1,:),r_Ib(2,:),r_Ib(3,:),'color',color_traj,'linewidth',1.5);
end
plot3(r_I_P1(1,:),r_I_P1(2,:),r_I_P1(3,:),'color',color_P1,'linewidth',0.5);
plot3(r_I_P2(1,:),r_I_P2(2,:),r_I_P2(3,:),'color',color_P2,'linewidth',0.5);
l_P1     = plot3(r1_vec(1,1),r1_vec(2,1),r1_vec(3,1),'markeredgecolor',color_P1,'markerfacecolor',color_P1,'marker','o','markersize',8,'color',color_P1);
l_P2     = plot3(r2_vec(1,1),r2_vec(2,1),r2_vec(3,1),'markeredgecolor',color_P2,'markerfacecolor',color_P2,'marker','o','markersize',8,'color',color_P2);
text_P1 = sprintf('   %s on %d-%d-%d',Planet_1,yyyy ,mm ,dd);
text_P2 = sprintf('   %s on %d-%d-%d',Planet_2,yyyyA,mmA,ddA);
AU = 1.49597871E8;
text(r1_vec(1,1)+AU/10,r1_vec(2,1)+AU/10,r1_vec(3,1)+AU/10,text_P1,'Color','w')%,'FontSize',16,'interpreter','latex');
text(r2_vec(1,1)+AU/10,r2_vec(2,1)+AU/10,r2_vec(3,1)+AU/10,text_P2,'Color','w')%,'FontSize',16,'interpreter','latex');
    
axis('equal')

l = legend([l_S,l_P1,l_P2,l],'Sun',sprintf('%s',Planet_1),sprintf('%s',Planet_2),'Interplanetary Transfer Orbit');
l.TextColor=[1 1 1];
end

function F = func_Lambert(X,PAR)
    J   = PAR(end);
    r1  = PAR(0*J+1:1*J)';
    r2  = PAR(1*J+1:2*J)';
    DTA = PAR(2*J+1:3*J)';
    Dt  = PAR(3*J+1:4*J)';
    mu  = PAR(end-1);
    
    a  = X(0*J+1:1*J);
    DE = X(1*J+1:2*J);
    p  = X(2*J+1:3*J);

    AU = 1.49597871E8;
    TU = 5022643;

    F = [(p/AU).*(a/AU).*(1-cos(DE)) - (r1/AU).*(r2/AU).*(1-cos(DTA));
        (Dt/TU).*((p/AU).*(mu/AU^3*TU^2)).^.5 - ((p/AU).*(a/AU).^3).^.5.*(DE-sin(DE)) - (r1/AU).*(r2/AU).*sin(DTA)
        ((r1/AU).*(r2/AU).*(1-cos(DTA)) - (p/AU).*((r1/AU)+(r2/AU))).*tan(DTA/2) + ((p/AU).^3.*(a/AU)).^.5.*sin(DE)];
end
function [r,R,mu,T,n,v,r_SOI] = func_parameters_Planet(Planet)

    % Planetary Parameters
    mu_Sun = 1.3271244e11; %[km^3/s^2]
    
    if strcmp(Planet,'Mercury')
        R = 2440;
        mu = 22032;
        r = 57.91e6;
    elseif strcmp(Planet,'Venus')
        R = 6052;
        mu = 324900;
        r = 108.2e6;
    elseif strcmp(Planet,'Earth')    
        R = 6378;
        r = 149.6E6;
        mu = 398600;
    elseif strcmp(Planet,'Mars')   
        R = 3396;
        mu = 42828;
        r = 227.9e6;
    elseif strcmp(Planet,'Jupiter')
        R = 71490;
        mu = 126686000;
        r = 778.6e6;
    elseif strcmp(Planet,'Saturn') 
        R = 60270;
        mu = 37931000;
        r = 1433e6;
    elseif strcmp(Planet,'Uranus') 
        R = 25560;
        mu = 5794000;
        r = 2872e6;
    elseif strcmp(Planet,'Neptune')
        R = 24760;
        mu = 6835100;
        r = 4495e6; 
    elseif strcmp(Planet,'Pluto')
        R = 1188;
        mu = 869.61;
        r = 5.90638E9;
    end
    
    T = 2*pi/sqrt(mu_Sun)*r^1.5; %[s]
    n = 2*pi/T;                  %[rad/s]
    v = sqrt(mu_Sun/r);          %[km/s]
    r_SOI = r*(mu/mu_Sun)^(2/5); %[km]
end