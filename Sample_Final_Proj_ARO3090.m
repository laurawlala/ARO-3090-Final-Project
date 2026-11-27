% Assuming a flyby: Mars
% Vehicle: SLS Block 1 - 5m fairing

% Note: Anything bracketed by [%--- code  %---] is to be manually input
% from the Prof's lambert code!

load('Porkchop_Candidates.mat','Candidate')
load('SecondLeg_Lambert_Result.mat','secondLeg')
k = secondLeg.k;

fbPlanet = 'Mars';

C3_max = 136.67;             % [km^2/s^2] - from the space vehicle
launch_date =  datetime(2039, 06, 06, 10, 45, 46, 'TimeZone', 'UTC');
jd_launch = juliandate(launch_date);
TOF1 = 246 + (19/24);        % [days]

% Earth -> Flyby first leg
%--------------------------------------------
C3_Earth_dep = Candidate(k).C3_dep;
v_inf_Earth_dep_vec = Candidate(k).vinf_dep_vec; % From porkchop sweep code
mag_Earth_dep = norm(v_inf_Earth_dep_vec)

C3_fb_arr = Candidate(k).C3_arr;
v_inf_fb_arr_vec = Candidate(k).vinf_arr_vec; % From porkchop sweep code
mag_fb_arr = norm(v_inf_fb_arr_vec)
%--------------------------------------------

% Flyby -> Saturn second leg
%--------------------------------------------
C3_fb_dep = secondLeg.C3_dep;
v_inf_fb_dep_vec = secondLeg.vinf_dep_vec;
mag_fb_dep = norm(v_inf_fb_dep_vec) 

C3_Saturn_arr = secondLeg.C3_arr; % Methinks this needs to be changed but its 1am and Im tireeed!!
v_inf_Saturn_arr_vec = secondLeg.vinf_arr_vec; % From porkchop sweep code
mag_Saturn_arr = norm(v_inf_Saturn_arr_vec)
%--------------------------------------------

tol = 1e-4
if abs(mag_fb_dep - mag_fb_arr) <= tol
    % Checking if vinf_fb_out_vec = -vinf_fb_in_vec
    disp('FEASIBLE! C3 in/out of flyby check good!')
else
    disp('NOT FEASIBLE! The C3 in and out of Flyby aint da same')
    return;
end

TOF2 = 2920 - TOF1;                % [days]

%--------------------------------------------
C3_Saturn_arr = 36.11; 
v_inf_Saturn_arr_vec = [-4.32099013; 0.61569784; 0.14696635];
mag_Saturn_arr = norm(v_inf_Saturn_arr_vec) % 4.367
%--------------------------------------------

[r,R_fb,mu,T,n,v,r_SOI,zp_min] = func_parameters_Planet(fbPlanet)

% Launch vehicle check

if C3_Earth_dep <= C3_max
    disp('FEASIBLE! Launch vehicle can support dep C3')
else
    disp('NAWT FEASIBLE!! Launch vehicle cannot support the departure C3')
    return;
end

vinf = (mag_fb_dep + mag_fb_arr)/2; % They should equal about the same so...
turn_angle = acosd(dot(v_inf_fb_dep_vec, v_inf_fb_arr_vec)/vinf^2); % is that the right way to dot product?
beta = (180 - turn_angle)/2;
e = 1/cosd(beta);
theta_inf = acosd(-1/e);

energy = vinf^2/2;
a = mu/(2*energy);
h = sqrt(a*mu*(e^2 - 1));
rp = a*(e - 1);
zp = rp - R_fb;

fprintf('Eccentricity = %.5f\n',e)
fprintf('Theta infinity = %.3f deg\n',theta_inf)
fprintf('Specific energy = %.5f km^2/s^2\n',energy)
fprintf('Semi-major axis = %.3f km\n',a)
fprintf('Angular momentum = %.3f km^2/s\n',h)
fprintf('Periapsis radius = %.3f km\n',rp)
fprintf('Periapsis altitude = %.3f km\n',zp)
fprintf('Minimum required altitude = %.3f km\n',zp_min)

if zp >= zp_min
    disp('Feasible Transfer')
else 
    disp ('NOT Feasible')
    return;
end



function [r,R_fb,mu,T,n,v,r_SOI,zp_min] = func_parameters_Planet(fbPlanet)
    % Planetary Parameters
    mu_Sun = 1.3271244e11; %[km^3/s^2]    
    if strcmp(fbPlanet,'Mercury')
        R_fb = 2440;
        mu = 22032;
        r = 57.91e6;
        zp_min = 150;
    elseif strcmp(fbPlanet,'Venus')
        R_fb = 6052;
        mu = 324900;
        r = 108.2e6;
        zp_min = 275;
    elseif strcmp(fbPlanet,'Earth')    
        R_fb = 6378;
        r = 149.6E6;
        mu = 398600;
        zp_min = 300;
    elseif strcmp(fbPlanet,'Mars')   
        R_fb = 3396;
        mu = 42828;
        r = 227.9e6;
        zp_min = 200;
    elseif strcmp(fbPlanet,'Jupiter')
        R_fb = 71490;
        mu = 126686000;
        r = 778.6e6;
        zp_min = 100000;
    elseif strcmp(fbPlanet,'Saturn') 
        R_fb = 60270;
        mu = 37931000;
        r = 1433e6;
        zp_min = 50000;
    elseif strcmp(fbPlanet,'Uranus') 
        R_fb = 25560;
        mu = 5794000;
        r = 2872e6;
        zp_min = 50000;
    elseif strcmp(fbPlanet,'Neptune')
        R_fb = 24760;
        mu = 6835100;
        r = 4495e6; 
        zp_min = 40000;
    elseif strcmp(fbPlanet,'Pluto')
        R_fb = 1188;
        mu = 869.61;
        r = 5.90638E9;
    end
    
    T = 2*pi/sqrt(mu_Sun)*r^1.5; %[s]
    n = 2*pi/T;                  %[rad/s]
    v = sqrt(mu_Sun/r);          %[km/s]
    r_SOI = r*(mu/mu_Sun)^(2/5); %[km]
end
