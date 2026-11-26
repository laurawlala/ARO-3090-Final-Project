% Launch Date

% TOF 1

% TOF 2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hohmann Transfer: Earth to Jupiter 
Planet_1 = 'Earth';
Plnet_2 = 'Jupiter';

% Fly By Date
yyyy = 2038;
mm = 1;
dd = 1;
TOF = 1388.8;

[V_inf_dep, V_inf_arr] = Lamberts_Interplanetary(Planet_1,Planet_2,yyyy,mm,dd,TOF); % V_inf at Arrival and Departure
V_inf_dep = V_inf_dep
V_inf_fb1 = V_inf_arr

% Hohmann Transfer: Jupiter to Saturn 
Planet_2= 'Jupiter';
Planet_3 = 'Saturn';

% Arrival Date
yyyy = 2038;
mm = 1;
dd = 1;
TOF = 1388.8;

[V_inf_dep, V_inf_arr] = Lamberts_Interplanetary(Planet_2,Planet_3,yyyy,mm,dd,TOF); % V_inf at Arrival and Departure
V_inf_fb2 = V_inf_dep
V_inf_arr = V_inf_arr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if abs(V_inf_fb1) == abs(V_inf_fb2)
    







































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rocket Equation

Isp = 300; % sec
m0 = 1000; % wet mass (kg)
g0 = 9.806E-3; % grav accel at sea level (km/sec^2)

%delta_m = m0(1-e^(-delta_Varr/(Isp*g0))); % Rocket Equation for mass of fuel (kg)







