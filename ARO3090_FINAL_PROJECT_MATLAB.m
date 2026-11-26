% Constants

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
V_inf_fb_1 = V_inf_arr

% Hohmann Transfer: Jupiter to Saturn 
Planet_1= 'Jupiter';
Planet_2 = 'Saturn';

% Arrival Date
yyyy = 2038;
mm = 1;
dd = 1;
TOF = 1388.8;

[V_inf_dep, V_inf_arr] = Lamberts_Interplanetary(Planet_1,Planet_2,yyyy,mm,dd,TOF); % V_inf at Arrival and Departure
V_inf_fb_2 = V_inf_dep
V_inf_arr = V_inf_arr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rocket Equation

Isp = 300; % sec
m0 = 1000; % wet mass (kg)
g0 = 9.806E-3; % grav accel at sea level (km/sec^2)

%delta_m = m0(1-e^(-delta_Varr/(Isp*g0))); % Rocket Equation for mass of fuel (kg)







