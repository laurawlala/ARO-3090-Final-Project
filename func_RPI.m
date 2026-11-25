function R_PI = func_RPI(Om,i,w)

% Computes the rotation matrix from Earth-Centered Inertial RF (I) to Perifocal RF (P).

ROm = [ cosd(Om),sind(Om), 0;
       -sind(Om),cosd(Om), 0;
               0,       0, 1];
Ri  = [1,       0,      0; 
       0, cosd(i),sind(i);
       0,-sind(i),cosd(i)];

Rw  = [ cosd(w),sind(w), 0;
       -sind(w),cosd(w), 0;
              0,      0, 1]; 
    
R_PI = Rw*Ri*ROm;