function [V_tot,V_MOND]=Vel_Newton_MOND(r_model,MStars,Mgas)


    % Convert units: kpc -> pc, surface density to consistent units
    R_pc = r_model * 1000;                 % Convert kpc to pc

    M_R = zeros(size(R_pc));
    M_H = zeros(size(R_pc));
    
    for i = 2:length(R_pc)
        integrand = 2 * pi * R_pc(1:i) .* MStars(1:i);
        M_R(i) = trapz(R_pc(1:i), integrand);
        integrand = 2 * pi * R_pc(1:i) .* Mgas(1:i);
        M_H(i) = trapz(R_pc(1:i), integrand);
    end
    G = 4.3009e-3;  % Gravitational constant [pc (km/s)^2 / M_sun]
    a0 = 1.2e-10;     % MOND a0 in m/s

    V_tot = sqrt(abs(G * (M_H+M_R) ./ R_pc));  % km/s


    % Newtonian acceleration
    a_N = G * (M_H+M_R) ./ (R_pc.^2);          % km^2/s^2 / pc
    

    % MOND correction: inverse interpolating function (standard form)
    y = a0*3.086e10./a_N;
    a_MOND = a_N.*sqrt(0.5 + sqrt(0.25 + y.^2));  
    % Rotational velocity
    V_MOND = sqrt(abs(a_MOND .* R_pc));       % km/s
    V_MOND(isnan(V_MOND)) = 0;
    V_tot(isnan(V_MOND)) = 0;


end
