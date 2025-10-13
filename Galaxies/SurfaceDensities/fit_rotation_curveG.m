function [V_tot,V_MOND,V_ele_G,V_ele,Zdist,StarProfile,GasZprofile,TAcel_x,TAcel_y]=fit_rotation_curveG(B_sigma,A_Sigma_star,A_Sigma_H,r_model,V_model,MStars,Mgas,sig_stars,sig_H,resolution,lpidx,GalaxyName)

    c = 299792.458; % Speed of light in vacuum (Km/s)
    sig_stars=sig_stars/c;
    sig_H=sig_H/c;

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

    V_tot = sqrt(G * (M_H+M_R) ./ R_pc);  % km/s


    % Newtonian acceleration
    a_N = G * (M_H+M_R) ./ (R_pc.^2);          % km^2/s^2 / pc
    

    % MOND correction: inverse interpolating function (standard form)
    y = a0*3.086e10./a_N;
    a_MOND = a_N.*sqrt(0.5 + sqrt(0.25 + y.^2));  
    % Rotational velocity
    V_MOND = sqrt(a_MOND .* R_pc);       % km/s

    % Electromagnetic G Force

    velr = V_model;
    velr(1) = 0;

    TAcel_x = zeros(size(R_pc));
    TAcel_y = zeros(size(R_pc));
    TAcel0_x = zeros(size(R_pc));
    TAcel0_y = zeros(size(R_pc));
    if lpidx==0
        plot(r_model,TAcel_x, 'r-')
        titles = append('Radial Aceleration: ',GalaxyName);
        title(titles);
        grid on;
        grid minor;
        xlabel('Radius (kpc)');
        ylabel('Aceleration (m/s^2)');
        hold on;
    else
        h = waitbar(0, 'Please wait...'); % Initialize waitbar with 0% progress and a message
    end

    for i = 2:(length(R_pc)-1)
        [gAcel_x,gAcel_y,gAcel0_x,gAcel0_y,Zdist,StarProfile,GasZprofile] = galaticRingGR_vec(sig_stars,sig_H,B_sigma,A_Sigma_star,A_Sigma_H,i,R_pc,velr,MStars,Mgas,resolution);
        TAcel_x = TAcel_x +  gAcel_x;
        TAcel_y = TAcel_y +  gAcel_y;
        TAcel0_x = TAcel0_x +  gAcel0_x;
        TAcel0_y = TAcel0_y +  gAcel0_y;
        if lpidx==0
            plot(r_model,TAcel_x, 'r-')
            plot(r_model(1:i),TAcel_x(1:i), 'b.', 'MarkerSize', 5);
            plot(r_model(1:i),TAcel0_x(1:i), 'g.', 'MarkerSize', 3);
            drawnow;
        else
            waitbar(i / length(R_pc), h, sprintf('Processing: %d: %d%%, ',lpidx,round(i/length(R_pc)*100))); % Update waitbar
        end
    end
    if lpidx==0
        p2=plot(r_model,TAcel_x, 'b-', 'MarkerSize', 8);
        p1=plot(r_model,TAcel0_x, 'g-', 'MarkerSize', 4);
        hold off;
        legend([p1 p2],{'Static Stars', 'Dynamic Stars + Gas'},'Location','southeast');
    else
        close(h)
    end
    TAcel_x(1)= 0;
    TAcel_y(1)= 0;
    TAcel0_x(1)= 0;

    TAcel_x(isnan(TAcel_x)) = TAcel0_x(isnan(TAcel_x));
    TAcel_x(TAcel_x>0) = 0;
    TAcel_x(isnan(TAcel_x)) = 0;

    V_ele = sqrt(abs(TAcel_x).*R_pc*3.086e10);
    V_ele_G = sqrt(abs(TAcel0_x).*R_pc*3.086e10);

end
