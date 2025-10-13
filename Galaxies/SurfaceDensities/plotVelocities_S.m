function figO = plotVelocities_S(gloc,r_data, V_obs, V_err,R_kpc, V_tot,V_MOND,r_model,V_model,V_ele,V_ele_G)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if (gloc>0)
        gain = V_model(gloc)/evel(gloc)
        evel = gain*evel;
    
        gain = V_model(gloc)/evelG(gloc)
        evelG = gain*evelG;
    end

    % Plot results
    errorbar(r_data, V_obs, V_err, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'M33 Data');
    hold on;
    plot(r_model, V_model, 'r-', 'LineWidth', 2, 'DisplayName', 'Best Fit');
    plot(R_kpc, V_tot, 'g-', 'LineWidth', 4, 'DisplayName', 'Stars+gas');
    plot(R_kpc, V_MOND, '-', 'LineWidth', 4, 'DisplayName', 'MOND');
    plot(r_model, V_ele, 'b+', 'LineWidth', 1, 'DisplayName', 'E.Vel');
    plot(r_model, V_ele_G, 'b.', 'LineWidth', 2, 'DisplayName', 'EG.Vel');

    xlabel('Radius (kpc)');
    ylabel('Velocity (km/s)');
    title('Rotation Curve');
    legend('Location', 'southwest','NumColumns',2);
    grid on;
    grid minor;
    maxy = V_obs + 3*V_err;
    maxy = max([maxy;V_tot;V_MOND;V_ele;V_ele_G]);
    xlim([0, max(R_kpc)]);
    ylim([0, maxy]); % Enforce physical limits
    hold off

end