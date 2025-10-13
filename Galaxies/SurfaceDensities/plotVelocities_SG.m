function maxy = plotVelocities_SG(r_data, V_obs, V_err,r_model,V_model,V_tot,V_MOND,V_ele,V_ele_G,GalaxyName)
%plotVelocities_SG Summary of this function goes here
%   Detailed explanation goes here

    % Plot results
    plot(r_model, V_model, 'r-', 'LineWidth', 2, 'DisplayName', 'Smooth');
    hold on;
    plot(r_model, V_tot, 'g-', 'LineWidth', 4, 'DisplayName', 'Stars+gas');
    plot(r_model, V_MOND, '-', 'LineWidth', 4, 'DisplayName', 'MOND');
    plot(r_model, V_ele, 'b+', 'LineWidth', 1, 'DisplayName', 'E.Vel');
    plot(r_model, V_ele_G, 'b.', 'LineWidth', 2, 'DisplayName', 'EG.Vel');
    errorbar(r_data, V_obs, V_err, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');

    xlabel('Radius (kpc)');
    ylabel('Velocity (km/s)');
    titles = append('Rotation Curve: ',GalaxyName);
    title(titles);

    legend('Location', 'southwest','NumColumns',2);
    grid on;
    grid minor;
    maxy = V_obs + 3*V_err;
    maxy = max([maxy;V_tot;V_MOND;V_ele;V_ele_G]);
    xlim([0, max(r_data)]);
    ylim([0, maxy]); % Enforce physical limits
    hold off

end