function [r_data,V_obs,V_err,R_kpc,r_model,V_model,SigmaS,SigmaG,MStars,Mgas] = SRotation_curves(Observations,GalaxyName)
    % Load M33 empirical data (replace with actual data)
    % Format: [Radius (kpc), Velocity (km/s), Error (km/s)]
    % M33 Rotation Curve Data from Corbelli et al. (2014)
    % Columns: [Radius (kpc), Velocity (km/s), Velocity Error (km/s), 
    %          HI Surface Density (M_sun/pc^2), Stellar Surface Density (M_sun/pc^2)]
    if (isempty(Observations))
    Observations = [
        0.24    37.3    6.2    7.12    316.59;
        0.28    37.9    5.5    7.09    292.13;
        0.46    47.1    3.9    6.76    230.46;
        0.64    53.5    2.8    6.38    184.39;
        0.73    55.1    3.3    6.32    164.93;
        0.82    58.5    2.5    6.22    151.18;
        1.08    66.2    3.3    5.84    124.13;
        1.22    69.4    2.8    5.52    111.64;
        1.45    74.6    5.0    5.75    93.78;
        1.71    77.9    2.9    6.12    77.00;
        1.87    81.7    2.5    6.46    68.51;
        2.20    86.8    3.0    6.98    57.80;
        2.28    90.1    2.1    7.07    55.47;
        2.69    94.4    3.9    6.27    44.91;
        2.70    95.4    1.8    6.26    44.68;
        3.12    99.2    2.6    6.33    35.99;
        3.18    98.7    2.0    6.29    34.89;
        3.53    101.3    2.4    6.77    29.14;
        3.66    101.5    1.3    6.83    27.25;
        4.15    106.3    2.8    6.69    21.17;
        4.64    109.4    1.4    6.47    16.45;
        5.13    108.8    1.6    6.19    12.78;
        5.62    107.3    2.7    6.13    9.93;
        6.11    108.2    2.2    5.94    7.72;
        6.60    109.8    1.6    5.62    6.00;
        7.09    110.1    1.7    5.14    4.66;
        7.57    111.1    1.6    4.35    3.64;
        8.06    113.0    2.6    3.20    2.83;
        8.55    113.9    2.4    2.41    2.20;
        9.04    115.1    3.2    1.78    1.71;
        9.53    116.3    3.1    1.26    1.33;
        10.02   119.1    2.7    0.99    1.04;
        10.51   121.0    2.7    0.82    0.96;
        10.99   121.5    2.1    0.72    0.89;
        11.48   118.6    2.0    0.65    0.82;
        11.97   118.7    1.9    0.57    0.76;
        12.46   117.2    1.9    0.48    0.70;
        12.95   116.2    2.0    0.42    0.65;
        13.44   118.3    2.4    0.36    0.60;
        13.93   119.0    2.4    0.34    0.55;
        14.41   121.3    2.9    0.33    0.51;
        14.90   121.4    3.7    0.33    0.47;
        15.39   120.3    3.4    0.26    0.44;
        15.88   121.9    4.2    0.23    0.40;
        16.37   126.3    7.2    0.22    0.37;
        16.86   126.3    4.5    0.21    0.34;
        17.35   127.2    5.5    0.15    0.32;
        17.84   126.2    4.7    0.15    0.29;
        18.32   124.2    12.2    0.11    0.27;
        18.81   127.2    13.3    0.14    0.25;
        19.30   120.2    14.2    0.15    0.23;
        19.79   121.8    13.2    0.17    0.22;
        20.28   136.0    19.1    0.12    0.20;
        20.77   128.3    14.3    0.10    0.18;
        21.26   127.4    15.4    0.09    0.17;
        21.74   120.1    13.1    0.09    0.16;
        22.23   112.2    14.9    0.08    0.15;
        22.72   119.6    13.4    0.07    0.13;
    ];
    end

    r_data = Observations(:,1);  % Radius [kpc]
    V_obs  = Observations(:,2);  % Observed velocity [km/s]
    V_err  = Observations(:,3);  % Error [km/s]
    

    SigmaS = [0;Observations(:,5)];         % Stellar surface density (M_sun/pc^2)
    SigmaG = [0;Observations(:,4)];         % Gas surface density (M_sun/pc^2)
    SigmaS(1) = 1.05*SigmaS(2);
    SigmaG(1) = SigmaG(2);

    npts = 200;
    R_kpc = [0.0001;r_data];          % Radius in kpc
    V_obs = [0;V_obs];
    V_err = [0;V_err];
    r_data = R_kpc;
    r_model = linspace(0.0001, 1.1*max(R_kpc), npts)'; % [kpc]
    V_model = V_obs;
    V_model(R_kpc>3) = smooth(V_obs(R_kpc>3),0.05);
    V_model(R_kpc<=3) = V_obs(R_kpc<=3);
    V_model(end) = V_obs(end);
    V_lin = [V_model;0];
    rmodlin = [r_data;10*max(r_data)];
    V_model = interp1(r_data,V_model,r_model,'spline');
    V_model_l = interp1(rmodlin,V_lin,r_model);
    V_model(floor(0.91*npts):npts)=V_model_l(floor(0.91*npts):npts);
    V_model(1) = V_obs(1);
    V_model(end) = V_model_l(end);
    V_model(r_model>3) = smooth(V_model(r_model>3),0.05,'lowess');
    V_model(r_model<=3) = V_model_l(r_model<=3);
    V_model(end) = V_model_l(end);
    V_model(npts/2:npts) = smooth(V_model(npts/2:npts),55,'sgolay',2);
    V_model(r_model>10) = smooth(V_model(r_model>10),55,'sgolay',2);
    V_model(end) = V_model_l(end);
    V_model(floor(3*npts/4):npts) = smooth(V_model(floor(3*npts/4):npts),25,'sgolay',1);
    V_model(end) = V_model_l(end);
%    V_model(r_model>3) = smooth(V_model(r_model>3),11,'sgolay',2);
    V_model = smooth(V_model,0.05,"lowess");
    V_model(end) = V_model_l(end);
    V_model(1) = V_obs(1);

    
    % Convert units: kpc -> pc, surface density to consistent units
    R_pc = R_kpc * 1000;                 % Convert kpc to pc
    

    % Get smooth estimations of the mass densities
    distp = 1000*r_model;
    distp(1) = 0;
    LSigmaS = log10(SigmaS);
    mz = LSigmaS(2);
    MStars = smooth(LSigmaS,0.05);
    MStars(1) = mz;
%    MStars = interp1(R_pc,MStars,distp,'spline');
    MStars = interp1(R_pc,MStars,distp);
    MStars(1) = mz;
    MStars(npts/2:npts) = smooth(MStars(npts/2:npts),95,'sgolay',1);
    MStars(1:floor(npts/4)) = smooth(MStars(1:floor(npts/4)),25,'sgolay',1);
    MStars(floor(npts/4):npts/2) = smooth(MStars(floor(npts/4):npts/2),25,'sgolay',1);
    MStars(1) = mz;
    MStars=smooth(MStars,25,'sgolay',2);
    MStars(1) = mz;
    MStars = 10.^MStars;

    mz = SigmaG(1);
    LSigmaG = SigmaG;
    LSigmaG = smooth(LSigmaG,0.05);
    LSigmaG(1)=mz;
    LSigmaG = log10(LSigmaG);
    mz = LSigmaG(1);
    Mgas = LSigmaG;
    lastgas = LSigmaG(end);
    Mgas = smooth(LSigmaG,0.10,'rlowess') ;
    Mgas(1)=mz;
    Mgas(end) = lastgas;
%    Mgas = interp1(R_pc,Mgas,distp,'spline');
    Mgas = interp1(R_pc,Mgas,distp);
    Mgas(1)=mz;
    Mgas(npts)=lastgas;
%    Mgas=smooth(Mgas,35,'sgolay',2);
%    Mgas(1) = mz;
    Mgas(npts/2:npts) = smooth(Mgas(npts/2:npts),95,'sgolay',1);
    Mgas(1:floor(npts/3)) = smooth(Mgas(1:floor(npts/3)),25,'sgolay',2);
    Mgas(floor(npts/3):npts/2) = smooth(Mgas(floor(npts/3):npts/2),35,'sgolay',2);
    Mgas(1) = mz;
    Mgas=smooth(Mgas,15,'sgolay',1);
    Mgas(1) = mz;
    Mgas = 10.^Mgas;

    % Plot the smoothed densities

    semilogy(R_kpc,SigmaS,'r.', 'LineWidth', 2, 'DisplayName', 'Stars')
    hold on
    semilogy(R_kpc,SigmaG,'b+', 'LineWidth', 2, 'DisplayName', 'Gas')
    semilogy(distp/1000,MStars,'r-', 'LineWidth', 1, 'DisplayName', 'Stars')
    semilogy(distp/1000,Mgas,'b-', 'LineWidth', 1, 'DisplayName', 'Gas')
    xlabel('Radius (kpc)');
    ylabel('Surface Density (M_sun/pc^2)');
    titles = append('Surface Density: ',GalaxyName);
    title(titles);
    legend('Location', 'northeast','NumColumns',2);
    grid on;
    grid minor;
    hold off

end
