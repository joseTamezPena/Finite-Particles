function [profile_data,Mstar,MLgain] = compute_surface_density(galaxy_name, massmodels_file, stellarmass_file,gasMassCorrection)
% COMPUTE_SURFACE_DENSITY Computes surface mass densities of stars and gas for a given galaxy.
%
% Inputs:
%   galaxy_name - String, name of the galaxy (e.g., 'UGCA442').
%   massmodels_file - String, path to MassModels_Lelli2016c.mrt file.
%   stellarmass_file - String, path to wise_ii.table1 file.
%
% Outputs:
%   Sigma_stars - Column vector, stellar surface mass density (Msun/pc^2) at each R.
%   Sigma_gas - Column vector, gas surface mass density (Msun/pc^2) at each R.
%   R - Column vector, radial positions (kpc).
%
% Assumptions:
%   - MassModels file format: ASCII with blocks per galaxy. Galaxy name on a line,
%     followed by data lines with 8 columns: R(kpc), Vobs(km/s), e_Vobs(km/s),
%     Vgas(km/s), Vdisk(Υ=1)(km/s), Vbul(Υ=1)(km/s), Sdisk(Lsun/pc^2), Sbul(Lsun/pc^2).
%     Lines starting with '#' are comments.
%   - Stellar mass file format: ASCII table with columns 'Galaxy' (string) and 'logMstar' (log10(Mstar/Msun)).
%   - Surface brightness Sdisk and Sbul are face-on values.
%   - Total luminosity L_total computed via trapezoidal integration: 2*pi * trapz(R, R .* (Sdisk + Sbul)).
%   - Υ = Mstar / L_total (Msun/Lsun).
%   - Sigma_stars = Υ * (Sdisk + Sbul).
%   - Sigma_gas approximated from Vgas using disk inversion formula: Sigma = (1/(2*pi*G*r)) * d(V^2 * r)/dr,
%     with units adjusted (r in pc, G = 4.30091e-3 pc Msun^{-1} (km/s)^2).
%   - Handles single galaxy block extraction.

% Constants
G = 4.30091e-3;  % pc Msun^{-1} (km s^{-1})^2

% Read stellar mass file
stellardata = readtable(stellarmass_file, 'Delimiter', '\t', 'ReadVariableNames', true);
gal_idx = strcmp(stellardata.Galaxy, galaxy_name);
if sum(gal_idx) == 0
    error('Galaxy %s not found in stellar mass file.', galaxy_name);
end
logMstar = stellardata.logMstar(gal_idx);
MLgain = stellardata.M_L(gal_idx);

Mstar = 10^logMstar;  % Total stellar mass in Msun

% Read mass models file and extract galaxy block
fid = fopen(massmodels_file, 'r');
if fid == -1
    error('Cannot open mass models file: %s', massmodels_file);
end

profile_data = [];
found_galaxy = false;
while ~feof(fid)
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end
    tline = strtrim(tline);
    if isempty(tline) || startsWith(tline, '#')
        continue;
    end
    parts = strsplit(tline);
    nameG = parts{1};
    found_galaxy = false;
    if strcmp(nameG, galaxy_name)
        found_galaxy = true;
    end
    if found_galaxy
%        parts = strsplit(tline);
%        length(parts)
        if length(parts) == 10
            profile_data(end+1, :) = str2double(parts);
        end
    end
end
fclose(fid);
if isempty(profile_data)
    if ~found_galaxy
        error('Galaxy %s not found in mass models file.', galaxy_name);
    else
        error('No data lines found for galaxy %s.', galaxy_name);
    end
end

% Assign columns: [Dis, R,Vobs, eVobs, Vgas, Vdisk, Vbul, Sdisk, Sbul]
R = profile_data(:, 3);  % kpc
Vobs = profile_data(:, 4);
eVobs = profile_data(:, 5);
Vgas = profile_data(:, 6);
Vdisk = profile_data(:, 7);
Vbul = profile_data(:, 8);
Sdisk = profile_data(:, 9);  % Lsun/pc^2
Sbul = profile_data(:, 10);  % Lsun/pc^2
if (R(1)>0.1)
    profile_data = [[0 0 0 0 0 0 0 0 0 0];profile_data];
    R = [R(1)/2;R];
    slp = 0.5;
    Vobs = [0.5*(slp*Vobs(1) + Vobs(2));Vobs];
    eVobs = [slp*eVobs(1);eVobs];
    Vgas = [slp*Vgas(1);Vgas];
    Vdisk = [slp*Vdisk(1);Vdisk];
    Vbul = [slp*Vbul(1);Vbul];
    Sdisk = [Sdisk(1);Sdisk];
    Sbul = [Sbul(1);Sbul];
end

% Compute total stellar surface brightness
S = Sdisk + Sbul;  % Lsun/pc^2

Sigma_stars = MLgain * S;  % Msun/pc^2

% Compute gas and star surface mass density from Newton apox and disk approximation
r_pc = R * 1000;  % pc
if length(R) < 2
    Sigma_gas = zeros(size(R));  % Cannot compute derivative
    Sigma_Stars2 = zeros(size(R));  % Cannot compute derivative
else
    V2r = smooth(Vgas,0.1);
    V2r(1) = Vgas(1);
    V2r(end) = Vgas(end);
    V2r = V2r .^ 2;
    V2r(1) = Vgas(1)^2;
    V2r = V2r.* r_pc;

    dV2r_dr = smoothAbsGrad(V2r,r_pc);

    rden = r_pc;
    rden(1) = 0.5*(r_pc(1) + r_pc(2));


    Sigma_gas = dV2r_dr ./ (2 * pi * G * rden);  % Msun/pc^2


    V2s = (smooth(Vdisk,0.1) .^ 2 + smooth(Vbul,0.1).^2);
    V2s(1) = (Vdisk(1)^ 2 + Vbul(1)^2);
    V2s = V2s.* r_pc;

    dV2s_dr = smoothAbsGrad(V2s,r_pc);

    Sigma_Stars2 = dV2s_dr ./ (2 * pi * G * rden);  % Msun/pc^2
end
Sigma_stars(Sigma_stars <= 0 ) = 0.1*min(Sigma_stars(Sigma_stars>0));
Sigma_Stars2(Sigma_Stars2 <=0 ) = 0.1*min(Sigma_Stars2(Sigma_Stars2>0));
Sigma_gas(Sigma_gas <= 0) = 0.1*min(Sigma_gas(Sigma_gas>0));
Sigma_gas(1) = Sigma_gas(2);
R = [R;1.1*max(R)];
Sigma_stars = [Sigma_stars;0.2*Sigma_stars(end)];
Sigma_Stars2 = [Sigma_Stars2;0.2*Sigma_Stars2(end)];
Sigma_gas = [Sigma_gas;0.2*Sigma_gas(end)];
Vobs = [Vobs;Vobs(end)];
eVobs = [eVobs;eVobs(end)];
Vgas = [Vgas;Vgas(end)];
Vdisk = [Vdisk;Vdisk(end)];
Vbul = [Vbul;Vbul(end)];
Sdisk = [Sdisk;Sdisk(end)];
Sbul = [Sbul;Sbul(end)];
profile_data = [profile_data;[0 0 0 0 0 0 0 0 0 0]];
profile_data(:,1)=R;
profile_data(:,2)=Sigma_stars;
profile_data(:,3)=Sigma_gas*gasMassCorrection;
profile_data(:,4)=Vobs;
profile_data(:,5)=eVobs;
profile_data(:,6)=Vgas;
profile_data(:,7)=Vdisk;
profile_data(:,8)=Vbul;
profile_data(:,9)=Sdisk;
profile_data(:,10)=Sbul;
profile_data = [profile_data Sigma_Stars2];
end