function [gacel_x,gacel_y,gacel0_x,gacel0_y,Zdist,StarProfile,GasZprofile] = galaticRingGR_vec(Sig_star,Sig_H,B_sigma,A_Sigma_star,A_Sigma_H,locIndex,distances,dVelocities,MStars,Mgas,nPoints)
    %galaticRingG computes the aceleration along the X axis of a ring of mater
    %rotating at the specified velocity. The aceleration is computed using only
    %electromagnetic atraction between neutral electric patches.
    %Sig_star = Velocity RMS of Stars
    %Sig_H = Velocity RMS of Gass
    %radius =  The ring radius (pc)
    %surfaceDensity =  the observed barionic surface density at the given radius (M_sun/pc^2)
    %ringDelta = how thick is the ring (pc)
    %deltStd = Change in Velocity across the delta radius (Km/s)
    %ringVelocity = the observed barionic mass velocity (km/s)
    %distances = the radial distances to compute thForce_M_sune force (pc)
    %dVelocities = Barionic velocities at each distances (km/s)
    %MStars = Stars density
    %Mgas = Gas density
    %nPoints = The number of integration points
    
    ringVelocity=dVelocities(locIndex);
    radius = distances(locIndex);
    ringDelta = distances(locIndex);
    if (locIndex>1)
        ringDelta = distances(locIndex)-distances(locIndex-1);
    end

    % Some Constants
    e = 1.602176634e-19; % elementary charge (C)
    M_sun = 1.988416e30;  % Sun mass (Kg)
    nparticles = 1.63e57; % Number of atoms in solar mass
    c = 299792.458; % Speed of light in vacuum (Km/s)
    E0 = 8.8541878188e-12; %vacuum electric permittivity (F⋅m−1)
    pc = 3.0857e16;   % 1 pc in meters (m)

    % Unitary Electric force between two charged particles (At 1m)
    Fe = e^2/(4*pi*E0);
    % Total unitary Aceleration of two Suns at 1 parsec
    SunAc = (Fe*(nparticles/pc)^2)/M_sun; % (m/s^2 per Sun at 1pc)


    thetapints = double(2*round(50 + 2*radius*nPoints/max(distances)));
    if thetapints>nPoints
        thetapints = nPoints;
    end

    theta = linspace(0, double(2*pi*(thetapints-1)./thetapints), thetapints);
    dteta = theta(2);
    theta = theta + 0.5*dteta;

    l_x = -ringVelocity*sin(theta)/c; % The x component of the velocity
    l_y = ringVelocity*cos(theta)/c; % The y component of the velocity    


    acel_x = zeros(size(distances));
    acel_y = zeros(size(distances));
    tacel_x = zeros(size(distances));
    tacel_y = zeros(size(distances));
    acel_x0 =  zeros(size(distances));
    acel_y0 =  zeros(size(distances));


    lambda_star = 2*pi*MStars(locIndex)*radius*ringDelta/thetapints; % Total # Suns in Stars per ring delta
    lambda_H = 2*pi*Mgas(locIndex)*radius*ringDelta/thetapints; % Total # Suns in GAS per ring delta
    smalldist = max(distances)/length(distances);
    smalldist = 0.1*min([smalldist,dteta*radius]);
    smalldist2 = smalldist^2;

    fanomaly0 = averageSamePatchAnomalityNoEscR_H(Sig_star,A_Sigma_star,B_sigma,0,0,0,0,0,0);    
    Zdist = linspace(0.0,max(distances),round(length(distances)/3));
    zwidth = Zdist(2)/2;
    Zdist = Zdist - zwidth;
    Zdist(1) = 0;

    totSmas = sum(MStars,"all");
    totSgas = sum(Mgas,"all");

    starRadius = 0.5*(distances'*MStars)/totSmas;
    gasRadius = 0.5*(distances'*Mgas)/totSgas;
    

    GasZprofile = (exp(-Zdist/(0.5*gasRadius)) + (sech(-Zdist/(gasRadius/16)).^2));
    GasZprofile(1) = 1;
    GasZprofile = GasZprofile/sum(GasZprofile);
    StarProfile = 2*(sech(-Zdist/(0.5*starRadius)).^2);
    StarProfile(1) = 1;
    StarProfile = StarProfile/sum(StarProfile);
    Zdist = Zdist(GasZprofile>0.01);
    lz = length(Zdist);

    GasZprofile = GasZprofile(1:lz);
    StarProfile = StarProfile(1:lz);
    GasZprofile = GasZprofile/sum(GasZprofile);
    StarProfile = StarProfile/sum(StarProfile);


    lambda_starZ = lambda_star.*StarProfile;
    lambda_HZ = lambda_H.*GasZprofile;
    
    for z = 1:lz
        for i = 1:length(distances)
    
            r=distances(i);
    
            r_y = dVelocities(i)/c; % rotation velocity in y direction
    
            if (i == 1)
                acel_x(i) = 0; % at center
                acel_y(i) = 0; % at center
            else
              % Integrate the net electrical anomality
                xpos = r - radius*cos(theta);
                ypos = -radius*sin(theta);
                dist2 = xpos.^2 + ypos.^2;

                distv2 = dist2;

                dist2 = dist2 + Zdist(z).^2;
                dist = sqrt(dist2);
    
                disvt = sqrt(distv2);
    
                urvect = [xpos./disvt; -ypos./disvt]';
                utvect = [ypos./disvt;xpos./disvt]';
                xvect = [xpos./disvt; ypos./disvt]';
                yvect = [-ypos./disvt;xpos./disvt]';

                dist(dist<smalldist) = smalldist;

                zcos = disvt./dist;
                if (z==1)
                    zcos = dist2./(disvt.*sqrt(zwidth^2+disvt.^2));
                end
                
                vr = (l_x'.*urvect(:,1) + l_y'.*urvect(:,2))';
                vt = (l_x'.*utvect(:,1) + l_y'.*utvect(:,2))';
                rv_r = (r_y'.*urvect(:,2))';
                rv_t = (r_y'.*utvect(:,2))';


                fanomaly_S_S = averageSamePatchAnomalityNoEscR_H(Sig_star,A_Sigma_star,B_sigma,vr,vt,0,rv_r,rv_t,0);
                fanomaly_S_G = averagePatchAnomalityNoEscR_H(Sig_star,Sig_H,A_Sigma_H,B_sigma,A_Sigma_star,vr,vt,0,rv_r,rv_t,0);   
                fanomaly_G_G = averageSamePatchAnomalityNoEscR_H(Sig_H,A_Sigma_H,B_sigma,vr,vt,0,rv_r,rv_t,0);
                fanomaly_G_S = averagePatchAnomalityNoEscR_H(Sig_H,Sig_star,A_Sigma_star,B_sigma,A_Sigma_H,vr,vt,0,rv_r,rv_t,0);

                fanomaly_S_S = reshape(fanomaly_S_S,[thetapints,3]);
                fanomaly_S_G = reshape(fanomaly_S_G,[thetapints,3]);
                fanomaly_G_G = reshape(fanomaly_G_G,[thetapints,3]);
                fanomaly_G_S = reshape(fanomaly_G_S,[thetapints,3]);

                vxanomality_S_S = fanomaly_S_S(:,1).*xvect(:,1) + fanomaly_S_S(:,2).*xvect(:,2);
                vyanomality_S_S = fanomaly_S_S(:,1).*yvect(:,1) + fanomaly_S_S(:,2).*yvect(:,2);

%                size(fanomaly_S_S)
%                size(vxanomality_S_S)
%                size(vyanomality_S_S)
    
                vxanomality_G_G = fanomaly_G_G(:,1).*xvect(:,1) + fanomaly_G_G(:,2).*xvect(:,2);
                vyanomality_G_G = fanomaly_G_G(:,1).*yvect(:,1) + fanomaly_G_G(:,2).*yvect(:,2);

                vxanomality_S_G = fanomaly_S_G(:,1).*xvect(:,1) + fanomaly_S_G(:,2).*xvect(:,2);
                vyanomality_S_G = fanomaly_S_G(:,1).*yvect(:,1) + fanomaly_S_G(:,2).*yvect(:,2);

    
                vxanomality_G_S = fanomaly_G_S(:,1).*xvect(:,1) + fanomaly_G_S(:,2).*xvect(:,2);
                vyanomality_G_S = fanomaly_G_S(:,1).*yvect(:,1) + fanomaly_G_S(:,2).*yvect(:,2);
    
                vxanomality0 = fanomaly0(1)*xvect(:,1) + fanomaly0(2)*xvect(:,2);
                vyanomality0 = fanomaly0(1)*yvect(:,1) + fanomaly0(2)*yvect(:,2);

                
                dist2(dist2<smalldist2) = smalldist2;
    
                integrand = zcos./dist2;
                

                integrand_x = integrand.*vxanomality_S_S';



                tacel_x(i) = MStars(i)*lambda_starZ(z)*trapz(integrand_x);
                integrand_x = integrand.*vxanomality_S_G';
                tacel_x(i) = tacel_x(i) + Mgas(i)*lambda_starZ(z)*trapz(integrand_x);
                integrand_x = integrand.*vxanomality_G_S';
                tacel_x(i) = tacel_x(i) + MStars(i)*lambda_HZ(z)*trapz(integrand_x);
                integrand_x = integrand.*vxanomality_G_G';
                tacel_x(i) = tacel_x(i) + Mgas(i)*lambda_HZ(z)*trapz(integrand_x);
                tacel_x(i) = tacel_x(i)/(Mgas(i) + MStars(i));
                
                acel_x(i) = acel_x(i) + tacel_x(i);

    
                integrand_y = integrand.*vyanomality_S_S';
                tacel_y(i) = MStars(i)*lambda_starZ(z)*trapz(integrand_y);
                integrand_y = integrand.*vyanomality_S_G';
                tacel_y(i) = tacel_y(i) + Mgas(i)*lambda_starZ(z)*trapz(integrand_y);
                integrand_y = integrand.*vyanomality_G_S';
                tacel_y(i) = tacel_y(i) + MStars(i)*lambda_HZ(z)*trapz(integrand_y);
                integrand_y = integrand.*vyanomality_G_G';
                tacel_y(i) = tacel_y(i) + Mgas(i)*lambda_HZ(z)*trapz(integrand_y);
                tacel_y(i) = tacel_y(i)/(Mgas(i) + MStars(i));

                acel_y(i) = acel_y(i) + tacel_y(i);

    
                integrand_x = integrand.*vxanomality0';
                acel_x0(i) = acel_x0(i) + (lambda_starZ(z)+lambda_HZ(z))*trapz(integrand_x);
                integrand_y = integrand.*vyanomality0';
                acel_y0(i) = acel_y0(i) + (lambda_starZ(z)+lambda_HZ(z))*trapz(integrand_y);
            end
        end

    end

    gacel_x = SunAc*acel_x;
    gacel_y = SunAc*acel_y;
    gacel0_x = SunAc*acel_x0;
    gacel0_y = SunAc*acel_y0;
end