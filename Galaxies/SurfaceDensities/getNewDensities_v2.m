function [Sigma_obs,Sigma_est,NSigma_S,NSigma_H,massAdjS,massAdjH,testwhere] = getNewDensities_v2(r_Kpc,Vobs,Vest,Sigma_S,Sigma_H,alpha,delay,ggain,testwhere)
%getNewDensities Estimates the new Star and Gas densities from estimated
%velocites
%   Detailed explanation goes here

    G = 4.30091e-3;  % pc Msun^{-1} (km s^{-1})^2

    if (delay > 0)
        delay = sum(r_Kpc<delay);
        if  (delay<2)
            delay = 2;
        end
        if  (delay>20)
            delay = 20;
        end
    end
    r_pc = 1000*r_Kpc;

    V2r = Vobs .^ 2;
    V2r = V2r.* r_pc;
    dV2r_dr = smoothAbsGrad(V2r,r_pc);
    
    rden = r_pc;
    rden(1) = 0.5*(r_pc(1) + r_pc(2));

    Sigma_obs = dV2r_dr ./ (2 * pi * G * rden);  % Msun/pc^2

    V2r = Vest .^ 2;
    V2r = V2r.* r_pc;
    dV2r_dr = smoothAbsGrad(V2r,r_pc);
    
    Sigma_est = dV2r_dr ./ (2 * pi * G * rden);  % Msun/pc^2

    Sigma_obs(1) = 0.5*(Sigma_obs(1) + Sigma_obs(2));
    Sigma_est(1) = 0.5*(Sigma_est(1) + Sigma_est(2));
    npts = length(Sigma_H);
    fistpoint = floor(1*npts/20);
    lastpoint = floor(19*npts/20);
    indx = 1:npts;

    if delay==0
        testwhere = (2*Sigma_H)>Sigma_S;
        if (sum(testwhere)>10)
            odata=Sigma_H;
            testwhere = indx(testwhere);
            testwhere = min(testwhere):max(testwhere);
            Sigma_H(testwhere) = 0.10*Sigma_H(testwhere) + 0.60*Sigma_S(testwhere) + 0.30*mean(Sigma_S(testwhere));
            Sigma_H(odata < Sigma_H) = odata(odata < Sigma_H);
            Sigma_H(Sigma_H<=0) = min(Sigma_H(Sigma_H>0));
            Sigma_H = log(Sigma_H);
            Sigma_H(testwhere) = smooth(Sigma_H(testwhere),0.50);
            Sigma_H = exp(Sigma_H);
        else
            testwhere = npts;
        end
        Sigma_H = log(Sigma_H);
        odata=Sigma_H;
        Sigma_H(fistpoint:lastpoint) = smooth(Sigma_H(fistpoint:lastpoint),0.25);
        Sigma_H = (Sigma_H + odata)/2;
        Sigma_H = exp(Sigma_H);
    end


    PSigma_error = (Sigma_est - Sigma_obs)./(1.0e-6 + Sigma_est);
    PSigma_error(2) = (PSigma_error(2) + PSigma_error(3))/2.0;
    PSigma_error(1) = 0.1*PSigma_error(1) + 0.9*PSigma_error(2);
    PSigma_error(end-1) = (PSigma_error(end-1) + PSigma_error(end-2))/2.0;
    PSigma_error(end) = (PSigma_error(end) + PSigma_error(end-1))/2.0;

%    mean(PSigma_error)
    if delay>0
        odata=PSigma_error;
        PSigma_error(fistpoint:lastpoint) = smooth(PSigma_error(fistpoint:lastpoint),0.05);
        PSigma_error = (odata + PSigma_error)/2;

    else
        alpha = 1.0;
        wheremean = (Sigma_S>Sigma_H) & (linspace(0,npts,npts) < (npts/2))';
        PSigma_error(:) = mean(PSigma_error(wheremean));
    end
    totmass = Sigma_S + Sigma_H;
    massAdj = PSigma_error;

    starratio = Sigma_S./totmass;
    

    % Adjusting Stars
    massAdjS = massAdj;
    massAdjH = massAdj;
    NSigma_S = Sigma_S - alpha*massAdj.*totmass.*starratio;
    NSigma_H = Sigma_H - alpha*massAdj.*totmass.*(1.0-starratio);
    starratio(starratio<=0) = min(starratio(starratio>0));

    if (delay > 0)
        madv=0.75;

        % The Star gain
        delay2 = ceil(delay/2);
        massAdjS(1:(end-delay2)) = PSigma_error((delay2+1):end);
        massAdjS(end-delay:end) = smooth(massAdjS(end-delay:end),0.5);
        massAdjS = 0.5*(massAdjS + massAdj);
        odata=massAdjS;    
        massAdjS(fistpoint:lastpoint)= smooth(massAdjS(fistpoint:lastpoint),0.05);
        massAdjS = (massAdjS + odata)/2;

        massAdjS = alpha*massAdjS;

    
        % The Gas gain
        massAdjH = massAdj;
        massAdjH((delay+1):end) = PSigma_error(1:(end-delay));
        massAdjH(1:(2*delay)) = smooth(massAdjH(1:(2*delay)),0.5);
        massAdjH = 0.5*(massAdjH + massAdj);
        odata=massAdjH;
        massAdjH(fistpoint:lastpoint)= smooth(massAdjH(fistpoint:lastpoint),0.05);
        massAdjH = (massAdjH + odata)/2;
        massAdjH = alpha*massAdjH;

        o_starratio = starratio;
        gasRatio = (1.0-starratio);
        o_gasRatio = gasRatio;

        if (ggain <= 1)
            starratio(testwhere) =  starratio(testwhere)*ggain;
            gasRatio = 1.0-starratio;
        else
            gasRatio(testwhere) =  gasRatio(testwhere)/ggain;
            starratio = 1.0-gasRatio;
        end

        calgain = max(abs(massAdjS),abs(massAdjH));
        calgain(calgain<=madv) = 1.0;
        calgain(calgain>madv) = madv./calgain(calgain>madv);
        calgain = smooth(calgain,0.10);


        % Correct Stars
        massAdjS = massAdjS.*calgain.*starratio;
        massAdjS(massAdjS > madv*o_starratio) = madv*o_starratio(massAdjS > madv*o_starratio);
        massAdjS(massAdjS < -madv*o_starratio) = -madv*o_starratio(massAdjS < -madv*o_starratio);

        NSigma_S = Sigma_S - massAdjS.*totmass;
        NSigma_S(NSigma_S<=0) = min(Sigma_S(Sigma_S>0));

        % Correct Gas
        massAdjH = massAdjH.*calgain.*gasRatio;
        massAdjH(massAdjH > madv*o_gasRatio) = madv*o_gasRatio(massAdjH > madv*o_gasRatio);
        massAdjH(massAdjH < -madv*o_gasRatio) = -madv*o_gasRatio(massAdjH < -madv*o_gasRatio);

        NSigma_H = Sigma_H - massAdjH.*totmass;
        NSigma_H(NSigma_H<=0) = min(Sigma_H(Sigma_H>0));
    
    end

    NSigma_S = log(NSigma_S);
    NSigma_S(floor(3*npts/5):npts) = smooth(NSigma_S(floor(3*npts/5):npts),55,'sgolay',1);
    NSigma_S(floor(2*npts/5):floor(4*npts/5)) = smooth(NSigma_S(floor(2*npts/5):floor(4*npts/5)),0.50);
    NSigma_S(fistpoint:lastpoint) = smooth(NSigma_S(fistpoint:lastpoint),0.10);
    NSigma_S = exp(NSigma_S);

    NSigma_H = log(NSigma_H);
    NSigma_H(floor(2*npts/4):npts) = smooth(NSigma_H(floor(2*npts/4):npts),55,'sgolay',1);
    NSigma_H(floor(2*npts/5):floor(4*npts/5)) = smooth(NSigma_H(floor(2*npts/5):floor(4*npts/5)),0.30);
    NSigma_H(fistpoint:lastpoint) = smooth(NSigma_H(fistpoint:lastpoint),0.10);
    NSigma_H = exp(NSigma_H);

end