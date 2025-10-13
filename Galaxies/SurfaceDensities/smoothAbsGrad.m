function DyDX = smoothAbsGrad(Y,X)
%smoothAbsGrad computes the abs gradient of velocities squared for density
%estimation
%   It estimates derivatives and then smooth the estimations.
    pts = length(Y);
    dV2r_dr = diff(Y)./diff(X);
    dV2r_dr2 = flip(diff(flip(Y))./diff(flip(X)));
    dsize = length(dV2r_dr2);
    dV2r_s = 0.5*(dV2r_dr(2:dsize)+dV2r_dr2(1:(dsize-1)));
    dV2r_dr = ([dV2r_dr(1);dV2r_s;dV2r_dr2(end)]);
    grdd = gradient(Y,X);
    DyDX = dV2r_dr;
    DyDX(2:(end-1))= 0.5*(dV2r_dr(2:(end-1)) + grdd(2:(end-1)));
    DyDX = abs(DyDX);

    DyDX(DyDX <=0 ) = min(DyDX(DyDX>0));
    DyDX = log(DyDX);
    DyDX(floor(3*pts/4):pts) = smooth(DyDX(floor(3*pts/4):pts),0.35,'lowess');
    DyDX(floor(pts/4):pts) = smooth(DyDX(floor(pts/4):pts),0.1);
    DyDX(floor(pts/20):floor(pts/3)) = smooth(DyDX(floor(pts/20):floor(pts/3)),0.05);
    DyDX(floor(pts/20):pts) = smooth(DyDX(floor(pts/20):pts),0.05);
    DyDX = smooth(DyDX,0.02);
    DyDX = exp(DyDX);

end