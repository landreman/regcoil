function [thetas,X,dXdtheta,d2Xdtheta2] = inverseTransformCurve(rmc,rms,ntheta)
    mmax = (length(rmc)-1);
    ms = 0:mmax;

    thetas = linspace(0,2*pi,ntheta+1);
    thetas(end) = [];
    
    X = zeros(ntheta,1);
    dXdtheta = zeros(ntheta,1);
    d2Xdtheta2 = zeros(ntheta,1);
    nmodes = length(rmc);
    for imn = 1:nmodes
        for itheta = 1:ntheta
           cosangle = cos(ms(imn)*thetas(itheta));
           sinangle = sin(ms(imn)*thetas(itheta));
           X(itheta) = X(itheta) + rmc(imn)*cosangle;
           dXdtheta(itheta) = dXdtheta(itheta) - ms(imn)*rmc(imn)*sinangle;
           d2Xdtheta2(itheta) = d2Xdtheta2(itheta) - ms(imn)^2*rmc(imn)*cosangle;
           if (ms(imn) ~= 0)
                X(itheta) = X(itheta) + rms(imn)*sinangle;
                dXdtheta(itheta) = dXdtheta(itheta) + ms(imn)*rms(imn)*cosangle;
                d2Xdtheta2(itheta) = d2Xdtheta2(itheta) - ms(imn)^2*rms(imn)*sinangle;
           end
       end
    end

end