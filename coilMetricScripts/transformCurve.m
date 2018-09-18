
function [rmc,rms] = transformCurve(thetas,X,mmax)
    ms = 0:mmax;
    nmodes = length(ms);
    rmc = zeros(nmodes,1);
    rms = zeros(nmodes,1);
    
    dtheta = thetas(2)-thetas(1);
    
    ntheta = length(thetas);
    for imn = 1:nmodes
        for itheta=1:ntheta
            rmc(imn) = rmc(imn)+X(itheta)*cos(ms(imn)*thetas(itheta));
            if (ms(imn) ~= 0)
                rms(imn) = rms(imn)+X(itheta)*sin(ms(imn)*thetas(itheta));
            end
        end
        if (ms(imn) ~= 0)
           rmc(imn) = rmc(imn)*dtheta/pi; 
           rms(imn) = rms(imn)*dtheta/pi;
        else
           rmc(imn) = rmc(imn)*dtheta/(2*pi);
           rms(imn) = rms(imn)*dtheta/(2*pi);
        end
    end
end