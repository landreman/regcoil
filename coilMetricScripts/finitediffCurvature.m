function curvature = finitediffCurvature(thetas_refined,X_refined,Y_refined,Z_refined)

    dtheta = thetas_refined(2)-thetas_refined(1);
    ntheta = length(thetas_refined);
    index = 1:ntheta;
    index_left = circshift(index,1,2);
    index_right = circshift(index,-1,2);
    
    dXdtheta = zeros(ntheta,1);
    dYdtheta = zeros(ntheta,1);
    dZdtheta = zeros(ntheta,1);
    d2Xdtheta2 = zeros(ntheta,1);
    d2Ydtheta2 = zeros(ntheta,1);
    d2Zdtheta2 = zeros(ntheta,1);
    for itheta=1:ntheta
       itheta_left = index_left(itheta);
       itheta_right = index_right(itheta);
       dXdtheta(itheta) = (X_refined(itheta_right)-X_refined(itheta_left))/(2*dtheta);
       dYdtheta(itheta) = (Y_refined(itheta_right)-Y_refined(itheta_left))/(2*dtheta);
       dZdtheta(itheta) = (Z_refined(itheta_right)-Z_refined(itheta_left))/(2*dtheta);
       d2Xdtheta2(itheta) = (X_refined(itheta_right)-2*X_refined(itheta)+X_refined(itheta_left))/(dtheta^2);
       d2Ydtheta2(itheta) = (Y_refined(itheta_right)-2*Y_refined(itheta)+Y_refined(itheta_left))/(dtheta^2);
       d2Zdtheta2(itheta) = (Z_refined(itheta_right)-2*Z_refined(itheta)+Z_refined(itheta_left))/(dtheta^2);
    end
    
%     figure()
%     plot(thetas_refined,dXdtheta)
%     xlabel('thetas_refined');
%     ylabel('dXdtheta');
%     
%     figure()
%     plot(thetas_refined,X_refined)
%     xlabel('thetas_refined');
%     ylabel('X_refined');
    
    curvature = ((dYdtheta.*d2Zdtheta2-dZdtheta.*d2Ydtheta2).^2 + (dZdtheta.*d2Xdtheta2-dXdtheta.*d2Zdtheta2).^2 ...
        + (dXdtheta.*d2Ydtheta2-dYdtheta.*d2Xdtheta2).^2).^(0.5);
    curvature = curvature./(dXdtheta.^2 + dYdtheta.^2 + dZdtheta.^2).^(1.5);    
end