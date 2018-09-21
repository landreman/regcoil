function curvature = fourierCurvature(thetas_refined,X_refined,Y_refined,Z_refined,mmax)

    [xmc,xms] = transformCurve(thetas_refined,X_refined,mmax);
    [ymc,yms] = transformCurve(thetas_refined,Y_refined,mmax);
    [zmc,zms] = transformCurve(thetas_refined,Z_refined,mmax);
    
    ntheta = length(X_refined);

    [thetas_transform,X_transform,dXdtheta,d2Xdtheta2] = inverseTransformCurve(xmc,xms,ntheta);
    [thetas_transform,Y_transform,dYdtheta,d2Ydtheta2] = inverseTransformCurve(ymc,yms,ntheta);
    [thetas_transform,Z_transform,dZdtheta,d2Zdtheta2] = inverseTransformCurve(zmc,zms,ntheta);

%     figure()
%     plot(thetas_refined,X_refined)
%     xlabel('theta')
%     ylabel('X')
%     
%     figure()
%     plot(thetas_transform,X_transform)
%     xlabel('theta')
%     ylabel('X transform')
    
    curvature = ((dYdtheta.*d2Zdtheta2-dZdtheta.*d2Ydtheta2).^2 + (dZdtheta.*d2Xdtheta2-dXdtheta.*d2Zdtheta2).^2 ...
        + (dXdtheta.*d2Ydtheta2-dYdtheta.*d2Xdtheta2).^2).^(0.5);
    curvature = curvature./(dXdtheta.^2 + dYdtheta.^2 + dZdtheta.^2).^(1.5);    
end