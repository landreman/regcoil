function [thetas_refined,X_refined,Y_refined,Z_refined] = splineArclength(X,Y,Z,ns_refined)
    ns = length(X);
    s = zeros(ns,1);
    for i=2:ns
       dX = X(i)-X(i-1);
       dY = Y(i)-Y(i-1);
       dZ = Z(i)-Z(i-1);
       s(i) = s(i-1) + sqrt(dX^2 + dY^2 + dZ^2); 
    end
    s_tot = s(end) + sqrt((X(end)-X(1))^2+(Y(end)-Y(1))^2+(Z(end)-Z(1))^2);
    s = s/s_tot;
    thetas = s*2*pi;

    %thetas_refined = linspace(thetas(1),thetas(end),ns_refined);
    thetas_refined = linspace(0,2*pi,ns_refined+1);
    thetas_refined(end) = [];
    
    % Interpolate on uniform grid
    X_interp = spline(thetas,X);
    Y_interp = spline(thetas,Y);
    Z_interp = spline(thetas,Z);
    X_refined = ppval(X_interp,thetas_refined);
    Y_refined = ppval(Y_interp,thetas_refined);
    Z_refined = ppval(Z_interp,thetas_refined);
    
%     figure()
%     plot(thetas,X)
%     hold on
%     plot(thetas_refined,X_refined)
%     xlabel('theta')
%     legend('X','X refined')
%     title('arclength')

end