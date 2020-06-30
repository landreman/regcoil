function [thetas_refined,X_refined,Y_refined,Z_refined] = splineThetas(thetas,X,Y,Z,ns_refined)
    
    [thetas,I] = sort(thetas);
    X = X(I);
    Y = Y(I);
    Z = Z(I);
    
    % Uniform theta grid 
    thetas_refined = linspace(thetas(1),thetas(end),ns_refined);
    
    figure()
    plot(thetas,X)
    title('spline thetas')
    
    % Interpolate on uniform grd
    X_interp = spline(thetas,X);
    Y_interp = spline(thetas,Y);
    Z_interp = spline(thetas,Z);
    X_refined = ppval(X_interp,thetas_refined);
    Y_refined = ppval(Y_interp,thetas_refined);
    Z_refined = ppval(Z_interp,thetas_refined);
    
    % Sort monotnically
%     [thetas_refined,I] = sort(thetas_refined);
%     X_refined = X_refined(I);
%     Y_refined = Y_refined(I);
%     Z_refined = Z_refined(I);
%     
%     figure()
%     plot(thetas,X)
%     hold on
%     plot(thetas_refined,X_refined)
%     xlabel('theta')
%     ylabel('X')
%     legend('X','X refined')
%     title('theta')
    
end