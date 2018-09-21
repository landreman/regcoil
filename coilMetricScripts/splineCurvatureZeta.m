
function [curvature_fdiff, curvature_spline] = splineCurvatureZeta(X,Y,Z,zetas,ns_refined)
    % Check that zetas is monotonic
    ns = length(X);
    if (zetas(2) > zetas(1)) % increasing
        for i = 3:ns
            if (zetas(i) < zetas(i-1))
                zetas(i) = zetas(i) + 2*pi;
            end
        end
    else % decreasing
        for i = 3:ns
            if (zetas(i) > zetas(i-1))
                zetas(i) = zetas(i) - 2*pi;
            end
        end
    end
    
    ns = length(X);
    % Construct constant arclength parameter
    s = zeros(ns,1);
    for i = 2:ns
        i_last = i - 1;
        dX = (X(i)-X(i_last));
        dY = (Y(i)-Y(i_last));
        dZ = (Z(i)-Z(i_last));
        s(i) = s(i-1) + sqrt(dX^2 + dY^2 + dZ^2);
    end
    
    figure(100)
    plot(zetas,s)
    hold on
    xlabel('zeta')
    ylabel('s')
    
    % Uniform theta grid 
    zetas_refined = linspace(0,2*pi,ns_refined+1);
    zetas_refined(end) = [];

    X_interp = spline(zetas,X);
    Y_interp = spline(zetas,Y);
    Z_interp = spline(zetas,Z);
    
    dXds = fnder(X_interp,1);
    d2Xds2 = fnder(X_interp,2);
    dYds = fnder(Y_interp,1);
    d2Yds2 = fnder(Y_interp,2);
    dZds = fnder(Z_interp,1);
    d2Zds2 = fnder(Z_interp,2);
    dXds = ppval(dXds,thetas_refined);
    dYds = ppval(dYds,thetas_refined);
    dZds = ppval(dZds,thetas_refined);
    d2Xds2 = ppval(d2Xds2,thetas_refined);
    d2Yds2 = ppval(d2Yds2,thetas_refined);
    d2Zds2 = ppval(d2Zds2,thetas_refined);
    curvature_spline = ((dYds.*d2Zds2-dZds.*d2Yds2).^2 + (dZds.*d2Xds2-dXds.*d2Zds2).^2 ...
        + (dXds.*d2Yds2-dYds.*d2Xds2).^2).^(0.5);
    curvature_spline = curvature_spline./(dXds.^2 + dYds.^2 + dZds.^2).^(1.5);
    
    X = ppval(X_interp,zetas_refined);
    Y = ppval(Y_interp,zetas_refined);
    Z = ppval(Z_interp,zetas_refined);
    dXds = zeros(ns_refined,1);
    dYds = zeros(ns_refined,1);
    dZds = zeros(ns_refined,1);
    d2Xds2 = zeros(ns_refined,1);
    d2Yds2 = zeros(ns_refined,1);
    d2Zds2 = zeros(ns_refined,1);
    for i=1:ns_refined
       if (i == ns_refined)
           ds = zetas_refined(2)-zetas_refined(1);
           dXds(i) = (X(i)-X(i-1))/ds;
           dYds(i) = (Y(i)-Y(i-1))/ds;
           dZds(i) = (Z(i)-Z(i-1))/ds;
       elseif (i == 1)
           ds = zetas_refined(2)-zetas_refined(1);
           dXds(i) = (X(i+1)-X(i))/ds;
           dYds(i) = (Y(i+1)-Y(i))/ds;
           dZds(i) = (Z(i+1)-Z(i))/ds;
       else
           i_next = i+1;
           i_last = i-1;
           ds = zetas_refined(3)-zetas_refined(1);
           dXds(i) = (X(i_next)-X(i_last))/ds;
           dYds(i) = (Y(i_next)-Y(i_last))/ds;
           dZds(i) = (Z(i_next)-Z(i_last))/ds;
       end
    end
    for i=1:ns_refined
        if (i == ns_refined)
            ds = zetas_refined(2)-zetas_refined(1);
            d2Xds2(i) = (X(i)-2*X(i-1)+X(i-2))/ds^2;
            d2Yds2(i) = (Y(i)-2*Y(i-1)+Y(i-2))/ds;
            d2Zds2(i) = (Z(i)-2*Z(i-1)+Z(i-2))/ds;
        elseif (i == 1)
            ds = zetas_refined(2)-zetas_refined(1);
            d2Xds2(i) = (X(i+2)-2*X(i+1)+X(i))/ds;
            d2Yds2(i) = (Y(i+1)-2*Y(i+1)+Y(i))/ds;
            d2Zds2(i) = (Z(i+1)-2*Z(i+1)+Z(i))/ds;
        else
            i_next = i+1;
            i_last = i-1;
            ds = zetas_refined(i)-zetas_refined(i_last);
            d2Xds2(i) = (X(i_next)-2*X(i)+X(i_last))/(ds^2);
            d2Yds2(i) = (Y(i_next)-2*Y(i)+Y(i_last))/(ds^2);
            d2Zds2(i) = (Z(i_next)-2*Z(i)+Z(i_last))/(ds^2);
        end
    end

    curvature_fdiff = ((dYds.*d2Zds2-dZds.*d2Yds2).^2 + (dZds.*d2Xds2-dXds.*d2Zds2).^2 ...
        + (dXds.*d2Yds2-dYds.*d2Xds2).^2).^(0.5);
    curvature_fdiff = curvature_fdiff./(dXds.^2 + dYds.^2 + dZds.^2).^(1.5);
end