function [curvature_fdiff,curvature_spline] = splineCurvatureArclength(X,Y,Z,ns_refined)
    ns = length(X);
    % Construct constant arclength parameter
    s = zeros(ns,1);
    ds = zeros(ns,1);
    for i = 2:ns
        i_last = i - 1;
        dX = (X(i)-X(i_last));
        dY = (Y(i)-Y(i_last));
        dZ = (Z(i)-Z(i_last));
        ds(i) = sqrt(dX^2 + dY^2 + dZ^2);
        s(i) = s(i-1) + sqrt(dX^2 + dY^2 + dZ^2);
    end
    
%     figure(800)
%     plot(s,X);
%     hold on
%     plot(s,Y);
%     plot(s,Z);
%     xlabel('s');
%     legend('X','Y','Z');
  
    X_interp = spline(s,X);
    Y_interp = spline(s,Y);
    Z_interp = spline(s,Z);
    % Uniform s grid 
    s_refined = linspace(0,s(end),ns_refined);
    ds = s_refined(2)-s_refined(1);
    
    % Spline curvature
    dXds = fnder(X_interp,1);
    d2Xds2 = fnder(X_interp,2);
    dYds = fnder(Y_interp,1);
    d2Yds2 = fnder(Y_interp,2);
    dZds = fnder(Z_interp,1);
    d2Zds2 = fnder(Z_interp,2);
    dXds = ppval(dXds,s_refined);
    dYds = ppval(dYds,s_refined);
    dZds = ppval(dZds,s_refined);
    d2Xds2 = ppval(d2Xds2,s_refined);
    d2Yds2 = ppval(d2Yds2,s_refined);
    d2Zds2 = ppval(d2Zds2,s_refined);
    curvature_spline = ((dYds.*d2Zds2-dZds.*d2Yds2).^2 + (dZds.*d2Xds2-dXds.*d2Zds2).^2 ...
        + (dXds.*d2Yds2-dYds.*d2Xds2).^2).^(0.5);
    curvature_spline = curvature_spline./(dXds.^2 + dYds.^2 + dZds.^2).^(1.5);
    
    % Finite diff curvature
    X = ppval(X_interp,s_refined);
    Y = ppval(Y_interp,s_refined);
    Z = ppval(Z_interp,s_refined);
    
%     figure(801)
%     plot(s_refined,X);
%     hold on
%     plot(s_refined,Y);
%     plot(s_refined,Z);
%     xlabel('s');
%     legend('X refined','Y refined','Z refined');
    
    dXds = zeros(ns_refined,1);
    dYds = zeros(ns_refined,1);
    dZds = zeros(ns_refined,1);
    d2Xds2 = zeros(ns_refined,1);
    d2Yds2 = zeros(ns_refined,1);
    d2Zds2 = zeros(ns_refined,1);
    for i=1:ns_refined
       if (i == ns_refined)
           dXds(i) = (X(i)-X(i-1))/ds;
           dYds(i) = (Y(i)-Y(i-1))/ds;
           dZds(i) = (Z(i)-Z(i-1))/ds;
       elseif (i == 1)
           dXds(i) = (X(i+1)-X(i))/ds;
           dYds(i) = (Y(i+1)-Y(i))/ds;
           dZds(i) = (Z(i+1)-Z(i))/ds;
       else
           i_next = i+1;
           i_last = i-1;
           dXds(i) = (X(i_next)-X(i_last))/(2*ds);
           dYds(i) = (Y(i_next)-Y(i_last))/(2*ds);
           dZds(i) = (Z(i_next)-Z(i_last))/(2*ds);
       end
    end
    for i=1:ns_refined
        if (i == ns_refined)
            d2Xds2(i) = (dXds(i)-dXds(i-1))/ds;
            d2Yds2(i) = (dYds(i)-dYds(i-1))/ds;
            d2Zds2(i) = (dZds(i)-dZds(i-1))/ds;
        elseif (i == 1)
            d2Xds2(i) = (dXds(i+1)-dXds(i))/ds;
            d2Yds2(i) = (dYds(i+1)-dYds(i))/ds;
            d2Zds2(i) = (dZds(i+1)-dZds(i))/ds;
        else
            i_next = i+1;
            i_last = i-1;
            d2Xds2(i) = (X(i_next)-2*X(i)+X(i_last))/(ds^2);
            d2Yds2(i) = (Y(i_next)-2*Y(i)+Y(i_last))/(ds^2);
            d2Zds2(i) = (Z(i_next)-2*Z(i)+Z(i_last))/(ds^2);
        end
    end
    curvature_fdiff = ((dYds.*d2Zds2-dZds.*d2Yds2).^2 + (dZds.*d2Xds2-dXds.*d2Zds2).^2 ...
        + (dXds.*d2Yds2-dYds.*d2Xds2).^2).^(0.5);
    curvature_fdiff = curvature_fdiff./(dXds.^2 + dYds.^2 + dZds.^2).^(1.5);

end