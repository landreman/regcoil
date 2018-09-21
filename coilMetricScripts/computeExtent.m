function extent = computeExtent(X,Y)
    phi = atan2(Y,X);
    phi_min = min(min(phi));
    phi_max = max(max(phi));
    extent = phi_max-phi_min;
end