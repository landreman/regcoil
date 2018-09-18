function coilLength=computeLength(X,Y,Z)
    ns = length(X);
    
    coilLength=0;
    for i=2:ns
        dX = X(i)-X(i-1);
        dY = Y(i)-Y(i-1);
        dZ = Z(i)-Z(i-1);
        coilLength = coilLength + sqrt(dX^2+dY^2+dZ^2);
    end
    coilLength = coilLength + sqrt((X(end)-X(1))^2+(Y(end)-Y(1))^2+(Z(end)-Z(1))^2);

end