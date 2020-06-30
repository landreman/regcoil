%% This is a script for computing coil metrics from a regcoil
%% output file. There are several options at the beginning of the file
%% for how to parameterize the filamentary curves and how to compute 
%% the curvarure. 

clear
clc

regcoil_out_filename = '/Users/elizabethpaul/Documents/Research/August_2019/20190812_summer_school_lab/regcoil/nfp_3/regcoil_out.rotating_ellipse.nc';
coilsPerHalfPeriod = 4;

%% Choice for parameterization angle of curve
angleChoice = 1; % equal arclength
%angleChoice = 2; % thetas

thetaShift = 400;

%% Number of gridpoints for interpolation
ns_refined = 400;

%% Resolution for FT
mmax = 15;

%% How to compute curvature
whichCurvature = 1; % FT
% whichCurvature = 2; % finite diff
    
% Read regcoil_out file:
filename = regcoil_out_filename;

rmnc_coil = ncread(filename,'rmnc_coil');
zmns_coil = ncread(filename,'zmns_coil');
xm_coil = double(ncread(filename,'xm_coil'));
xn_coil = double(ncread(filename,'xn_coil'));

nfp = double(ncread(filename,'nfp'));
net_poloidal_current_Amperes = ncread(filename,'net_poloidal_current_Amperes');
theta = ncread(filename,'theta_coil');
theta = circshift(theta,thetaShift);
for itheta=1:(length(theta)-1)
    if (theta(itheta) > theta(itheta+1))
       theta(itheta+1) = theta(itheta+1)+2*pi;
    end
end

nzeta = double(ncread(filename,'nzeta_coil'));
nzetal=nzeta*nfp;
zetal = linspace(0,2*pi,nzetal+1);
zetal(end)=[];
[zetal_2D, theta_2D] = meshgrid(zetal,theta);
potential0 = ncread(filename,'current_potential');
potential1 = potential0(:,:,end);
potential1 = circshift(potential1,thetaShift,1);

potential = kron(ones(1,nfp),potential1) + kron(((1:nfp)-1)*net_poloidal_current_Amperes/nfp,ones(numel(theta),nzeta));
potential = potential / net_poloidal_current_Amperes * nfp;

contours = linspace(0,nfp,1+coilsPerHalfPeriod*2*nfp);
contours(end)= [];
dc = contours(2)-contours(1);
contours = contours + 0.5*dc;

contours_theta = cell(coilsPerHalfPeriod,1);
contours_zeta = cell(coilsPerHalfPeriod,1);
contours_x = cell(coilsPerHalfPeriod,1);
contours_y = cell(coilsPerHalfPeriod,1);
contours_z = cell(coilsPerHalfPeriod,1);
contours_dxdtheta = cell(coilsPerHalfPeriod,1);
contours_dydtheta = cell(coilsPerHalfPeriod,1);
contours_dzdtheta = cell(coilsPerHalfPeriod,1);
contours_dxdzeta = cell(coilsPerHalfPeriod,1);
contours_dydzeta = cell(coilsPerHalfPeriod,1);
contours_dzdzeta = cell(coilsPerHalfPeriod,1);
coils_x = cell(coilsPerHalfPeriod,1);
coils_y = cell(coilsPerHalfPeriod,1);
coils_z = cell(coilsPerHalfPeriod,1);
for j=1:coilsPerHalfPeriod
    this_contour = contours(j+2*coilsPerHalfPeriod);
    C = contourc(zetal,theta,potential,[this_contour,this_contour]);
    N = C(2,1);
    if N ~= size(C,2)-1
        fprintf('It appears there are multiple disconnected contours. This program presently cannot handle this.\n')
        size(C);
    end
    this_zeta = C(1,2:end)';
    this_theta = C(2,2:end)';
    contours_zeta{j} = [this_zeta; this_zeta(1)];
    contours_theta{j}  = [this_theta; this_theta(1)];
    contours_x{j} = zeros(size(contours_theta{j}));
    contours_y{j} = zeros(size(contours_theta{j}));
    contours_z{j} = zeros(size(contours_theta{j}));
    contours_dxdtheta{j} = zeros(size(contours_theta{j}));
    contours_dydtheta{j} = zeros(size(contours_theta{j}));
    contours_dzdtheta{j} = zeros(size(contours_theta{j}));
    contours_dxdzeta{j} = zeros(size(contours_theta{j}));
    contours_dydzeta{j} = zeros(size(contours_theta{j}));
    contours_dzdzeta{j} = zeros(size(contours_theta{j}));
end

x = zeros(size(theta_2D));
y = zeros(size(theta_2D));
z = zeros(size(theta_2D));
dxdtheta = zeros(size(theta_2D));
dydtheta = zeros(size(theta_2D));
dzdtheta = zeros(size(theta_2D));
dxdzeta = zeros(size(theta_2D));
dydzeta = zeros(size(theta_2D));
dzdzeta = zeros(size(theta_2D));

mnmax_coil = length(xm_coil);
for i = 1:mnmax_coil
    angle = xm_coil(i)*theta_2D - xn_coil(i)*zetal_2D;
    angle2 = zetal_2D; 

    x = x + rmnc_coil(i)*cos(angle).*cos(angle2);
    y = y + rmnc_coil(i)*cos(angle).*sin(angle2);
    z = z + zmns_coil(i)*sin(angle);

    for j=1:coilsPerHalfPeriod

        angle = xm_coil(i)*contours_theta{j} - xn_coil(i)*contours_zeta{j};
        angle2 = contours_zeta{j};

        contours_x{j} = contours_x{j} + rmnc_coil(i)*cos(angle).*cos(angle2);
        contours_y{j} = contours_y{j} + rmnc_coil(i)*cos(angle).*sin(angle2);
        contours_z{j} = contours_z{j} + zmns_coil(i)*sin(angle);

        contours_dxdtheta{j} = contours_dxdtheta{j} - xm_coil(i)*rmnc_coil(i)*sin(angle).*cos(angle2);
        contours_dydtheta{j} = contours_dydtheta{j} - xm_coil(i)*rmnc_coil(i)*sin(angle).*sin(angle2);
        contours_dzdtheta{j} = contours_dzdtheta{j} + xm_coil(i)*zmns_coil(i)*cos(angle);

        contours_dxdzeta{j} = contours_dxdzeta{j} + xn_coil(i)*rmnc_coil(i)*sin(angle).*cos(angle2) ...
            - rmnc_coil(i)*cos(angle).*sin(angle2);
        contours_dydzeta{j} = contours_dydzeta{j} + xn_coil(i)*rmnc_coil(i)*sin(angle).*sin(angle2) ...
            + rmnc_coil(i)*cos(angle).*cos(angle2);
        contours_dzdzeta{j} = contours_dzdzeta{j} - xn_coil(i)*zmns_coil(i)*cos(angle);

    end
end

extent = zeros(coilsPerHalfPeriod,1);
coilLength = zeros(coilsPerHalfPeriod,1);
length_fdiff = zeros(coilsPerHalfPeriod,1);
mean_curv_s = zeros(coilsPerHalfPeriod,1);
max_curv_s = zeros(coilsPerHalfPeriod,1);
mean_curv_theta = zeros(coilsPerHalfPeriod,1);
max_curv_theta = zeros(coilsPerHalfPeriod,1);

for j=1:coilsPerHalfPeriod
    
    X = contours_x{j}(:,1);
    Y = contours_y{j}(:,1);
    Z = contours_z{j}(:,1);
    thetas = contours_theta{j};
    % Remove repeated entries
    X(end) = [];
    Y(end) = [];
    Z(end) = [];
    thetas(end) = [];
    
    coilLength(j) = computeLength(X,Y,Z);
    extent(j) = computeExtent(X,Y);
    
    if angleChoice == 1
        [thetas_refined,X_refined,Y_refined,Z_refined] = splineArclength(X,Y,Z,ns_refined);
    end
    if angleChoice == 2
        [thetas_refined,X_refined,Y_refined,Z_refined] = splineThetas(thetas,X,Y,Z,ns_refined);
    end
    
    if whichCurvature == 1 % FT
        [xmc,xms] = transformCurve(thetas_refined,X_refined,mmax);
        [ymc,yms] = transformCurve(thetas_refined,Y_refined,mmax);
        [zmc,zms] = transformCurve(thetas_refined,Z_refined,mmax);

        [thetas_transform,X_transform,dXdtheta,d2Xdtheta2] = inverseTransformCurve(xmc,xms,ns_refined);
        [thetas_transform,Y_transform,dYdtheta,d2Ydtheta2] = inverseTransformCurve(ymc,yms,ns_refined);
        [thetas_transform,Z_transform,dZdtheta,d2Zdtheta2] = inverseTransformCurve(zmc,zms,ns_refined);

        curvature = fourierCurvature(thetas_refined,X_refined,Y_refined,Z_refined,mmax);
%         figure()
%         plot(thetas_refined,X_refined);
%         hold on
%         plot(thetas_transform,X_transform);
%         xlabel('theta')
%         legend('X refined','X transform');
%         title('arclength');
    end
    
    if whichCurvature==2 % Finite diff
        curvature = finitediffCurvature(thetas_refined,X_refined,Y_refined,Z_refined);
    end
    
    mean_curv_s(j) = mean(curvature);
    max_curv_s(j) = max(curvature);
end

% Find closest coil-coil approach
min_dist = zeros(coilsPerHalfPeriod,1);

for i=1:coilsPerHalfPeriod
    min_dist(i) = 100;
    for j=1:coilsPerHalfPeriod
        if (i ~= j)
            Xj = contours_x{j}(:,1);
            Yj = contours_y{j}(:,1);
            Zj = contours_z{j}(:,1);
            Xi = contours_x{i}(:,1);
            Yi = contours_y{i}(:,1);
            Zi = contours_z{i}(:,1);
            ntj = length(Xj);
            nti = length(Xi);
            for iv=1:nti
                for jv=1:ntj
                    X_diff = Xi(iv)-Xj(jv);
                    Y_diff = Yi(iv)-Yj(jv);
                    Z_diff = Zi(iv)-Zj(jv);
                    dist = sqrt(X_diff^2+Y_diff^2+Z_diff^2);
                    if (dist < min_dist(i))
                       min_dist(i) = dist; 
                    end
                end
            end
        end
    end
end
fprintf('min. coil-coil dist: %d [m]\n', min(min_dist));
fprintf('mean toroidal extent: %d [rad]\n', mean(extent));
fprintf('max toroidal extent: %d [rad]\n', max(extent));
fprintf('mean length: %d [m]\n', mean(coilLength));
fprintf('max length: %d [m]\n', max(coilLength));
if (angleChoice==1)
    fprintf('mean_curvature (arclength): %d [m^{-1}]\n', mean(mean_curv_s));
    fprintf('max curvature (arclength): %d [m^{-1}]\n', max(max_curv_s));
end
if (angleChoice==2)
    fprintf('mean_curvature (theta): %d [m^{-1}]\n', mean(mean_curv_theta));
    fprintf('max curvature (theta): %d [m^{-1}]\n', max(max_curv_theta));
end