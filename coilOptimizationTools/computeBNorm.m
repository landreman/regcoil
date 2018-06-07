%% This script assumes that REGCOIL has been run with save_level = 0
%% For a given regcoil output file, the relative normal field
%% (accounting for finite coil effects) is computed

clear all;
clc;

%% Global variables
global nfp;
global mu_0;
global coilsPerHalfPeriod;
global Bnormal_from_plasma_current;
global nx; 
global ny;
global nz;
global N;
global r_plasma;
global dtheta_plasma;
global dzeta_plasma;
global Ntheta_plasma;
global Nzeta_plasma;

%% Input parameters
regcoil_out_filename = '/Users/elizabethpaul/Documents/Research/Spring_2018/20180605_QI_REGCOIL/offset_1.6/regcoil_scan/chi2_B_2.55/regcoil_out.chi2_B_2.55.nc';
nescin_filename = '/Users/elizabethpaul/Documents/Research/Spring_2018/20180605_QI_REGCOIL/offset_1.6/regcoil_scan/chi2_B_2.55/nescin.offset';
wout_filename = '/Users/elizabethpaul/Documents/Research/Spring_2018/20180605_QI_REGCOIL/offset_1.6/regcoil_scan/chi2_B_2.55/wout_35.nc';
thetaShift = 10;
ilambda = 2;
coilsPerHalfPeriod = 6;


mu_0 = 4*pi * 1.0e-7;

%% Read regcoil_out file:
filename = regcoil_out_filename;
nfp = double(ncread(filename,'nfp'));
net_poloidal_current_Amperes = ncread(filename,'net_poloidal_current_Amperes');
normal_plasma = ncread(filename,'normal_plasma');
norm_normal_plasma = ncread(filename,'norm_normal_plasma');
r_plasma = ncread(filename,'r_plasma');
theta_plasma = ncread(filename,'theta_plasma');
zeta_plasma = ncread(filename,'zeta_plasma');
Ntheta_plasma = ncread(filename,'ntheta_plasma');
Nzeta_plasma = ncread(filename,'nzeta_plasma');
dtheta_plasma = theta_plasma(2)-theta_plasma(1);
dzeta_plasma = zeta_plasma(2)-zeta_plasma(1);

Bnormal_from_plasma_current_1 = ncread(filename,'Bnormal_from_plasma_current');
Bnormal_from_net_coil_currents_1 = ncread(filename,'Bnormal_from_net_coil_currents');
Bnormal_from_net_coil_currents = [];
for jfp=1:nfp
    Bnormal_from_plasma_current = horzcat(Bnormal_from_plasma_current,Bnormal_from_plasma_current_1);
    Bnormal_from_net_coil_currents = horzcat(Bnormal_from_net_coil_currents,Bnormal_from_net_coil_currents_1);
    N = horzcat(N,norm_normal_plasma);
end
normal_plasma_1 = reshape(normal_plasma(1,:,:),Ntheta_plasma,nfp*Nzeta_plasma);
normal_plasma_2 = reshape(normal_plasma(2,:,:),Ntheta_plasma,nfp*Nzeta_plasma);
normal_plasma_3 = reshape(normal_plasma(3,:,:),Ntheta_plasma,nfp*Nzeta_plasma);

nx = normal_plasma_1./N;
ny = normal_plasma_2./N;
nz = normal_plasma_3./N;

Ntheta_plasma = ncread(filename,'ntheta_plasma');
Nzeta_plasma = ncread(filename,'nzeta_plasma');

theta = ncread(filename,'theta_coil');
zeta = ncread(filename,'zeta_coil');
theta = circshift(theta,thetaShift);
% Make sure theta is increasing
for itheta=1:(length(theta)-1)
    if (theta(itheta) > theta(itheta+1))
       theta(itheta+1) = theta(itheta+1)+2*pi;
    end
end

nzeta = double(ncread(filename,'nzeta_coil'));
potential0 = ncread(filename,'current_potential');
potential1 = potential0(:,:,ilambda);
potential1 = circshift(potential1,thetaShift,1);
potential = potential1 / net_poloidal_current_Amperes * nfp;
potential_3 = transpose([potential-1,potential,potential+1]);
zeta_3 = cat(1,zeta-2*pi/nfp,zeta,zeta+2*pi/nfp);

contours = linspace(0,1,coilsPerHalfPeriod*2+1);
contours(end) = [];
dc = contours(2) - contours(1);
contours = contours + dc/2;

contour_zeta=cell(2*coilsPerHalfPeriod*nfp,1);
contour_theta=cell(2*coilsPerHalfPeriod*nfp,1);
numCoils = 2*coilsPerHalfPeriod;
for j=1:numCoils
    C = contourc(zeta_3,theta,transpose(potential_3),[contours(j),contours(j)]);
    ncont = C(2,1);
    if ncont ~= size(C,2)-1
        fprintf('It appears there are multiple disconnected contours. This program presently cannot handle this.\n')
        size(C);
    end
    this_zeta = C(1,2:end)';
    this_theta = C(2,2:end)';
    for jfp=1:nfp
        d = 2*pi*(jfp-1)/nfp;
        contour_zeta{numCoils*(jfp-1)+j} = this_zeta+d;
        contour_theta{numCoils*(jfp-1)+j} = this_theta;
    end
end

%% Read surface from nescin file:
filename = nescin_filename;
fid = fopen(filename,'r');
search_string = '------ Current Surface';
while true
line = fgetl(fid);
    if strncmp(line,search_string,numel(search_string))
        break
    end
end
line = fgetl(fid); %Number of fourier modes in table
line = fgetl(fid);
mnmax_nescin = sscanf(line,'%d');
line = fgetl(fid); %Table of fourier coefficients
line = fgetl(fid); %m,n,crc2,czs2,crs2,czc2
xm_nescin = zeros(mnmax_nescin,1);
xn_nescin = zeros(mnmax_nescin,1);
rmnc_nescin = zeros(mnmax_nescin,1);
zmns_nescin = zeros(mnmax_nescin,1);
for i=1:mnmax_nescin
    line = fgetl(fid);
    data = sscanf(line,'%d %d %g %g %g %g %g %g');
    xm_nescin(i) = data(1);
    xn_nescin(i) = data(2);
    rmnc_nescin(i) = data(3);
    zmns_nescin(i) = data(4);
end
fclose(fid);

contour_R = cell(nfp*numCoils,1);
contour_Z = cell(nfp*numCoils,1);
contour_X = cell(nfp*numCoils,1);
contour_Y = cell(nfp*numCoils,1);
for j =1:nfp*numCoils
   contour_R{j} = zeros(size(contour_theta{j}));
   contour_Z{j} = zeros(size(contour_theta{j}));
   contour_X{j} = zeros(size(contour_theta{j}));
   contour_Y{j} = zeros(size(contour_theta{j}));
end

for i = 1:mnmax_nescin
    for j=1:numCoils*nfp
        angle = xm_nescin(i)*contour_theta{j} + xn_nescin(i)*contour_zeta{j}*nfp;
        angle2 = contour_zeta{j};
        contour_R{j} = contour_R{j} + rmnc_nescin(i)*cos(angle);
        contour_Z{j} = contour_Z{j} + zmns_nescin(i)*sin(angle);
    end
end

for j=1:numCoils*nfp
   contour_X{j} = contour_R{j}.*cos(contour_zeta{j});
   contour_Y{j} = contour_R{j}.*sin(contour_zeta{j});
end

numCoils = coilsPerHalfPeriod*2*nfp;
current = net_poloidal_current_Amperes / numCoils;
chi2B(contour_X,contour_Y,contour_Z,current)

function [Bx,By,Bz] = biotSavart(x_coil,y_coil,z_coil,current)
    global mu_0;
    global r_plasma;
    global Ntheta_plasma;
    global Nzeta_plasma;
    global nfp;
    
    x_plasma = zeros(Ntheta_plasma,Nzeta_plasma*nfp);
    y_plasma = zeros(Ntheta_plasma,Nzeta_plasma*nfp);
    z_plasma = zeros(Ntheta_plasma,Nzeta_plasma*nfp);
    x_plasma(:,:) = r_plasma(1,:,:);
    y_plasma(:,:) = r_plasma(2,:,:);
    z_plasma(:,:) = r_plasma(3,:,:);
    
    figure(4)
    surf(x_plasma,y_plasma,z_plasma)
    hold on
    plot3(x_coil,y_coil,z_coil)

    Nt = length(x_coil);
    Bx = 0*x_plasma;
    By = 0*x_plasma;
    Bz = 0*x_plasma;
    
    unit = ones(size(x_plasma));
    x_coil_prev = x_coil(1);
    y_coil_prev = y_coil(1);
    z_coil_prev = z_coil(1);
    for it=2:Nt
        rx = x_plasma - x_coil(it)*unit;
        ry = y_plasma - y_coil(it)*unit;
        rz = z_plasma - z_coil(it)*unit;
        dx_coil = x_coil(it) - x_coil_prev;
        dy_coil = y_coil(it) - y_coil_prev;
        dz_coil = z_coil(it) - z_coil_prev;
        r = ((x_coil(it)*unit-x_plasma).^2+(y_coil(it)*unit-y_plasma).^2+(z_coil(it)*unit-z_plasma).^2).^(.5);
        Bx = Bx + (dy_coil*rz-dz_coil*ry)./(r.^3);
        By = By + (dz_coil*rx-dx_coil*rz)./(r.^3);
        Bz = Bz + (dx_coil*ry-dy_coil*rx)./(r.^3);
        
        x_coil_prev = x_coil(it);
        y_coil_prev = y_coil(it);
        z_coil_prev = z_coil(it);
    end
    % Joint first and last points
    dx_coil = x_coil(1) - x_coil_prev;
    dy_coil = y_coil(1) - y_coil_prev;
    dz_coil = z_coil(1) - z_coil_prev;
    r = ((x_coil(1)*unit-x_plasma).^2+(y_coil(1)*unit-y_plasma).^2+(z_coil(1)*unit-z_plasma).^2).^(.5);
    Bx = Bx + (dy_coil*rz-dz_coil*ry)./(r.^3);
    By = By + (dz_coil*rx-dx_coil*rz)./(r.^3);
    Bz = Bz + (dx_coil*ry-dy_coil*rx)./(r.^3);
    Bx = Bx*current*mu_0/(4*pi);
    By = By*current*mu_0/(4*pi);
    Bz = Bz*current*mu_0/(4*pi);
end

% Compute integrated field error on plasma surface due to theta_coils,
% zeta_coils with current_coils
function result = chi2B(contours_x,contours_y,contours_z,current)
    global nfp;
    global coilsPerHalfPeriod;
    global Bnormal_from_plasma_current;
    global nx;
    global ny;
    global nz;
    
    Bx = 0;
    By = 0;
    Bz = 0;
    for icoil=1:2*coilsPerHalfPeriod*nfp
       [this_Bx,this_By,this_Bz] = biotSavart(contours_x{icoil},contours_y{icoil},contours_z{icoil},current);
       Bx = Bx + this_Bx;
       By = By + this_By;
       Bz = Bz + this_Bz;
    end
    B_normal = (Bx.*nx+By.*ny+Bz.*nz);
    B_mag = sqrt(Bx.*Bx + By.*By + Bz.*Bz);
    B_normal = B_normal + Bnormal_from_plasma_current;
    mean_relative_err = mean(mean(abs(B_normal./B_mag)));
    max_relative_err = max(max(abs(B_normal./B_mag)));
    
    fprintf('mean_relative_err = %d',mean_relative_err);
    fprintf('max_relative_err = %d',max_relative_err);
    figure(6)
    contourf(B_normal)
    hold on
    title('Bnormal [T]')
    xlabel('zeta')
    ylabel('theta')
    colorbar() 
end

