function regcoil()

% This matlab script performs the same steps to the fortran program. The
% fortran and matlab versions are completely independent of each other. For
% identical inputs, they should give identical outputs to within roundoff
% error.

clear

symmetry_option = 3;
% 1 = sines only
% 2 = cosines only
% 3 = both sines and cosines

%load_bnorm = true;
load_bnorm = false;

bnorm_filename = '/Users/mattland/Box Sync/MATLAB/bnorm.d23p4_tm';

% This next value will be over-written if a VMEC equilibrium is used:
net_poloidal_current_Amperes = 1.4;
%net_toroidal_current_Amperes = 0.3;
net_toroidal_current_Amperes = 0.3;

% Resolution parameters:
% **********************************

ntheta_plasma = 32;
ntheta_coil   = 33;
nzeta_plasma = 34;
nzeta_coil   = 35;
mpol_coil  = 15;
ntor_coil  = 8;

%{
ntheta_plasma = 35;
ntheta_coil   = 34;
nzeta_plasma = 33;
nzeta_coil   = 32;
mpol_coil  = 16;
ntor_coil  = 15;
%}
%{
ntheta_plasma = 32;
ntheta_coil   = 32;
nzeta_plasma = 32;
nzeta_coil   = 32;
mpol_coil  = 8;
ntor_coil  = 8;
%}
%{
ntheta_plasma = 128;
ntheta_coil   = 128;
nzeta_plasma = 128;
nzeta_coil   = 128;
mpol_coil  = 32;
ntor_coil  = 32;
%}
% Options for the shape of the plasma surface:
% **********************************
geometry_option_plasma = 2;
R0_plasma = 3.0;
a_plasma = 1.0;
nfp_imposed = 1;
%woutFilename = 'C:\Users\landreman\Box Sync\MATLAB\20150601-01 Sfincs version 3\equilibria\wout_w7x_standardConfig.nc';
%woutFilename = '/Users/mattland/Box Sync/MATLAB/wout_d23p4_tm.nc';
woutFilename = 'equilibria/wout_d23p4_tm.nc';

% Options for the shape of the coil surface:
% **********************************
geometry_option_coil = 3;
R0_coil = 3.0;
a_coil = 1.7;
separation = 0.35;
%nescin_filename = 'nescin.w7x_standardConfig_separation0.3';
%nescin_filename = '/Users/mattland/Box Sync/MATLAB/nescin.w7x_winding_surface_from_Drevlak';
nescin_filename = 'equilibria/nescin.w7x_winding_surface_from_Drevlak';

% Options for the regularization parameter:
% **********************************
nlambda = 5;
lambda_min = 1e-30;
lambda_max = 1;

% Plotting options:
% **********************************

plot_results = true;
%plot_results = false;

max_nlambda_for_contour_plots = 18;

%plot3DFigure = true;
plot3DFigure = false;

plotGrids = true;
%plotGrids = false;

plotVectors = true;
%plotVectors = false;

%stopAfterInitialPlots = true;
stopAfterInitialPlots = false;

figureOffset = 0;


% Options related to checking fortran version
% *******************************************

compareToFortran = true;
%compareToFortran = false;

%fortranNcFilename = 'C:\Users\landreman\Box Sync\MATLAB\bdistrib_out.compareToMatlab.nc';
%fortranNcFilename = '/Users/mattland/regcoil/examples/compareToMatlab1/regcoil_out.compareToMatlab1.nc';
fortranNcFilename = 'examples/compareToMatlab2/regcoil_out.compareToMatlab2.nc';

fortranComparisonThreshhold_abs = 1e-11;

% *************************************************************************
% *************************************************************************
% End of input parameters.
% *************************************************************************
% *************************************************************************

mu0 = 4*pi*(1e-7);


    function compareVariableToFortran(variableName, varargin)
        % Specify 'abs' as an argument to compare the absolute values.
        % This is useful for the singular vectors, which are only defined
        % up to a sign in practice.
        if ~ compareToFortran
            return
        end
        try
            fortranVariable = double(ncread(fortranNcFilename,variableName));
        catch
            fprintf(['*** Variable ',variableName,' does not exist in the fortran output file.\n'])
            return
        end
        matlabVariable = eval(variableName);
        assignin('base',[variableName,'_m'],matlabVariable)
        assignin('base',[variableName,'_f'],fortranVariable)
        if isrow(matlabVariable)
            matlabVariable = matlabVariable(:);
        end
        if isrow(fortranVariable)
            fortranVariable = fortranVariable(:);
        end
        try
            % Next lines will cause an exception if sizes are different:
            if nargin>1 && strcmp(varargin{1},'abs')
                differences = abs(abs(matlabVariable) - abs(fortranVariable)) > fortranComparisonThreshhold_abs;
            else
                differences = abs(matlabVariable - fortranVariable) > fortranComparisonThreshhold_abs;
                %differences = (abs(matlabVariable - fortranVariable) > fortranComparisonThreshhold_abs) && ;
            end
            if any(any(any(differences)))
                fprintf(['*** Variable ',variableName,' is the same size Matlab and fortran but differs in value. max|diff|=%g\n'],max(max(max(differences))))
            else
                fprintf(['    Variable ',variableName,' is the same in Matlab and fortran.\n'])
            end
        catch
            fprintf(['*** Variable ',variableName,' is a different size between Matlab and fortran.\n'])
        end
    end

compareVariableToFortran('ntheta_plasma')
compareVariableToFortran('ntheta_coil')
compareVariableToFortran('nzeta_plasma')
compareVariableToFortran('nzeta_coil')
compareVariableToFortran('geometry_option_plasma')
compareVariableToFortran('geometry_option_coil')

% *********************************************
% Set up range of lambda to try
% *********************************************

lambda = zeros(nlambda,1);
lambda(2:end) = lambda_min * exp((0:(nlambda-2))/(nlambda-2)*log(lambda_max/lambda_min));

compareVariableToFortran('nlambda')
compareVariableToFortran('lambda')

% *********************************************
% Initialize the plasma surface:
% *********************************************

switch geometry_option_plasma
    case {0,1}
        % Plain axisymmetric circular torus
        nfp = nfp_imposed;
        mnmax = 2;
        xm = [0,1];
        xn = [0,0];
        rmnc = [R0_plasma; a_plasma];
        zmns = [0; a_plasma];
        whichSurface = 2;
        Rmajor_p = R0_plasma;
        
    case {2}
        % Load flux surface info from VMEC
        filename = woutFilename;
        ns = double(ncread(filename,'ns'));
        Rmajor_p = double(ncread(filename,'Rmajor_p'));
        nfp = double(ncread(filename,'nfp'));
        xm = double(ncread(filename,'xm'));
        xn = double(ncread(filename,'xn'));
        xm_nyq = double(ncread(filename,'xm_nyq'));
        xn_nyq = double(ncread(filename,'xn_nyq'));
        rmnc = double(ncread(filename,'rmnc'));
        zmns = double(ncread(filename,'zmns'));
        bmnc = double(ncread(filename,'bmnc'));
        mnmax = double(ncread(filename,'mnmax'));
        mnmax_nyq = double(ncread(filename,'mnmax_nyq'));
        
        whichSurface = ns; % Choose the outermost surface
        % Discard the other surfaces:
        rmnc = rmnc(:,whichSurface);
        zmns = zmns(:,whichSurface);
        
    otherwise
        error('Invalid setting for geometry_option_plasma')
end

switch geometry_option_plasma
    case {2}
        % BNORM scales B_n by curpol=(2*pi/nfp)*bsubv(m=0,n=0)
        % where bsubv is the extrapolation to the last full mesh point of
        % bsubvmnc.  Let's undo this scaling now.
        bsubvmnc = ncread(woutFilename,'bsubvmnc');
        bsubv00 = 1.5*bsubvmnc(1,end) - 0.5*bsubvmnc(1,end-1);
        curpol = 2*pi/nfp*bsubv00;  % /1 since nfp=1.
        
        bvco = ncread(woutFilename,'bvco');
        net_poloidal_current_Amperes = (2*pi/mu0) * (1.5*bvco(end) - 0.5*bvco(end-1));
        fprintf('New value for net_poloidal_current_Amperes: %g\n',net_poloidal_current_Amperes)
    otherwise
        curpol = 1;
end

compareVariableToFortran('net_poloidal_current_Amperes')
compareVariableToFortran('net_toroidal_current_Amperes')

nzetal_plasma = nzeta_plasma * nfp;
nzetal_coil   = nzeta_coil   * nfp;

theta_plasma = linspace(0,2*pi,ntheta_plasma+1);
theta_plasma(end) = [];
zeta_plasma = linspace(0,2*pi/nfp,nzeta_plasma+1);
zeta_plasma(end) = [];
zetal_plasma = linspace(0,2*pi,nzetal_plasma+1);
zetal_plasma(end) = [];
[zetal_plasma_2D, theta_plasma_2D] = meshgrid(zetal_plasma, theta_plasma);

x = zeros(ntheta_plasma,nzetal_plasma);
y = zeros(ntheta_plasma,nzetal_plasma);
z = zeros(ntheta_plasma,nzetal_plasma);

dxdtheta = zeros(ntheta_plasma,nzetal_plasma);
dydtheta = zeros(ntheta_plasma,nzetal_plasma);
dzdtheta = zeros(ntheta_plasma,nzetal_plasma);

dxdzeta = zeros(ntheta_plasma,nzetal_plasma);
dydzeta = zeros(ntheta_plasma,nzetal_plasma);
dzdzeta = zeros(ntheta_plasma,nzetal_plasma);

for i=1:mnmax
    angle = xm(i)*theta_plasma_2D-xn(i)*zetal_plasma_2D;
    angle2 = zetal_plasma_2D;
    
    x = x + rmnc(i)*cos(angle).*cos(angle2);
    y = y + rmnc(i)*cos(angle).*sin(angle2);
    z = z + zmns(i)*sin(angle);
    
    dxdtheta = dxdtheta - xm(i)*rmnc(i)*sin(angle).*cos(angle2);
    dydtheta = dydtheta - xm(i)*rmnc(i)*sin(angle).*sin(angle2);
    dzdtheta = dzdtheta + xm(i)*zmns(i)*cos(angle);
    
    dxdzeta = dxdzeta + xn(i)*rmnc(i)*sin(angle).*cos(angle2) ...
        - rmnc(i)*cos(angle).*sin(angle2);
    dydzeta = dydzeta + xn(i)*rmnc(i)*sin(angle).*sin(angle2) ...
        + rmnc(i)*cos(angle).*cos(angle2);
    dzdzeta = dzdzeta - xn(i)*zmns(i)*cos(angle);
end

Nx = dydzeta .* dzdtheta - dzdzeta .* dydtheta;
Ny = dzdzeta .* dxdtheta - dxdzeta .* dzdtheta;
Nz = dxdzeta .* dydtheta - dydzeta .* dxdtheta;
norm_normal_plasma = sqrt(Nx.*Nx + Ny.*Ny + Nz.*Nz);
norm_normal_plasma = norm_normal_plasma(:,1:nzeta_plasma);
dtheta_plasma = theta_plasma(2)-theta_plasma(1);
dzeta_plasma = zeta_plasma(2)-zeta_plasma(1);
area_plasma = sum(sum(norm_normal_plasma)) * dtheta_plasma * dzeta_plasma * nfp;

r_plasma = zeros(3, ntheta_plasma, nzetal_plasma);
drdtheta_plasma = zeros(3, ntheta_plasma, nzetal_plasma);
drdzeta_plasma = zeros(3, ntheta_plasma, nzetal_plasma);
normal_plasma = zeros(3, ntheta_plasma, nzetal_plasma);

r_plasma(1,:,:) = x;
r_plasma(2,:,:) = y;
r_plasma(3,:,:) = z;
drdtheta_plasma(1,:,:) = dxdtheta;
drdtheta_plasma(2,:,:) = dydtheta;
drdtheta_plasma(3,:,:) = dzdtheta;
drdzeta_plasma(1,:,:) = dxdzeta;
drdzeta_plasma(2,:,:) = dydzeta;
drdzeta_plasma(3,:,:) = dzdzeta;
normal_plasma(1,:,:) = Nx;
normal_plasma(2,:,:) = Ny;
normal_plasma(3,:,:) = Nz;

compareVariableToFortran('nfp')
compareVariableToFortran('theta_plasma')
compareVariableToFortran('zeta_plasma')
compareVariableToFortran('zetal_plasma')

compareVariableToFortran('r_plasma')
compareVariableToFortran('drdtheta_plasma')
compareVariableToFortran('drdzeta_plasma')
compareVariableToFortran('normal_plasma')
compareVariableToFortran('norm_normal_plasma')
compareVariableToFortran('area_plasma')

% *********************************************
% Initialize the coil surface:
% *********************************************

    function [theta, zeta, zetal, theta_2D, zetal_2D, r, drdtheta, drdzeta, d2rdtheta2, d2rdthetadzeta, d2rdzeta2, normal, norm_normal, area] ...
            = initSurface(ntheta, nzeta, geometry_option, R0, a, separation, nescin_filename)
        
        nzetal = nzeta*nfp;
        theta = linspace(0,2*pi,ntheta+1);
        theta(end) = [];
        zeta = linspace(0,2*pi/nfp,nzeta+1);
        zeta(end) = [];
        zetal = linspace(0,2*pi,nzetal+1);
        zetal(end) = [];
        
        dtheta = theta(2)-theta(1);
        dzeta  = zeta(2)-zeta(1);
        [zetal_2D, theta_2D] = meshgrid(zetal, theta);
        
        x = zeros(size(theta_2D));
        y = zeros(size(theta_2D));
        z = zeros(size(theta_2D));
        dxdtheta = zeros(size(theta_2D));
        dydtheta = zeros(size(theta_2D));
        dzdtheta = zeros(size(theta_2D));
        dxdzeta = zeros(size(theta_2D));
        dydzeta = zeros(size(theta_2D));
        dzdzeta = zeros(size(theta_2D));
        if false
            d2xdtheta2 = zeros(size(theta_2D));
            d2ydtheta2 = zeros(size(theta_2D));
            d2zdtheta2 = zeros(size(theta_2D));
            d2xdthetadzeta = zeros(size(theta_2D));
            d2ydthetadzeta = zeros(size(theta_2D));
            d2zdthetadzeta = zeros(size(theta_2D));
            d2xdzeta2 = zeros(size(theta_2D));
            d2ydzeta2 = zeros(size(theta_2D));
            d2zdzeta2 = zeros(size(theta_2D));
        end
        
        switch(geometry_option)
            case {0,1}
                if geometry_option == 0
                    R0_to_use = Rmajor_p;
                else
                    R0_to_use = R0;
                end
                
                x = (R0_to_use + a * cos(theta_2D)) .* cos(zetal_2D);
                y = (R0_to_use + a * cos(theta_2D)) .* sin(zetal_2D);
                z = a * sin(theta_2D);
                
                dxdtheta = -a * sin(theta_2D) .* cos(zetal_2D);
                dydtheta = -a * sin(theta_2D) .* sin(zetal_2D);
                dzdtheta = a * cos(theta_2D);
                
                dxdzeta = -(R0_to_use + a * cos(theta_2D)) .* sin(zetal_2D);
                dydzeta =  (R0_to_use + a * cos(theta_2D)) .* cos(zetal_2D);
                dzdzeta = zeros(size(theta_2D));
                
                if false
                    d2xdtheta2 = -a * cos(theta_2D) .* cos(zetal_2D);
                    d2ydtheta2 = -a * cos(theta_2D) .* sin(zetal_2D);
                    d2zdtheta2 = -a * sin(theta_2D);

                    d2xdthetadzeta =  a * sin(theta_2D) .* sin(zetal_2D);
                    d2ydthetadzeta = -a * sin(theta_2D) .* cos(zetal_2D);
                    d2zdthetadzeta = zeros(size(theta_2D));

                    d2xdzeta2 = -(R0_to_use + a * cos(theta_2D)) .* cos(zetal_2D);
                    d2ydzeta2 = -(R0_to_use + a * cos(theta_2D)) .* sin(zetal_2D);
                    d2zdzeta2 = zeros(size(theta_2D));
                end
            case 2
                error('geometry_option = 2 is not yet implemented for coil and outer surfaces.')
                
            case 3
                % Read coil surface from nescin file
                
                fid = fopen(nescin_filename,'r');
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
                fprintf('  Reading %d modes from nescin file %s\n',mnmax_nescin,nescin_filename)
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
                % Done reading nescin file.
                                
                for i = 1:mnmax_nescin
                    angle = xm_nescin(i)*theta_2D + xn_nescin(i)*zetal_2D*nfp;
                    angle2 = zetal_2D;
                    
                    x = x + rmnc_nescin(i)*cos(angle).*cos(angle2);
                    y = y + rmnc_nescin(i)*cos(angle).*sin(angle2);
                    z = z + zmns_nescin(i)*sin(angle);
                    
                    dxdtheta = dxdtheta - xm_nescin(i)*rmnc_nescin(i)*sin(angle).*cos(angle2);
                    dydtheta = dydtheta - xm_nescin(i)*rmnc_nescin(i)*sin(angle).*sin(angle2);
                    dzdtheta = dzdtheta + xm_nescin(i)*zmns_nescin(i)*cos(angle);
                    
                    dxdzeta = dxdzeta - nfp*xn_nescin(i)*rmnc_nescin(i)*sin(angle).*cos(angle2) ...
                        - rmnc_nescin(i)*cos(angle).*sin(angle2);
                    dydzeta = dydzeta - nfp*xn_nescin(i)*rmnc_nescin(i)*sin(angle).*sin(angle2) ...
                        + rmnc_nescin(i)*cos(angle).*cos(angle2);
                    dzdzeta = dzdzeta + nfp*xn_nescin(i)*zmns_nescin(i)*cos(angle);
                    
                end
                
            otherwise
                error('Invalid geometry_option')
        end
        
        Nx = dydzeta .* dzdtheta - dzdzeta .* dydtheta;
        Ny = dzdzeta .* dxdtheta - dxdzeta .* dzdtheta;
        Nz = dxdzeta .* dydtheta - dydzeta .* dxdtheta;
        
        r = zeros(3, ntheta, nzetal);
        drdtheta = zeros(3, ntheta, nzetal);
        drdzeta = zeros(3, ntheta, nzetal);
        normal = zeros(3, ntheta, nzetal);
        
        r(1,:,:) = x;
        r(2,:,:) = y;
        r(3,:,:) = z;
        drdtheta(1,:,:) = dxdtheta;
        drdtheta(2,:,:) = dydtheta;
        drdtheta(3,:,:) = dzdtheta;
        drdzeta(1,:,:) = dxdzeta;
        drdzeta(2,:,:) = dydzeta;
        drdzeta(3,:,:) = dzdzeta;
        normal(1,:,:) = Nx;
        normal(2,:,:) = Ny;
        normal(3,:,:) = Nz;
        
        if false
            d2rdtheta2 = zeros(3, ntheta, nzetal);
            d2rdthetadzeta = zeros(3, ntheta, nzetal);
            d2rdzeta2 = zeros(3, ntheta, nzetal);
            
            d2rdtheta2(1,:,:) = d2xdtheta2;
            d2rdtheta2(2,:,:) = d2ydtheta2;
            d2rdtheta2(3,:,:) = d2zdtheta2;
            
            d2rdthetadzeta(1,:,:) = d2xdthetadzeta;
            d2rdthetadzeta(2,:,:) = d2ydthetadzeta;
            d2rdthetadzeta(3,:,:) = d2zdthetadzeta;
            
            d2rdzeta2(1,:,:) = d2xdzeta2;
            d2rdzeta2(2,:,:) = d2ydzeta2;
            d2rdzeta2(3,:,:) = d2zdzeta2;
        else
            d2rdtheta2 = 0;
            d2rdthetadzeta = 0;
            d2rdzeta2 = 0;
        end
        
        norm_normal = sqrt(Nx.*Nx + Ny.*Ny + Nz.*Nz);
        norm_normal = norm_normal(:,1:nzeta);
        dtheta = theta(2)-theta(1);
        dzeta  = zeta(2)-zeta(1);
        area = sum(sum(norm_normal)) * dtheta * dzeta * nfp;
    end

tic
fprintf('Initializing coil surface.\n')
[theta_coil, zeta_coil, zetal_coil, theta_coil_2D, zetal_coil_2D, r_coil, drdtheta_coil, drdzeta_coil, ...
    d2rdtheta2_coil, d2rdthetadzeta_coil, d2rdzeta2_coil, normal_coil, norm_normal_coil, area_coil] ...
    = initSurface(ntheta_coil, nzeta_coil, geometry_option_coil, R0_coil, a_coil, separation, nescin_filename);
fprintf('Done. Took %g seconds.\n',toc)

compareVariableToFortran('theta_coil')
compareVariableToFortran('zeta_coil')
compareVariableToFortran('zetal_coil')
compareVariableToFortran('r_coil')
compareVariableToFortran('drdtheta_coil')
compareVariableToFortran('drdzeta_coil')
compareVariableToFortran('normal_coil')
compareVariableToFortran('norm_normal_coil')

% *********************************************
% Make 3D figure of surfaces
% *********************************************

%r_plasma = ncread(fortranNcFilename,'r_plasma');

if plot3DFigure
    r_plasma_toplot = r_plasma;
    r_coil_toplot = r_coil;
    
    % "Rotate" in theta so the seam in the plot is on the bottom
    nshift = round(ntheta_plasma*0.25);
    r_plasma_toplot = circshift(r_plasma_toplot, [0,nshift,0]);
    nshift = round(ntheta_coil*0.25);
    r_coil_toplot = circshift(r_coil_toplot, [0,nshift,0]);
    
    % Close surfaces for plotting:
    r_plasma_toplot(:,:,end+1) = r_plasma_toplot(:,:,1);
    r_plasma_toplot(:,end+1,:) = r_plasma_toplot(:,1,:);
    
    % For coil and outer surfaces, close in theta, but don't bother closing
    % in zeta:
    r_coil_toplot(:,end+1,:) = r_coil_toplot(:,1,:);
    
    
    
    mask = zetal_coil < 0.7*nfp;
    r_coil_toplot = r_coil_toplot(:,:,mask);
        
    figure(1 + figureOffset)
    clf
    set(gcf,'Color','w')
    faceColor = [1,0,0];
    surf(squeeze(r_plasma_toplot(1,:,:)),squeeze(r_plasma_toplot(2,:,:)),squeeze(r_plasma_toplot(3,:,:)),'FaceColor',faceColor,'EdgeColor','none')
    hold on
    if plotGrids
        plot3(squeeze(r_plasma_toplot(1,:,:)),squeeze(r_plasma_toplot(2,:,:)),squeeze(r_plasma_toplot(3,:,:)),'.r')
    end
    daspect([1,1,1])
    %shading interp
    axis vis3d
    hold on
    
    %faceColor = [1,0,1];
    faceColor = [0,1,0];
    surf(squeeze(r_coil_toplot(1,:,:)),squeeze(r_coil_toplot(2,:,:)),squeeze(r_coil_toplot(3,:,:)),'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.75)
    
    faceColor = [0,0,1];
    %surf(squeeze(r_outer_toplot(1,:,:)),squeeze(r_outer_toplot(2,:,:)),squeeze(r_outer_toplot(3,:,:)),'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.75)
    
    if plotGrids
        plot3(squeeze(r_coil_toplot(1,:,:)),squeeze(r_coil_toplot(2,:,:)),squeeze(r_coil_toplot(3,:,:)),'.m')
        %plot3(squeeze(r_outer_toplot(1,:,:)), squeeze(r_outer_toplot(2,:,:)), squeeze(r_outer_toplot(3,:,:)),'.b')
    end
    %surf(X_coil,Y_coil,Z_coil,'FaceColor',faceColor,'EdgeColor','none')
    %surf(X_coil,Y_coil,Z_coil,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.75)
    light
    lighting gouraud
    zoom(1.6)
    campos([  574.9370 -457.0244  424.3304])
    camva(1.0271)
    axis off
    
end

if stopAfterInitialPlots
    return
end

% *********************************************
% Set up Fourier arrays
% *********************************************

    function [mnmax, xm, xn] = setupFourierArrays(mpol,ntor)
        % xm is non-negative, while xn can be negative
        % xn is the rapidly increasing variable.
        
        % When xm=0, xn=1..ntor.
        % When xm>0, xn=-ntor..ntor.
        mnmax = ntor + mpol*(2*ntor+1);
        
        xm = zeros(mnmax,1);
        xn = zeros(mnmax,1);
        xn(1:ntor) = 1:ntor;
        nextIndex = ntor+1;
        for m = 1:mpol
            indices = nextIndex:(nextIndex+2*ntor);
            xm(indices) = m;
            xn(indices) = (-ntor):ntor;
            nextIndex = nextIndex + 2*ntor+1;
        end
        
    end

[mnmax_coil, xm_coil, xn_coil] = setupFourierArrays(mpol_coil, ntor_coil);
xn_coil = xn_coil * nfp;

compareVariableToFortran('mpol_coil')
compareVariableToFortran('ntor_coil')
compareVariableToFortran('mnmax_coil')
compareVariableToFortran('xm_coil')
compareVariableToFortran('xn_coil')
compareVariableToFortran('symmetry_option')

% *********************************************
% Load BNORM data.
% *********************************************

Bnormal_from_plasma_current = zeros(ntheta_plasma,nzeta_plasma);
[zeta_plasma_2D, theta_plasma_2D] = meshgrid(zeta_plasma, theta_plasma);
if load_bnorm
    fid = fopen(bnorm_filename,'r');
    if fid<0
        error('Unable to open BNORM file %s.\n',bnorm_filename)
    end
    while ~ feof(fid);
        [fileline,count] = fscanf(fid,'%f %f %f\n',3);
        if count == 3
            mm  = fileline(1);
            nn  = fileline(2);
            amp = fileline(3);
            Bnormal_from_plasma_current = Bnormal_from_plasma_current + amp*sin(mm*theta_plasma_2D + nfp*nn*zeta_plasma_2D);
        end
    end
else
end
Bnormal_from_plasma_current = Bnormal_from_plasma_current * curpol;
Bnormal_from_plasma_current_1D = reshape(Bnormal_from_plasma_current, [ntheta_plasma*nzeta_plasma,1]);
compareVariableToFortran('Bnormal_from_plasma_current')

% *********************************************
% Compute h
% *********************************************
tic
fprintf('Computing h.\n')
h = zeros(ntheta_plasma,nzeta_plasma);

zeta_coil_indices = 1:nzeta_coil;
G_drdtheta_minus_I_drdzeta_x = squeeze(net_poloidal_current_Amperes * drdtheta_coil(1,:,:) - net_toroidal_current_Amperes * drdzeta_coil(1,:,:));
G_drdtheta_minus_I_drdzeta_y = squeeze(net_poloidal_current_Amperes * drdtheta_coil(2,:,:) - net_toroidal_current_Amperes * drdzeta_coil(2,:,:));
G_drdtheta_minus_I_drdzeta_z = squeeze(net_poloidal_current_Amperes * drdtheta_coil(3,:,:) - net_toroidal_current_Amperes * drdzeta_coil(3,:,:));

d_x = reshape(G_drdtheta_minus_I_drdzeta_x(:,zeta_coil_indices) / (2*pi), [ntheta_coil*nzeta_coil,1]);
d_y = reshape(G_drdtheta_minus_I_drdzeta_y(:,zeta_coil_indices) / (2*pi), [ntheta_coil*nzeta_coil,1]);
d_z = reshape(G_drdtheta_minus_I_drdzeta_z(:,zeta_coil_indices) / (2*pi), [ntheta_coil*nzeta_coil,1]);

for itheta_plasma = 1:ntheta_plasma
    for izeta_plasma = 1:nzeta_plasma
        adx = r_plasma(1,itheta_plasma,izeta_plasma) - squeeze(r_coil(1,:,:));
        ady = r_plasma(2,itheta_plasma,izeta_plasma) - squeeze(r_coil(2,:,:));
        adz = r_plasma(3,itheta_plasma,izeta_plasma) - squeeze(r_coil(3,:,:));
        adr2 = adx.*adx + ady.*ady + adz.*adz;
        dr32 = adr2 .* sqrt(adr2);
        tempMatrix = ( ...
              G_drdtheta_minus_I_drdzeta_x.*ady*normal_plasma(3,itheta_plasma,izeta_plasma) ...
            + G_drdtheta_minus_I_drdzeta_y.*adz*normal_plasma(1,itheta_plasma,izeta_plasma) ...
            + G_drdtheta_minus_I_drdzeta_z.*adx*normal_plasma(2,itheta_plasma,izeta_plasma) ...
            - G_drdtheta_minus_I_drdzeta_z.*ady*normal_plasma(1,itheta_plasma,izeta_plasma) ...
            - G_drdtheta_minus_I_drdzeta_x.*adz*normal_plasma(2,itheta_plasma,izeta_plasma) ...
            - G_drdtheta_minus_I_drdzeta_y.*adx*normal_plasma(3,itheta_plasma,izeta_plasma)) ./ dr32;
        
        h(itheta_plasma,izeta_plasma) = sum(sum(tempMatrix));
    end
end

dtheta_coil = theta_coil(2)-theta_coil(1);
dzeta_coil = zeta_coil(2)-zeta_coil(1);
h = h * (dtheta_coil*dzeta_coil*mu0/(8*pi*pi));
Bnormal_from_net_coil_currents = h ./ norm_normal_plasma;
Bnormal_from_net_coil_currents_1D = reshape(Bnormal_from_net_coil_currents, [ntheta_plasma*nzeta_plasma,1]);
fprintf('Done. Took %g seconds.\n',toc)
compareVariableToFortran('Bnormal_from_net_coil_currents')

% ***********************************************
% Compute the basis functions and f on the (theta,zeta) grids.
% ***********************************************

switch symmetry_option
    case {1,2}
        num_basis_functions = mnmax_coil;
    case {3}
        num_basis_functions = mnmax_coil * 2;
    otherwise
        error('Invalid value for symmetry_option')
end
basis_functions = zeros(ntheta_coil*nzeta_coil, num_basis_functions);
f_x = zeros(ntheta_coil*nzeta_coil, num_basis_functions);
f_y = zeros(ntheta_coil*nzeta_coil, num_basis_functions);
f_z = zeros(ntheta_coil*nzeta_coil, num_basis_functions);

fprintf('Computing Fourier functions and f.\n')
tic
[zeta_coil_2D, theta_coil_2D] = meshgrid(zeta_coil,theta_coil);
zeta_coil_indices = 1:nzeta_coil;
switch symmetry_option
    case {1}
        % sines only
        for imn = 1:mnmax_coil
            angle = xm_coil(imn)*theta_coil_2D - xn_coil(imn)*zeta_coil_2D;
            cosangle = cos(angle);
            sinangle = sin(angle);
            basis_functions(:,imn) = reshape(sinangle, [ntheta_coil*nzeta_coil,1]);
            f_x(:,imn) = reshape(cosangle.*squeeze(xn_coil(imn)*drdtheta_coil(1,:,zeta_coil_indices) + xm_coil(imn)*drdzeta_coil(1,:,zeta_coil_indices)), [ntheta_coil*nzeta_coil,1]);
            f_y(:,imn) = reshape(cosangle.*squeeze(xn_coil(imn)*drdtheta_coil(2,:,zeta_coil_indices) + xm_coil(imn)*drdzeta_coil(2,:,zeta_coil_indices)), [ntheta_coil*nzeta_coil,1]);
            f_z(:,imn) = reshape(cosangle.*squeeze(xn_coil(imn)*drdtheta_coil(3,:,zeta_coil_indices) + xm_coil(imn)*drdzeta_coil(3,:,zeta_coil_indices)), [ntheta_coil*nzeta_coil,1]);
        end
    case {2}
        % cosines only
        for imn = 1:mnmax_coil
            angle = xm_coil(imn)*theta_coil_2D - xn_coil(imn)*zeta_coil_2D;
            cosangle = cos(angle);
            sinangle = sin(angle);
            basis_functions(:,imn) = reshape(cosangle, [ntheta_coil*nzeta_coil,1]);
            f_x(:,imn) = reshape(-sinangle.*squeeze(xn_coil(imn)*drdtheta_coil(1,:,zeta_coil_indices) + xm_coil(imn)*drdzeta_coil(1,:,zeta_coil_indices)), [ntheta_coil*nzeta_coil,1]);
            f_y(:,imn) = reshape(-sinangle.*squeeze(xn_coil(imn)*drdtheta_coil(2,:,zeta_coil_indices) + xm_coil(imn)*drdzeta_coil(2,:,zeta_coil_indices)), [ntheta_coil*nzeta_coil,1]);
            f_z(:,imn) = reshape(-sinangle.*squeeze(xn_coil(imn)*drdtheta_coil(3,:,zeta_coil_indices) + xm_coil(imn)*drdzeta_coil(3,:,zeta_coil_indices)), [ntheta_coil*nzeta_coil,1]);
        end
    case {3}
        % Both sines and cosines
        for imn = 1:mnmax_coil
            angle = xm_coil(imn)*theta_coil_2D - xn_coil(imn)*zeta_coil_2D;
            cosangle = cos(angle);
            sinangle = sin(angle);
            basis_functions(:,imn)            = reshape(sinangle, [ntheta_coil*nzeta_coil,1]);
            basis_functions(:,imn+mnmax_coil) = reshape(cosangle, [ntheta_coil*nzeta_coil,1]);
            temparr = squeeze(xn_coil(imn)*drdtheta_coil(1,:,zeta_coil_indices) + xm_coil(imn)*drdzeta_coil(1,:,zeta_coil_indices));
            f_x(:,imn)            = reshape( cosangle.*temparr, [ntheta_coil*nzeta_coil,1]);
            f_x(:,imn+mnmax_coil) = reshape(-sinangle.*temparr, [ntheta_coil*nzeta_coil,1]);
            temparr = squeeze(xn_coil(imn)*drdtheta_coil(2,:,zeta_coil_indices) + xm_coil(imn)*drdzeta_coil(2,:,zeta_coil_indices));
            f_y(:,imn)            = reshape( cosangle.*temparr, [ntheta_coil*nzeta_coil,1]);
            f_y(:,imn+mnmax_coil) = reshape(-sinangle.*temparr, [ntheta_coil*nzeta_coil,1]);
            temparr = squeeze(xn_coil(imn)*drdtheta_coil(3,:,zeta_coil_indices) + xm_coil(imn)*drdzeta_coil(3,:,zeta_coil_indices));
            f_z(:,imn)            = reshape( cosangle.*temparr, [ntheta_coil*nzeta_coil,1]);
            f_z(:,imn+mnmax_coil) = reshape(-sinangle.*temparr, [ntheta_coil*nzeta_coil,1]);
        end
end
fprintf('Done. Took %g sec.\n',toc)
        


% *********************************************
% Compute g
% *********************************************
%return

fprintf('Computing inductance\n')
tic

zeta_plasma_indices = 1:nzeta_plasma;
inductance = zeros(ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil);
for itheta_coil = 1:ntheta_coil
    for izeta_coil = 1:nzeta_coil
        index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil;
        for l_coil = 0:(nfp-1)
            izetal_coil = izeta_coil + l_coil*nzeta_coil;
            dx = r_plasma(1,:,zeta_plasma_indices) - r_coil(1,itheta_coil,izetal_coil);
            dy = r_plasma(2,:,zeta_plasma_indices) - r_coil(2,itheta_coil,izetal_coil);
            dz = r_plasma(3,:,zeta_plasma_indices) - r_coil(3,itheta_coil,izetal_coil);
            dr2 = dx.*dx + dy.*dy + dz.*dz;
            denominator = dr2 .* sqrt(dr2);
            temp = (normal_plasma(1,:,zeta_plasma_indices)*normal_coil(1,itheta_coil,izetal_coil) ...
                +   normal_plasma(2,:,zeta_plasma_indices)*normal_coil(2,itheta_coil,izetal_coil) ...
                +   normal_plasma(3,:,zeta_plasma_indices)*normal_coil(3,itheta_coil,izetal_coil) ...
                - (3./dr2) .* (dx .* normal_plasma(1,:,zeta_plasma_indices) + dy .* normal_plasma(2,:,zeta_plasma_indices) + dz .* normal_plasma(3,:,zeta_plasma_indices)) ...
                .* (dx * normal_coil(1,itheta_coil,izetal_coil) + dy * normal_coil(2,itheta_coil,izetal_coil) + dz * normal_coil(3,itheta_coil,izetal_coil))) ./ denominator;
            inductance(:,index_coil) = inductance(:,index_coil) + ...
                reshape(temp, [ntheta_plasma*nzeta_plasma,1]);
        end
    end
end
inductance = inductance * (mu0/(4*pi));
fprintf('Done. Took %g sec.\n',toc)

compareVariableToFortran('inductance')

tic1 = tic;
g = (dtheta_coil * dzeta_coil) * inductance * basis_functions;
fprintf('Matmul: %g\n',toc(tic1))


compareVariableToFortran('g')
%return
% *********************************************
% Compute matrices and RHS for the normal equations:
% *********************************************

norm_normal_plasma_vec = reshape(norm_normal_plasma,[ntheta_plasma*nzeta_plasma,1]);
norm_normal_coil_vec   = reshape(norm_normal_coil,  [ntheta_coil*nzeta_coil,    1]);
diag_inv_norm_normal_plasma = diag(1./norm_normal_plasma_vec);
diag_inv_norm_normal_coil   = diag(1./norm_normal_coil_vec);

tic
fprintf('Computing RHS_B and RHS_K.\n')
RHS_B = (-dtheta_plasma*dzeta_plasma)*((Bnormal_from_plasma_current_1D + Bnormal_from_net_coil_currents_1D)' * g)';
RHS_K = (dtheta_coil*dzeta_coil)*((d_x ./ norm_normal_coil_vec)' * f_x + (d_y ./ norm_normal_coil_vec)' * f_y + (d_z ./ norm_normal_coil_vec)' * f_z)';
fprintf('Done. Took %g sec.\n',toc)

compareVariableToFortran('RHS_B')
compareVariableToFortran('RHS_K')

tic
fprintf('Computing matrix_B.\n')
matrix_B = (dtheta_plasma*dzeta_plasma)*( (g') * diag_inv_norm_normal_plasma * g );
fprintf('Done. Took %g sec.\n',toc)
compareVariableToFortran('matrix_B')

tic
fprintf('Computing matrix_K.\n')
matrix_K = (dtheta_coil*dzeta_coil)*( ...
    (f_x') * diag_inv_norm_normal_coil * f_x ...
    +(f_y') * diag_inv_norm_normal_coil * f_y ...
    +(f_z') * diag_inv_norm_normal_coil * f_z );
fprintf('Done. Took %g sec.\n',toc)
compareVariableToFortran('matrix_K')


% *********************************************
% Solve the system for each lambda:
% *********************************************

single_valued_current_potential_mn = zeros(num_basis_functions, nlambda);
single_valued_current_potential_thetazeta = zeros(ntheta_coil, nzeta_coil, nlambda);
current_potential = zeros(ntheta_coil, nzeta_coil, nlambda);
[zeta_coil_2D, theta_coil_2D] = meshgrid(zeta_coil, theta_coil);
chi2_B = zeros(nlambda,1);
chi2_K = zeros(nlambda,1);
Bnormal_total = zeros(ntheta_plasma, nzeta_plasma, nlambda);
K2 = zeros(ntheta_coil, nzeta_coil, nlambda);

for ilambda=1:nlambda
    fprintf('Solving system for lambda = %g  (%d of %d)\n',lambda(ilambda), ilambda, nlambda)
    
    tic
    matrix = matrix_B + lambda(ilambda) * matrix_K;
    RHS    = RHS_B    + lambda(ilambda) * RHS_K;
    %matrix = matrix_B - lambda(ilambda) * matrix_K;
    %RHS    = RHS_B    - lambda(ilambda) * RHS_K;
    fprintf('  Summing matrices: %g sec.\n',toc)
    
    tic
    solution = matrix \ RHS;
    fprintf('  Solve: %g sec.\n',toc)
    
    tic
    single_valued_current_potential_mn(:,ilambda) = solution;
    this_single_valued_current_potential_thetazeta = reshape(basis_functions*solution, [ntheta_coil,nzeta_coil]);
    single_valued_current_potential_thetazeta(:,:,ilambda) = this_single_valued_current_potential_thetazeta;
    current_potential(:,:,ilambda) = this_single_valued_current_potential_thetazeta ...
        + zeta_coil_2D * (net_poloidal_current_Amperes/(2*pi)) ...
        + theta_coil_2D * (net_toroidal_current_Amperes/(2*pi));
    
    this_Bnormal = Bnormal_from_plasma_current + Bnormal_from_net_coil_currents ...
        + reshape(g*solution, [ntheta_plasma,nzeta_plasma]) ./ norm_normal_plasma;
    Bnormal_total(:,:,ilambda) = this_Bnormal;
    chi2_B(ilambda) = nfp*dtheta_plasma*dzeta_plasma*sum(sum(this_Bnormal .* this_Bnormal .* norm_normal_plasma));
    
    K_difference_x = d_x - f_x*solution;
    K_difference_y = d_y - f_y*solution;
    K_difference_z = d_z - f_z*solution;
    this_K2_over_N = reshape(K_difference_x.*K_difference_x + K_difference_y.*K_difference_y + K_difference_z.*K_difference_z, [ntheta_coil, nzeta_coil]) ./(norm_normal_coil);
    K2(:,:,ilambda) = this_K2_over_N ./ norm_normal_coil;
    chi2_K(ilambda) = nfp*dtheta_coil*dzeta_coil*sum(sum(this_K2_over_N));
    
    fprintf('  Diagnostics: %g sec.\n',toc)
    fprintf('  chi2_B: %g,   chi2_K: %g\n',chi2_B(ilambda),chi2_K(ilambda))
end

compareVariableToFortran('single_valued_current_potential_mn')
compareVariableToFortran('single_valued_current_potential_thetazeta')
compareVariableToFortran('current_potential')
compareVariableToFortran('Bnormal_total')
compareVariableToFortran('K2')
compareVariableToFortran('chi2_B')
compareVariableToFortran('chi2_K')

%return


% *********************************************
% Done with the main calculation.
% Now plot results.
% *********************************************

if ~ plot_results
    return
end

figure(2)
clf
numRows=2;
numCols=3;

subplot(numRows,numCols,1)
loglog(chi2_K, chi2_B,'o-')
xlabel('chi2 K')
ylabel('chi2 B')

subplot(numRows,numCols,2)
loglog(lambda, chi2_B,'o-')
xlabel('lambda')
ylabel('chi2 B')

subplot(numRows,numCols,3)
semilogy(lambda, chi2_B,'o-')
xlabel('lambda')
ylabel('chi2 B')

subplot(numRows,numCols,5)
loglog(lambda, chi2_K,'o-')
xlabel('lambda')
ylabel('chi2 K')

subplot(numRows,numCols,6)
semilogy(lambda, chi2_K,'o-')
xlabel('lambda')
ylabel('chi2 K')

% ***********************************************************************
% Plot single-valued part of the current potential

figure(3)
clf
numContours = 25;

numPlots = min([max_nlambda_for_contour_plots,nlambda]);
ilambda_to_plot = unique(round(linspace(1,nlambda,numPlots)));
numPlots = numel(ilambda_to_plot);

numCols=ceil(sqrt(numPlots));
numRows=ceil(numPlots / numCols);

for iplot = 1:numPlots
    subplot(numRows,numCols,iplot)
    contourf(zeta_coil_2D, theta_coil_2D, single_valued_current_potential_thetazeta(:,:,ilambda_to_plot(iplot)), numContours,'EdgeColor','none')
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title(['Single valued current potential for lambda=',num2str(lambda(ilambda_to_plot(iplot)))])
end
    

% ***********************************************************************
% Plot full current potential

figure(4)
clf

for iplot = 1:numPlots
    subplot(numRows,numCols,iplot)
    contourf(zeta_coil_2D, theta_coil_2D, current_potential(:,:,ilambda_to_plot(iplot)), numContours,'EdgeColor','none')
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title(['Current potential for lambda=',num2str(lambda(ilambda_to_plot(iplot)))])
end
    
% ***********************************************************************
% Plot K^2

figure(5)
clf

for iplot = 1:numPlots
    subplot(numRows,numCols,iplot)
    contourf(zeta_coil_2D, theta_coil_2D, K2(:,:,ilambda_to_plot(iplot)), numContours,'EdgeColor','none')
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title(['K^2 for lambda=',num2str(lambda(ilambda_to_plot(iplot)))])
end

% ***********************************************************************
% Plot B_normal

figure(6)
clf

numCols=ceil(sqrt(numPlots+2));
numRows=ceil((numPlots+2) / numCols);

subplot(numRows,numCols,1)
contourf(zeta_plasma_2D, theta_plasma_2D, Bnormal_from_plasma_current, numContours,'EdgeColor','none')
colorbar
xlabel('zeta')
ylabel('theta')
title('Bnormal from plasma current')

subplot(numRows,numCols,2)
contourf(zeta_plasma_2D, theta_plasma_2D, Bnormal_from_net_coil_currents, numContours,'EdgeColor','none')
colorbar
xlabel('zeta')
ylabel('theta')
title('Bnormal from net coil currents')

for iplot = 1:numPlots
    subplot(numRows,numCols,iplot+2)
    contourf(zeta_plasma_2D, theta_plasma_2D, Bnormal_total(:,:,ilambda_to_plot(iplot)), numContours,'EdgeColor','none')
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title(['Total Bnormal for lambda=',num2str(lambda(ilambda_to_plot(iplot)))])
end

end

%    stringForTop = ['Singular vectors v of the transfer matrix: coil surface (threshold=',num2str(pseudoinverse_thresholds(whichThreshold)),')'];
%    annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
%        'Interpreter','none','VerticalAlignment','bottom',...
%        'FontSize',11,'LineStyle','none','String',stringForTop);
