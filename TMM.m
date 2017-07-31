% This code calculates Transfer Matrix Method, which was described by Raymond 
% C. Rumpf at the course Computational EM. I also used code of
% Sathyanarayan Rao from https://www.mathworks.com/matlabcentral/fileexchange/47637-transmittance-and-reflectance-spectra-of-multilayered-dielectric-stack-using-transfer-matrix-method
% 
clear all, %close all, %clc;
 

I = eye(2);
% units 
degrees = pi/180;

%% parameters
% device
UR = [ 1 1 1 1 1 1 1 1 1 1]; % array of permeabilities in each layer 
ER = [ 3 1 3 1 3 1 3 1 3 1 ]; % array of permittivities in each layer  
L  = [ 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5 0.25 0.5].*1e-6; 
N  = 10;

% external
er1 = 1.0; % reflection side medium
ur1 = 1.0;
er2 = 5.0;  % transition side medium
ur2 = 1.0;
ni  = sqrt(er1); % refrctive index of incidence media

%source
lam0  =(400:1:800)*1e-9;%free space wavelength 
k0    =(2*pi)./lam0;
theta = 17 * degrees; %elevation angle 
phi   = 12 * degrees; %azimuthal angle 
pte   = 1/sqrt(2);
ptm   = 1i*pte;

%% TMM algorithm
% Transverse Wave Vectors
kx = ni*sin(theta)*cos(phi);
ky = ni*sin(theta)*sin(phi);
kz = sqrt(1-(kx*kx)-(ky*ky));

% homogeneous Gap Medium Parameters
Wh = I;
% Qh = [kx*ky 1-kx*kx; ky*ky-1 -kx*ky]; % Q
Qh=[ kx*ky 1-(kx*kx) ;(ky*ky)-1  -kx*ky];
Omh = 1i*kz*Wh; % Omega
Vh = Qh * (Omh^-1);

% initialaze Global Scattering Matrix
Sg11 = zeros(2); 
Sg12 = I;
Sg21 = I;
Sg22 = zeros(2);

%% reflection side

krz = sqrt(ur1*er1 - (kx*kx) - (ky*ky));
Qr = 1/ur1 * [kx*ky ur1*er1-kx*kx; ky*ky-ur1*er1 -kx*ky];
Omr = 1i*krz*I;
Vr = Qr*(Omr^-1);
Ar = I + (Vh^-1)*Vr; %!!!!!!!!!!!! check out this place later
Br = I - (Vh^-1)*Vr; 
Sr11 = -Ar^-1 * Br;
Sr12 = 2 * Ar^-1;
Sr21 = 0.5 * (Ar - Br*Ar^-1 * Br);
Sr22 = Br*Ar^-1;

% connect external reflection region to device

[Sg11,Sg12,Sg21,Sg22] = star_product(Sr11,Sr12,Sr21,Sr22,Sg11,Sg12,Sg21,Sg22);


for j = 1:length(lam0)
    % reflection + each layer for every wavelength
    for i = 1:N
        kz = sqrt(ER(i)*UR(i) - (kx*kx) - (ky*ky));
        Q = 1/UR(i) * [kx*ky UR(i)*ER(i)-kx*kx; ky*ky-UR(i)*ER(i) -kx*ky];
        Om = 1i*kz*I;
        V = Q*(Om^-1);
        A = I + (V^-1)*Vh;
        B = I - (V^-1)*Vh;
        X = expm(Om*k0(j)*L(i));
        D = A - X*B*A^-1*X*B;
        S11 = D^-1*(X*B*A^-1*X*A-B);
        S22 = S11;
        S12 = D^-1*X*(A-B*A^-1*B);
        S21 = S12;
        % update gloabal S-matrix
        
        [Sg11,Sg12,Sg21,Sg22] = star_product(Sg11,Sg12,Sg21,Sg22,S11,S12,S21,S22);
    end
    % now we add transmission
    ktz = sqrt(ur2*er2 - (kx*kx)-(ky*ky));
    Qt = 1/ur2 * [kx*ky ur2*er2-(kx*kx); (ky*ky)-ur2*er2 -kx*ky];
    Omt = 1i*ktz*I;
    Vt = Qt*(Omt^-1);
    At = I + (Vh^-1)*Vt;
    Bt = I - (Vh^-1)*Vt;
    St11 = Bt*At^-1;
    St12 = 0.5*(At-Bt*At^-1*Bt);
    St21 = 2*At^-1;
    St22 = -(At^-1)*Bt;
    
    [Sf11,Sf12,Sf21,Sf22] = star_product(Sg11,Sg12,Sg21,Sg22,St11,St12,St21,St22);
    
    % Source
    Kinc = k0(j)*ni*[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
    ncn = [0 0 -1]; % surface normal
    if theta == 0
        aTE = [0 1 0];
    else
        aTE = cross(Kinc,ncn)./norm(cross(Kinc,ncn));
    end
    aTM = cross(aTE,Kinc)./norm(cross(aTE,Kinc));
    P = pte*aTE + ptm*aTM;
    cinc = [P(1); P(2)];
    % tranmitted and reflected fields
    Er=Sf11*cinc;
    Et=Sf21*cinc;
    Erx=Er(1);
    Ery=Er(2);
    Etx=Et(1);
    Ety=Et(2);
    
    % Longitudial field components
    Erz = -(kx*Erx + ky*Ery)/krz;
    Etz = -(kx*Etx + ky*Ety)/ktz;
    
    % Transmittance and reflectance
    R = abs(Erx).^2 + abs(Ery).^2 + abs(Erz).^2;
    T = (abs(Etx)^2 + abs(Ety)^2 + abs(Etz)^2)*real((ur1*ktz)/(ur2*krz));
    
    Tx(j) = abs(T);
    Rx(j) = abs(R);
    
end

% plot result
figure(1)
plot(lam0,Tx,lam0,Rx,'r','linewidth',2)
axis([400e-9 800e-9 -0.05 1.05])
xlabel('\lambda','fontsize',14)
ylabel('T (in blue) or R (in red)','fontsize',14)
title('Transmittance (T) and Reflectance(R)','fontsize',14)
fh = figure(1);
set(fh, 'color', 'white');

function [SAB11,SAB12,SAB21,SAB22] = star_product (SA11,SA12,SA21,SA22,SB11,SB12,SB21,SB22)
% calculates star procudtion for two scattering matrix
    I = eye(2);
    D = SA12 * (I - SB11*SA22)^-1;
    F = SB21 * (I - SA22*SB11)^-1;

    SAB11 = SA11 + D*SB11*SA21;
    SAB12 = D*SB12;
    SAB21 = F*SA21;
    SAB22 = SB22 + F*SA22*SB12;
end