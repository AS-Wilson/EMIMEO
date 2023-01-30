function T_R_A=tmm(lam,theta,n_layers,d_layers,n_input,n_output,pol)

% This MATLAB function calculates:
% Transmittance (Transmitted power divided by Input power);
% Reflectance (Reflected power divided by incident power),
% and Absorption (Absorbed power divided by incident power).
% R + T + A = 1 because of energy conservation.
% The output T_R_A is a vector with three values:
% T_R_A(1) = Transmittance
% T_R_A(2) = Reflectance
% T_R_A(3) = Absorption
%
% INPUT PARAMETERS %
% lam: wavelength in meters (this is a scalar quantity). Example: lam=532e-9
% theta: angle in degrees (this is a scalar quantity). Example: theta=30
% n_layers: VECTOR of (complex) refractive indices in each layer excluding the input and
% output semi-infinite substrates. The first index refers to the layer
% adjacent to n_input. IMPORTANT NOTE: The imaginary part of the refractive index
% must have positive sign to indicate damping!!!
% d_layers: VECTOR of layers' thicknesses in meters
% n_input: (real-valued) index of input medium. This is a scalar quantity
% n_output: (real-valued) index of output medium. This is a scalar quantity
% pol: 0=TE (s-polarization), 1=TM (p-polarization)
%%%%%%%%%%%%%%%%%%%%
%
%%%% SETUP %%%%%%%
%
% Incident plane wave
%                   \      |
%                    \     |
%                     \    |   Input Region  (n_input)
%                      \   |
%                       \  |
%                        \ |
%       --------------------------------------------
%   d(1)                  n(1)
%       --------------------------------------------
%   d(2)                  n(2)
%       --------------------------------------------                                                |
%                         ....
%       --------------------------------------------
%   d(N)                  n(N)
%       --------------------------------------------
%
%                     Output Region (n_output)
%
%


N_layers=length(d_layers);
eta0=120*pi; % Vacuum impedance
k0=2*pi/lam; % Wavenumber in vacuo
sin_theta=sin(theta*pi/180);  % sine of angle of incidence

p_input=cos(theta*pi/180)/eta0*n_input;
p_output=sqrt(k0^2*n_output^2-k0^2*n_input^2*sin_theta^2)/k0/eta0;
q_input=cos(theta*pi/180)*eta0/n_input;
q_output=eta0*sqrt(k0^2*n_output^2-k0^2*n_input^2*sin_theta^2)/k0/n_output^2;

% Matrix inizialization
M_te=diag(ones(2,1)); M_tm=M_te;

if pol==1
    for ii=1:N_layers
        p=sqrt(k0^2*n_layers(ii)^2-k0^2*n_input^2*sin_theta^2)/k0/eta0;
        q=eta0*sqrt(k0^2*n_layers(ii)^2-k0^2*n_input^2*sin_theta^2)/k0/n_layers(ii)^2;
        tm(1,1)=cos(d_layers(ii)*sqrt(k0^2*n_layers(ii)^2-k0^2*n_input^2*sin_theta^2));
        tm(2,2)=tm(1,1);
        tm(1,2)=-1i/q*sin(d_layers(ii)*sqrt(k0^2*n_layers(ii)^2-k0^2*n_input^2*sin_theta^2));
        tm(2,1)=q^2*tm(1,2);
        M_tm=M_tm*tm;
    end
    r_tm=((M_tm(1,1)+M_tm(1,2)*q_output)*q_input-(M_tm(2,1)+M_tm(2,2)*q_output))/...
        ((M_tm(1,1)+M_tm(1,2)*q_output)*q_input+(M_tm(2,1)+M_tm(2,2)*q_output));
    t_tm=2*q_input/...
        ((M_tm(1,1)+M_tm(1,2)*q_output)*q_input+(M_tm(2,1)+M_tm(2,2)*q_output));
    
    T_tm=real(abs(t_tm)^2*n_input/n_output*cos(asin(sin(theta*pi/180)*n_input/n_output))/cos(theta*pi/180));
    R_tm=abs(r_tm)^2; A_tm=1-R_tm-T_tm;
    
    T_R_A(1)=T_tm;
    T_R_A(2)=R_tm;
    T_R_A(3)=A_tm;
    
elseif pol==0
    
    for ii=1:N_layers
        p=sqrt(k0^2*n_layers(ii)^2-k0^2*n_input^2*sin_theta^2)/k0/eta0;
        q=eta0*sqrt(k0^2*n_layers(ii)^2-k0^2*n_input^2*sin_theta^2)/k0/n_layers(ii)^2;
        te(1,1)=cos(d_layers(ii)*sqrt(k0^2*n_layers(ii)^2-k0^2*n_input^2*sin_theta^2));
        te(2,2)=te(1,1);
        te(1,2)=-1i/p*sin(d_layers(ii)*sqrt(k0^2*n_layers(ii)^2-k0^2*n_input^2*sin_theta^2));
        te(2,1)=p^2*te(1,2);
        M_te=M_te*te;
    end
    r_te=((M_te(1,1)+M_te(1,2)*p_output)*p_input-(M_te(2,1)+M_te(2,2)*p_output))/...
        ((M_te(1,1)+M_te(1,2)*p_output)*p_input+(M_te(2,1)+M_te(2,2)*p_output));
    t_te=2*p_input/...
        ((M_te(1,1)+M_te(1,2)*p_output)*p_input+(M_te(2,1)+M_te(2,2)*p_output));
    
    T_te=real(abs(t_te)^2*p_output/p_input); R_te=abs(r_te)^2; A_te=1-R_te-T_te;

    % results are stored in the T_R_A vector
    
    T_R_A(1)=T_te;
    T_R_A(2)=R_te;
    T_R_A(3)=A_te;
end