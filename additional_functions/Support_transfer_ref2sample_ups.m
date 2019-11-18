function [Dat_out, SXgrid, SYgrid, SZgrid, q_sam] = Support_transfer_ref2sample_ups(Dat_in, M, Sam, options)

% Script to map support or reconstruction of reflection from detector conjugated space to the sample space. 
% F.Hofmann Oxford August 2018 

%% Inputs
%% Reflection details & data in detector conjugated frame. 
% Dat  - 3D data file in detector conjugated space as returned by Jesse's phasing program. 

% detector details
% M.gam_spec - gamma detector angles as given in spec file
% M.del_spec - delta detector angle as given in spec file
% M.dd - sample to detector distance in m
% M.pix - detector pixel size in m - this is effective, so larger than actual if binnig is used...

% sample details: 
% M.the_spec - sample stack theta angle as given in spec file...
% M.chi_spec - sample stack chi angle as given in spec file...
% M.phi_spec - sample stack phi angle as given in spec file...
% NOTE: to get everything in lab Coordinates as in Jesse's code set: M.the_spec = 0; M.chi_spec = 90; M.phi_spec = 0; 
% NOTE: to get reconstruction in sample coordinates that are aligned with lab coordiantes when the sample is not tilted give the angles noted in the spec file for M.the_spec, M.chi_spec and M.phi_spec

% experiment/scan details: 
% M.inc_the_spec - scanning step increment in theta in degrees. direction same as in spec. 
% M.lam - wavelength in m.

%% Output details:
% Sam.N1 -  Size of the output image in sample coordinate space in number of pixels. Must be integer and even. Y direction
% Sam.N2 -  Size of the output image in sample coordinate space in number of pixels. Must be integer and even. X direction. 
% Sam.N3 -  Size of the output image in sample coordinate space in number of pixels. Must be integer and even. Z direction
% Sam.pix - pixel size in sample coordinates in m
% Sam.ups - factor by which diffraction pattern is zero-padded to give finer resoltuion in real-space. 

%% Options: 
% 'none' - just transfers from det to sample space
% 'centred' - centre the centre of mass of the amplitude of the input

%% Check the input
if nargin == 3
    options = 'none';
end

%% Assume 34IDE coordinate frame...
% flip detector data so it's the right way around for convention used (q1 up, q2 along x axis)
%M.dat = flip(flip(M.dat,2),1);
Dat_in = flip(flip(Dat_in,2),1);

% do Fourier padding...
pads = round(size(Dat_in)./2 * (Sam.ups-1));

Dat_in_ft = fftshift(fftn(Dat_in));
Dat_in_ft_pad = padarray(Dat_in_ft, pads, 0, 'both'); % padarray(A, [1,2,2], 0, 'both')
Dat = ifftn(ifftshift(Dat_in_ft_pad));

% get size of the reflection matrix in detector conjugated space.
[M.N1, M.N2, M.N3] = size(Dat);

% make grids of N1, N2 and N3 integers that address the pixels in the 3D volume
% N grid for the master reflection...
[M.N1grid, M.N2grid, M.N3grid] = ndgrid(-(M.N1-1)/2:(M.N1-1)/2, ...
                                          -(M.N2-1)/2:(M.N2-1)/2, ...
                                          -(M.N3-1)/2:(M.N3-1)/2);

% N grid for the sample space                                     
[Sam.N1grid, Sam.N2grid, Sam.N3grid] = ndgrid(-(Sam.N1-1)/2:(Sam.N1-1)/2, ...
                                          -(Sam.N2-1)/2:(Sam.N2-1)/2, ...
                                          -(Sam.N3-1)/2:(Sam.N3-1)/2);
                                      
% meshgrid of coordinates in the sample space:                                       
SXgrid = Sam.pix .* Sam.N2grid;                        
SYgrid = Sam.pix .* Sam.N1grid;
SZgrid = Sam.pix .* Sam.N3grid;

% convert angles to right handed convention used in F.Hofmann J. Synch Rad 2017 paper. Also convert to radians
M.phi = M.phi_spec*pi/180; 
M.chi = (90-M.chi_spec)*pi/180; 
M.the = M.the_spec*pi/180; 
M.gam = -M.gam_spec*pi/180; 
M.del = M.del_spec*pi/180; 
M.inc_the = M.inc_the_spec*pi/180;

%% make coordinates...
% compute q and q1, q2, q3 reciprocal space vectors in lab coordinates
s0_lab = 2*pi/M.lam*[0 0 1]';
M.s_lab = roty(M.del)*rotx(M.gam)*s0_lab;
M.q_lab = M.s_lab - s0_lab; % q vector for centre of detector.

M.q1_n_lab = roty(M.del)*rotx(M.gam)*[0 1 0]'; % unit vector along q1 direction in lab coordiantes [0 1 0]
M.dq1 = 2*pi/M.lam*M.pix/M.dd;
M.q2_n_lab = roty(M.del)*rotx(M.gam)*[1 0 0]'; % [1 0 0]
M.dq2 = 2*pi/M.lam*M.pix/M.dd;

q3_temp = (roty(M.inc_the)*M.q_lab)-M.q_lab;
M.dq3 = norm(q3_temp);
M.q3_n_lab = q3_temp./M.dq3;
clear q3_temp

% make r1, r2, r3 and dr1, dr2, dr3 in lab coordiantes. This follows the treatment in Berenguer et al. PRB 88, 144101 (2013).
M.Wq = dot(cross(M.q1_n_lab, M.q2_n_lab), M.q3_n_lab)* M.dq1 * M.N1 * M.dq2  * M.N2 * M.dq3 * M.N3;
M.r1_lab = 2*pi*cross(M.N2*M.dq2*M.q2_n_lab, M.N3*M.dq3*M.q3_n_lab)./M.Wq;
M.r2_lab = 2*pi*cross(M.N3*M.dq3*M.q3_n_lab, M.N1*M.dq1*M.q1_n_lab)./M.Wq;
M.r3_lab = 2*pi*cross(M.N1*M.dq1*M.q1_n_lab, M.N2*M.dq2*M.q2_n_lab)./M.Wq;

% rotate Master r1, r2, r3 into sample coordinates. Sample xyz coordinates
% are aligned with lab coordiantes when the sample stage angles are: M.the_spec = 0; M.chi_spec = 90; M.phi_spec = 0; 
M.Rsam = roty(M.the)*rotz(M.chi)*rotx(M.phi); % rotation matrix to rotate a vector in sample coordiantes into lab coordinates:
M.r1_sam = M.Rsam'*M.r1_lab;
M.r2_sam = M.Rsam'*M.r2_lab;
M.r3_sam = M.Rsam'*M.r3_lab;

% rotate q vector for detector centre into sample coordinates: 
q_sam = M.Rsam'*M.q_lab;

%% map from conjugated space to sample space:
% interpolate the master reconstruction from the master reflection
% conjugated detector frame into the slave reflection detector conjugated
% frame: 

% make the r_m matrix:
M.r_m = [M.r1_sam, M.r2_sam, M.r3_sam];

% make T matrix: 
T = inv(M.r_m);

% centre if required:
com = centerOfMass(abs(Dat))-(size(Dat)+1)/2

% map coordinates of sample space points into fractional indicies of master
% reflection space:
Sam.nm1 = T(1,1)*SXgrid + T(1,2)*SYgrid + T(1,3)*SZgrid;
Sam.nm2 = T(2,1)*SXgrid + T(2,2)*SYgrid + T(2,3)*SZgrid;
Sam.nm3 = T(3,1)*SXgrid + T(3,2)*SYgrid + T(3,3)*SZgrid;

if strcmp(options, 'none')
    % interpolate reflection data in the detector conjugated frame to get data in the orthogonal sample coordiante frame. 
    Sam.dat = interp3(M.N2grid, M.N1grid, M.N3grid, Dat, Sam.nm2, Sam.nm1, Sam.nm3, 'linear',  0); % make any values outside master data zero.
elseif strcmp(options, 'centred')
    % interpolate reflection data in the detector conjugated frame to get data in the orthogonal sample coordiante frame. 
    Sam.dat = interp3(M.N2grid-com(2), M.N1grid-com(1), M.N3grid-com(3), Dat, Sam.nm2, Sam.nm1, Sam.nm3, 'linear',  0); % make any values outside master data zero.
end


% compute conjugated reflection so that the output is consistent with the
% mapping to lab coordinates done in Jesse Clark's code. 
Dat_out = conj_reflect(Sam.dat)*(Sam.ups^3); % multiply by up sampling factor so that the magnitude is the same before and after upsampling...

end

% little function to work out the centre of mass of an array
function com = centerOfMass(A)
As = size(A);
Nd = ndims(A);
M = sum(A(:));
com = zeros(1,Nd);
if M==0
    com = [];
else
    for ii = 1:Nd
        shp = ones(1,Nd);
        shp(ii) = As(ii);
        rep = As;
        rep(ii) = 1;
        ind = repmat(reshape(1:As(ii),shp),rep);
        com(ii) = sum(ind(:).*A(:))./M;
    end
end
end