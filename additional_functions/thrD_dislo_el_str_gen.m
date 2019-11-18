function [epsi, rot] = thrD_dislo_el_str_gen(A,B,C,P,b,nu,h)

%written by F.Hofmann on 13/7/09
%function to carry out central difference numerical differentiation of the
%3D dislocation triangular loop Barnett solution to give elastic strain and
%lattice rotation

%modified by F. Hofmann on 5/11/18

%this script uses the "disp_dislo_tri_ABC_gen" function.

%inputs:
%A,B,C - three column vectors defining the nodes of the dislocation loop
%P - column vector with coordinates of the points at which the displacement
%is evaluated
%b - column vector with 3 burgers vector components
%h - scalar which gives the step used for the numerical differentiation

%generate the peturbed positions of P required:
P_xp = P + [h 0 0]';
P_xn = P - [h 0 0]';
P_yp = P + [0 h 0]';
P_yn = P - [0 h 0]';
P_zp = P + [0 0 h]';
P_zn = P - [0 0 h]';

P_all(1:3, 1+0*size(P,2):1*size(P,2)) = P; %displacement at point
P_all(1:3, 1+1*size(P,2):2*size(P,2)) = P + [ h 0 0]'*ones(1,size(P,2)); % xp
P_all(1:3, 1+2*size(P,2):3*size(P,2)) = P + [-h 0 0]'*ones(1,size(P,2)); % xn
P_all(1:3, 1+3*size(P,2):4*size(P,2)) = P + [0  h 0]'*ones(1,size(P,2)); % yp
P_all(1:3, 1+4*size(P,2):5*size(P,2)) = P + [0 -h 0]'*ones(1,size(P,2)); % yn
P_all(1:3, 1+5*size(P,2):6*size(P,2)) = P + [0 0  h]'*ones(1,size(P,2)); % zp
P_all(1:3, 1+6*size(P,2):7*size(P,2)) = P + [0 0 -h]'*ones(1,size(P,2)); % zn

%generate the perturbed displacement fields:
u_all = disp_dislo_tri_ABC_gen(A,B,C,P_all,b,nu);

%compute the numerical derivatives (comments behind give HLT notation of terms)
du_dx = (u_all(1:3, 1+1*size(P,2):2*size(P,2)) - u_all(1:3, 1+2*size(P,2):3*size(P,2)))./(2*h); 
du_dy = (u_all(1:3, 1+3*size(P,2):4*size(P,2)) - u_all(1:3, 1+4*size(P,2):5*size(P,2)))./(2*h);
du_dz = (u_all(1:3, 1+5*size(P,2):6*size(P,2)) - u_all(1:3, 1+6*size(P,2):7*size(P,2)))./(2*h);


%Assemble elastic strain tensor
%epsi = [dux_dx              0.5*(duy_dx+dux_dy) 0.5*(dux_dz+duz_dx) ;
%        0.5*(duy_dx+dux_dy) duy_dy              0.5*(duz_dy+duy_dz) ;
%        0.5*(dux_dz+duz_dx) 0.5*(duz_dy+duy_dz) duz_dz             ];

epsi(:,1,1) = du_dx(1,:);
epsi(:,1,2) = 0.5.*(du_dy(1,:)+du_dx(2,:));
epsi(:,1,3) = 0.5.*(du_dz(1,:)+du_dx(3,:));

epsi(:,2,1) = 0.5.*(du_dx(2,:)+du_dy(1,:));
epsi(:,2,2) = du_dy(2,:);
epsi(:,2,3) = 0.5.*(du_dz(2,:)+du_dy(3,:));

epsi(:,3,1) = 0.5.*(du_dx(3,:)+du_dz(1,:));
epsi(:,3,2) = 0.5.*(du_dy(3,:)+du_dz(2,:));
epsi(:,3,3) = du_dz(3,:);


%Assemble elastic strain tensor
%rot = [(0.5*(duz_dy-duy_dz));
%       (0.5*(dux_dz-duz_dx));
%       (0.5*(duy_dx-dux_dy))];
 
rot(:,1) = 0.5.*(du_dy(3,:)-du_dz(2,:));
rot(:,2) = 0.5.*(du_dz(1,:)-du_dx(3,:));
rot(:,3) = 0.5.*(du_dx(2,:)-du_dy(1,:));

end

        