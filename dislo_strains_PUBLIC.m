clear all
close all

%% Material parameters
% material and simulation proprerties
C = [0 -5E-5 0]'; %not sure if this is the best of all places to define C. C is the "3rd" point in every dislocation triangle...
nu = 0.28; %tungsten poisson ratio. 
YM = 410*10^9; %Young's Modulus in Pa
a0 = 3.165*10^-10; %tungsten lattice constant in m
% for W properties cite Bolef, Di and J. de Klerk, Elastic constants of single crystal Mo and W between 77 and 500K, JAP, 33, 2311-2314, 1962
% and Featherstone, F.H. and J.R. Neighbours, Elastic constants of tantalum, tungsten and molybdenum, Phys Rev. 130, 1324-1333, 1963

%% load experimental data 
load('recD_sam_str_rot.mat')


%% load node data from Nicks dislocation finding:

%plot the dislocation positions...
load('segmentation190123.mat'); % variable in this is called tmp1
tmp1 = segmentation;
nn = 1; % reduced resolution along dislocation line...

% extract dislocation node coordinates and add pints outside object where needed. 
ddd = tmp1(1,4); d1xyzt = cell2mat([ddd{:}]);  d1_vect = (d1xyzt(end,:) - d1xyzt(1,:)).*10; 
                        d1xyz = [d1xyzt(1,:)-d1_vect ; d1xyzt(1:nn:end-1,:); d1xyzt(end,:) ;d1xyzt(end,:)+d1_vect];
                        
ddd = tmp1(2,4); d2xyzt = cell2mat([ddd{:}]);  d2_vect = (d2xyzt(end,:) - d2xyzt(1,:)).*10; 
                        d2xyz = [d2xyzt(1,:)-d2_vect ; d2xyzt(1:nn:end-1,:); d2xyzt(end,:) ; d2xyzt(end,:)+d2_vect];
                        
ddd = tmp1(3,4); d3xyzt = cell2mat([ddd{:}]);  d3xyz = [d3xyzt(1,:) + [0 -1000 0]; d3xyzt(1:nn:end-1,:) ; d3xyzt(end,:)];

ddd = tmp1(4,4); d4xyzt = cell2mat([ddd{:}]);  d4xyz = [d4xyzt(1:nn:end-1,:); d4xyzt(end,:); [0 -1000 0] + d4xyzt(end,:)];

ddd = tmp1(5,4); d5xyzt = cell2mat([ddd{:}]);  d5xyz = [d5xyzt(1:nn:end-1,:); d5xyzt(end,:); [0 -1000 0] + d5xyzt(end,:)];

% dislo 3, 4 and 5 don't quite meet in the same place - fix this...
d3xyz(end,:) = d4xyz(1,:); 

% transform dislocation positions from pixel to sample coordiante positions: 
mid_pix = [240./2 240./2 240./2]; % middle pixel position in x, y, z directions 
pix_size = 5*10^-9; %voxel size in m. 
d1xyz_pos = (d1xyz - ones(size(d1xyz,1),1)*mid_pix).*pix_size;
d2xyz_pos = (d2xyz - ones(size(d2xyz,1),1)*mid_pix).*pix_size;
d3xyz_pos = (d3xyz - ones(size(d3xyz,1),1)*mid_pix).*pix_size;
d4xyz_pos = (d4xyz - ones(size(d4xyz,1),1)*mid_pix).*pix_size;
d5xyz_pos = (d5xyz - ones(size(d5xyz,1),1)*mid_pix).*pix_size;

% make sample morphology from reflection amplitude data
for iii = 1: size(Sam_red,1)
amp4d(:,:,:,iii) = Sam_red(iii,1).a_rr;
end
amp_max = squeeze(max(amp4d,[],4));
amp_all = squeeze(mean(amp4d,4));
amp_min = squeeze(min(amp4d,[],4));
clear amp4d
mask_all = amp_all; mask_all(mask_all<0.3) = 0; mask_all(mask_all>0) = 1;

% make rn list of dislocation nodes: 
rn = [d1xyz_pos; d2xyz_pos; d3xyz_pos; d4xyz_pos; d5xyz_pos]+0.1*pix_size; %offset every position by 0.1 od a pixel for numerical stability
rn(:,4) = 7; %make all the nodes fixed. 

% work out dislocation burgers vectors
% compute crystal basis vectors: 
% reflection order:
%10-1,
%-10-1
%1-10
%0-1-1
%110
%01-1
%using the q directions work out the directions of crystal 100, 010 and 001 axes. Call these n100, n010, n001 respectively. 
n100 = Sam_red(1,1).q_sam - Sam_red(2,1).q_sam;  n100  = n100 ./norm(n100);
n100a = Sam_red(3,1).q_sam + Sam_red(5,1).q_sam; n100a = n100a./norm(n100a); %compute as a check...
n010 = Sam_red(5,1).q_sam - Sam_red(3,1).q_sam;  n010  = n010 ./norm(n010);
n010a = Sam_red(6,1).q_sam - Sam_red(4,1).q_sam; n010a = n010a./norm(n010a); %compute as a check...
n001 = -Sam_red(1,1).q_sam - Sam_red(2,1).q_sam;  n001  = n001 ./norm(n001);
n001a = -Sam_red(4,1).q_sam - Sam_red(6,1).q_sam; n001a = n001a./norm(n001a); %compute as a check...
n001b = cross(n100, n010); % just to check this 001 cross 010 gives 001...

bv1_dir = -0.5 * [1 1 1]';   bv1 = (bv1_dir(1,1)* n100 + bv1_dir(2,1)* n010 + bv1_dir(3,1)* n001)* a0; bv1n = bv1./norm(bv1);
bv2_dir = 0.5 * [-1 1 1]';  bv2 = (bv2_dir(1,1)* n100 + bv2_dir(2,1)* n010 + bv2_dir(3,1)* n001)* a0; bv2n = bv2./norm(bv2);
bv3_dir = -[1 0 0]';         bv3 = (bv3_dir(1,1)* n100 + bv3_dir(2,1)* n010 + bv3_dir(3,1)* n001)* a0; bv3n = bv3./norm(bv3);
bv4_dir = -0.5 * [1 1 1]';   bv4 = (bv4_dir(1,1)* n100 + bv4_dir(2,1)* n010 + bv4_dir(3,1)* n001)* a0; bv4n = bv4./norm(bv4);
bv5_dir = -0.5 * [1 -1 -1]'; bv5 = (bv5_dir(1,1)* n100 + bv5_dir(2,1)* n010 + bv5_dir(3,1)* n001)* a0; bv5n = bv5./norm(bv5);

% work out the segments for the dislocations
Nd1 = size(d1xyz_pos,1);
Nd2 = size(d2xyz_pos,1);
Nd3 = size(d3xyz_pos,1);
Nd4 = size(d4xyz_pos,1);
Nd5 = size(d5xyz_pos,1);

links1(:, 1:2) = [(1:Nd1-1)', (2:Nd1)'];                links1(:,3:5) = ones(Nd1-1,1)*bv1';
links2(:, 1:2) = [(1:Nd2-1)', (2:Nd2)']+Nd1;            links2(:,3:5) = ones(Nd2-1,1)*bv2';
links3(:, 1:2) = [(1:Nd3-1)', (2:Nd3)']+Nd1+Nd2;        links3(:,3:5) = ones(Nd3-1,1)*bv3';
links4(:, 1:2) = [(1:Nd4-1)', (2:Nd4)']+Nd1+Nd2+Nd3;    links4(:,3:5) = ones(Nd4-1,1)*bv4';
links5(:, 1:2) = [(1:Nd5-1)', (2:Nd5)']+Nd1+Nd2+Nd3+Nd4;links5(:,3:5) = ones(Nd5-1,1)*bv5';

links = [links1 ; links2 ; links3; links4; links5];


% Plot dislocation positions overlaid on sample morphology...
figure;
amp2plot = amp_all;
pat = patch(isosurface(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9,amp2plot,0.2));
isonormals(Xgrid, Ygrid, Zgrid,amp2plot,pat);
set(pat,'FaceColor','yellow','EdgeColor','none');
alpha(pat,0.1); hold on;
daspect([1,1,1]); view([0 0 1]); axis equal vis3d xy on;
camlight;  hold on; camup([0 1 0]); view(-180, -55);
%plot coordiantes:
quiver3(-490,-490,-490,500./0.9,0,0, 'r', 'LineWidth', 1.5); quiver3(-490,-490,-490,0,500./0.9,0, 'g', 'LineWidth', 1.5); quiver3(-490,-490,-490,0,0,500./0.9, 'b', 'LineWidth', 1.5);
%plot dislocations: 
plot3((d1xyz_pos(2:end-1,1))*1E9, (d1xyz_pos(2:end-1,2))*1E9, (d1xyz_pos(2:end-1,3))*1E9, 'm-', 'LineWidth', 3); hold on
plot3((d2xyz_pos(2:end-1,1))*1E9, (d2xyz_pos(2:end-1,2))*1E9, (d2xyz_pos(2:end-1,3))*1E9, 'c-', 'LineWidth', 3); hold on
plot3((d3xyz_pos(2:end,1))*1E9, (d3xyz_pos(2:end,2))*1E9, (d3xyz_pos(2:end,3))*1E9, 'r-', 'LineWidth', 3); hold on
plot3((d4xyz_pos(1:end-1,1))*1E9, (d4xyz_pos(1:end-1,2))*1E9, (d4xyz_pos(1:end-1,3))*1E9, 'g-', 'LineWidth', 3); hold on
plot3((d5xyz_pos(1:end-1,1))*1E9, (d5xyz_pos(1:end-1,2))*1E9, (d5xyz_pos(1:end-1,3))*1E9, 'b-', 'LineWidth', 3); hold on
%plot burgers vectors:
bvpl = 120;
d1pp = (d1xyz_pos(round(end/2),:))*1E9;
quiver3(d1pp(1,1),d1pp(1,2),d1pp(1,3),bv1n(1,1)*bvpl,bv1n(2,1)*bvpl,bv1n(3,1)*bvpl, 'k', 'LineWidth', 1.5, 'MaxHeadSize',15);
d2pp = (d2xyz_pos(round(end*0.75),:))*1E9;
quiver3(d2pp(1,1),d2pp(1,2),d2pp(1,3),bv2n(1,1)*bvpl,bv2n(2,1)*bvpl,bv2n(3,1)*bvpl, 'k', 'LineWidth', 1.5, 'MaxHeadSize',15);
d3pp = (d3xyz_pos(round(end/2),:))*1E9;
quiver3(d3pp(1,1),d3pp(1,2),d3pp(1,3),bv3n(1,1)*bvpl,bv3n(2,1)*bvpl,bv3n(3,1)*bvpl, 'k', 'LineWidth', 1.5, 'MaxHeadSize',15);
d4pp = (d4xyz_pos(round(end/2),:))*1E9;
quiver3(d4pp(1,1),d4pp(1,2),d4pp(1,3),bv4n(1,1)*bvpl,bv4n(2,1)*bvpl,bv4n(3,1)*bvpl, 'k', 'LineWidth', 1.5, 'MaxHeadSize',15);
d5pp = (d5xyz_pos(round(end/2),:))*1E9;
quiver3(d5pp(1,1),d5pp(1,2),d5pp(1,3),bv5n(1,1)*bvpl,bv5n(2,1)*bvpl,bv5n(3,1)*bvpl, 'k', 'LineWidth', 1.5, 'MaxHeadSize',15);
% add labels to the axes...
%xlabel('X'); ylabel('Y'); zlabel('Z')
axis([-650 650 -650 650 -650 650])
axis off
% pp = fill3([-77.5 -77.5 -77.5 -77.5],[-450 490 490 -450],[420 420 -420 -420], 'g') % option to plot a section plane...
% alpha(pp,0.2); 

%% calculate the anticipated dislocation strains using linear elasticity and the Barnett dislocation triangle solution
% use Xgrid, Ygrid and Zgrid that's also used for plotting of experimental results....

% 3D to save time only evaluate at points inside the crystal...

mask_max = amp_max;
mask_max(mask_max>0.3) = 1; mask_max(mask_max<1) = 0;
Xgr_sm = Xgrid(mask_max==1);
Ygr_sm = Ygrid(mask_max==1);
Zgr_sm = Zgrid(mask_max==1);

% generate the matrices in which the results will be stored...
eps_new = zeros(size(Xgr_sm(:), 1), 3, 3);
rot_new = zeros(size(Xgr_sm(:), 1), 3);

% calculate dislocation strain fields
ccc_max = size(links,1);
ccc = 0;
h = 0.005*pix_size; % step for numerical differentiation is one hundredth of the distance between points. 

for linkno = 1 : size(links,1);
    %counter
    format shortG
    ccc = ccc+1; ccc/ccc_max
    %define geometry
    A = double(rn(links(linkno,1), 1:3)');
    B = double(rn(links(linkno,2), 1:3)');

    %define burgers vector
    bv = double(links(linkno,3:5)');
    
    % new calculation approach
    % define empty matrix for strains of this element...
   
    tic
    P = [Xgr_sm(:), Ygr_sm(:), Zgr_sm(:)]';
    [eps_link_new, rot_link_new] = thrD_dislo_el_str_gen(A,B,C,P,bv,nu,h);
    
    eps_new = eps_new + eps_link_new;
    rot_new = rot_new + rot_link_new;
    toc
    
end

% 3D sort strain into correct matricies: 
eps_pred_3D = zeros(size(Xgrid,1), size(Xgrid,2), size(Xgrid,3), 3, 3);
rot_pred_3D = zeros(size(Xgrid,1), size(Xgrid,2), size(Xgrid,3), 3);

for iii = 1:3
    for jjj = 1:3
        eps_tmp = zeros(size(Xgrid,1), size(Xgrid,2), size(Xgrid,3));
        eps_tmp(mask_max==1) = eps_new(:,iii,jjj);
        eps_pred_3D(:,:,:,iii,jjj) = eps_tmp;
        clear eps_tmp
    end
    rot_tmp = zeros(size(Xgrid,1), size(Xgrid,2), size(Xgrid,3));
    rot_tmp(mask_max==1) = rot_new(:,iii);
    rot_pred_3D(:,:,:,iii) = rot_tmp;
    clear rot_tmp
end

% save the outputs of the calculation
save('eps_rot_pred.mat', 'Xgrid', 'Ygrid', 'Zgrid', 'eps_pred_3D', 'rot_pred_3D')

%% video of lattice strains and rotations predicted from dislocation model...

v = VideoWriter('3D_model_strain_rot_tensor_smaller_scale.mp4', 'MPEG-4');
v.Quality = 95;
open(v)
calim = [-0.0008 0.0008];

figure;
%colormap(RdYeBlue);
for aaa = 60:184 % size(Sam(iii,1).oph_rr, 2) This gives correct range for 5 nm real-space pixel size... 
P(1) = subplot(3,3,1); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_pred_3D(:,aaa,:,1,1))); view([0 0 1])
%str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(Sam_red(1,1).eps_waf(:,aaa,:,1,1))); view([0 0 1])
shading flat; axis equal; axis([-480 480 -480 480]); %use axis([-500 500 -500 500]) to make video...
caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{xx}')
axis off

P(2) = subplot(3,3,2); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_pred_3D(:,aaa,:,1,2))); view([0 0 1])
%str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(Sam_red(1,1).eps_waf(:,aaa,:,1,2))); view([0 0 1])
shading flat; axis equal; axis([-480 480 -480 480]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{xy}')
axis off

P(3) = subplot(3,3,3); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_pred_3D(:,aaa,:,1,3))); view([0 0 1])
%str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(Sam_red(1,1).eps_waf(:,aaa,:,1,3))); view([0 0 1])
shading flat; axis equal; axis([-480 480 -480 480]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{xz}')
axis off

P(5) = subplot(3,3,5); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_pred_3D(:,aaa,:,2,2))); view([0 0 1])
%str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(Sam_red(1,1).eps_waf(:,aaa,:,2,2))); view([0 0 1])
shading flat; axis equal; axis([-480 480 -480 480]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{yy}')
axis off

P(6) = subplot(3,3,6); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_pred_3D(:,aaa,:,2,3))); view([0 0 1])
%str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(Sam_red(1,1).eps_waf(:,aaa,:,2,3))); view([0 0 1])
shading flat; axis equal; axis([-480 480 -480 480]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{yz}')
axis off

P(9) = subplot(3,3,9); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_pred_3D(:,aaa,:,3,3))); view([0 0 1])
%str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(Sam_red(1,1).eps_waf(:,aaa,:,3,3))); view([0 0 1])
shading flat; axis equal; axis([-480 480 -480 480]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{zz}')
axis off

P(4) = subplot(3,3,4); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(rot_pred_3D(:,aaa,:,3))); view([0 0 1])
%str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(Sam_red(1,1).rot_waf(:,aaa,:,3))); view([0 0 1])
shading flat; axis equal; axis([-480 480 -480 480]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\omega_{z}')
axis off

P(7) = subplot(3,3,7); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(rot_pred_3D(:,aaa,:,2))); view([0 0 1])
%str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(Sam_red(1,1).rot_waf(:,aaa,:,2))); view([0 0 1])
shading flat; axis equal; axis([-480 480 -480 480]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\omega_{y}')
axis off

P(8) = subplot(3,3,8); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(rot_pred_3D(:,aaa,:,1))); view([0 0 1])
%str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(Sam_red(1,1).rot_waf(:,aaa,:,1))); view([0 0 1])
shading flat; axis equal; axis([-480 480 -480 480]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\omega_{x}')
hold on
patch([-250 250 250 -250],[-430 -430 -450 -450], [0 0 0])
axis off


h=colorbar;
set(h, 'Position', [.9 .18 .03 .68])
pos1 = get(P(1), 'Position');
for i=1:9
    pos=get(P(i), 'Position');
    set(P(i), 'Position', [0.7*pos(1)-0.5*pos1(1) 1.1*pos(2)-0*pos1(2) 1.4*pos(3) 1.4*pos(4)]);
end


suptitle(['Predicted strain and rotation. Xpos: ', num2str(Xgrid(1,aaa,1)*10^9), ' nm'])

% make video 
frame = print('-RGBImage', '-r300');
writeVideo(v,frame) ;
end
close(v)
