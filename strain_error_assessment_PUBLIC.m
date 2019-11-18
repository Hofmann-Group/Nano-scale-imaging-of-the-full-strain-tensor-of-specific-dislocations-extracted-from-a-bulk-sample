%% script to look at errors between measured and simulated strains

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

%% load simulated strain fields
load('eps_rot_pred.mat')

%% make sample mask from reflection data
for iii = 1: size(Sam_red,1)
amp4d(:,:,:,iii) = Sam_red(iii,1).a_rr;
end
amp_max = squeeze(max(amp4d,[],4));
amp_all = squeeze(mean(amp4d,4));
amp_min = squeeze(min(amp4d,[],4));
clear amp4d
mask_all = amp_all; mask_all(mask_all<0.3) = 0; mask_all(mask_all>0) = 1;
mask_max = amp_max; mask_max(mask_max>0.3) = 1; mask_max(mask_max<1) = 0;

% Make from model predictions such that very large strains in pixels on the dislocation line in the model are excluded from the calculation
mask_mod = mask_max.*0+1;
for ii = 1:3
    for jj = 1:3
        pred_eps  = squeeze(eps_pred_3D(:,:,:,ii,jj));
        mask_mod(abs(pred_eps)>=0.01)=0; %exclude any pixels with strain greater than 1%.
    end
end

% Make mask that includes full crystal less a 30 nm thick surface layer....
contr = 6; %how many pixels to remove from the outside of mask_max...
[xx,yy,zz] = ndgrid(-contr:contr);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= contr;
mask_max_ero = imerode(mask_max,nhood).*mask_mod;


% Make mask that includes only material outside dislocaitions 30 nm radius around dislocations
mask_min = amp_min; mask_min(mask_min<0.1) = 0; mask_min(mask_min>0) = 1;
contr = 6; %how many pixels to remove from the outside of mask_max...
[xx,yy,zz] = ndgrid(-contr:contr);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= contr;
mask_no_dislo = imerode(mask_min,nhood).*mask_mod;

% Make a mask that only contains strain within a 30 nm radius of dislocation lines. 
mask_min = amp_min; mask_min(mask_min<0.1) = 0; mask_min(mask_min>0) = 1;
mask_min_inv = 1-mask_min;
contr = 6; %how many pixels to remove from the outside of mask_max...
[xx,yy,zz] = ndgrid(-contr:contr);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= contr;
mask_max_ero = imerode(mask_max,nhood);
mask_skinny_dislo = mask_max_ero.*mask_min_inv;
mask_skinny_dislo(190:220, 105:135, :) = 0;
mask_skinny_dislo(60:80, 70:85, :) = 0;
mask_skinny_dislo(:, 150:165, 60:80) = 0;
mask_skinny_dislo(:, 120:130, 70:80) = 0;
mask_fat_dislo = ((1-imerode(1-mask_skinny_dislo,nhood)).*mask_max_ero).*mask_mod;


%% compare strains - exp vs predicted strains in all voxels within mask_max
% visualise what the region within which strains are compared looks like...
figure;
amp2plot = 1-mask_fat_dislo;
pat = patch(isosurface(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9,amp2plot,0.3));
isonormals(Xgrid, Ygrid, Zgrid,amp2plot,pat);
set(pat,'FaceColor','red','EdgeColor','none');
alpha(pat,0.1); hold on;
%xlabel('X'), ylabel('Y'), zlabel('Z')
pat = patch(isosurface(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9,amp_max,0.3));
isonormals(Xgrid, Ygrid, Zgrid,amp2plot,pat);
set(pat,'FaceColor','yellow','EdgeColor','none');
alpha(pat,0.1); hold on;
daspect([1,1,1]); view([0 0 1]); axis equal vis3d xy on;
camlight;  hold on; camup([0 1 0]); view(-180, -55);
%plot coordiantes:
quiver3(-490,-490,-490,500./0.9,0,0, 'r', 'LineWidth', 1.5); quiver3(-490,-490,-490,0,500./0.9,0, 'g', 'LineWidth', 1.5); quiver3(-490,-490,-490,0,0,500./0.9, 'b', 'LineWidth', 1.5);

%% Pick for which volume to compare strains...
% measure average difference between predicted and measured strains for different crystal regions:
%filt_mask = mask_max_ero; % whole crystal including dislocation, less a 30 nm thick near-surface layer. 
%filt_mask = mask_no_dislo; % material of the crystal where there are no dislocations & exclusing surface layer. excluded 30 nm surface layer and 30 nm radius around dislocaitons. 
filt_mask = mask_fat_dislo; %mask such that only the strain within pipes of ~30 nm radius around the dislocations is counted. 

for ii = 1:3
    for jj = 1:3        
        exp_eps  = squeeze(Sam_red(1).eps_waf(:,:,:,ii,jj));
        pred_eps = squeeze(eps_pred_3D(:,:,:,ii,jj));
        
        A = abs(pred_eps(filt_mask==1));
        eps_mean(ii,jj) = mean(A(:));
        
        B = abs(pred_eps(filt_mask==1)-exp_eps(filt_mask==1));
        eps_diff(ii,jj) = mean(B(:));
        

    end
    
    exp_rot  = squeeze(Sam_red(1).rot_waf(:,:,:,ii));
    pred_rot = squeeze(rot_pred_3D(:,:,:,ii));
    
    C = abs(pred_rot(filt_mask==1));
    rot_mean(ii,1) = mean(C(:));
    
    D = abs(pred_rot(filt_mask==1)-exp_rot(filt_mask==1));
    rot_diff(ii,1) = mean(D(:));
    
end
format short e
eps_mean
eps_diff
rot_mean
rot_diff


%% Comparison of strain around circular paths around dislocations:

% plot a slice in the YZ plane through the crystal...
p_sl = 105 % pick which slice to consider. 
figure
subplot(1,2,1)
pcolor(squeeze(Zgrid(:,p_sl,:).*10^9), squeeze(Ygrid(:,p_sl,:).*10^9), squeeze(Sam_red(1,1).eps_waf(:, p_sl,:,1,1)));
axis equal tight
shading flat
caxis([-0.002 0.002])
subplot(1,2,2)
pcolor(squeeze(Zgrid(:,p_sl,:).*10^9), squeeze(Ygrid(:,p_sl,:).*10^9), squeeze(eps_pred_3D(:,p_sl,:,1,1)));
axis equal tight
shading flat
caxis([-0.002 0.002])
suptitle(['xposition: ', num2str(Xgrid(1,p_sl,1).*10^9), 'nm'])
hold on


ang = [0:10:360]./360.*2*pi; % angular positions for angular plot
ang_deg = [0:10:360];
r = ang.*0 + 30; % radial position for angular plot
[p1, p2] = pol2cart(ang,r);

% dislocations for which the comparison of strain around a circular path was done...
% x position: -77.5, dislocation positions: 
dpos(:,2) = [-77.5; 172.5; -52.5]; % dislocation 2 position 
dpos(:,4) = [-77.5; -322.5; 197.5]; % dislocation 4 position 
dpos(:,5) = [-77.5; -267.5; 102.5]; % dislocation 5 position 

plot(p1(1:25)-52.5, p2(1:25)+172.5, 'k')

% Pick for which dislocation the strain variation along a circular pathe should be plotted...
dislo2plot = 5;
for ii = 1:3
    for jj = 1:3
        eps_cont_exp(:,ii,jj) = interp3(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9, squeeze(Sam_red(1,1).eps_waf(:,:,:,ii,jj)), ...
            p1.*0+dpos(1,dislo2plot), p2+dpos(2,dislo2plot), p1+dpos(3,dislo2plot));
        
        eps_cont_sim(:,ii,jj) = interp3(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9, squeeze(eps_pred_3D(:,:,:,ii,jj)), ...
            p1.*0+dpos(1,dislo2plot), p2+dpos(2,dislo2plot), p1+dpos(3,dislo2plot));
        
        eps_cont_diff(ii,jj) = mean(squeeze(abs(eps_cont_exp(:,ii,jj)-eps_cont_sim(:,ii,jj))));
        
    end
    rot_cont_exp(:,ii) = interp3(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9, squeeze(Sam_red(1,1).rot_waf(:,:,:,ii)), ...
            p1.*0+dpos(1,dislo2plot), p2+dpos(2,dislo2plot), p1+dpos(3,dislo2plot));
        
    rot_cont_sim(:,ii) = interp3(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9, squeeze(rot_pred_3D(:,:,:,ii)), ...
            p1.*0+dpos(1,dislo2plot), p2+dpos(2,dislo2plot), p1+dpos(3,dislo2plot));
        
    rot_cont_diff(ii,1) = mean(squeeze(abs(rot_cont_exp(:,ii)-rot_cont_sim(:,ii))));
end

eps_cont_diff
rot_cont_diff

figure; 
axlim = [0, 360, -1.6e-3, 1.6e-3]
lw = 1.5;
subplot(3,3,1);
plot(ang_deg, squeeze(eps_cont_exp(:,1,1)), 'r', 'LineWidth', lw); hold on; plot(ang_deg, squeeze(eps_cont_sim(:,1,1)), 'b', 'LineWidth', lw); axis(axlim); title('\epsilon_{xx}')
subplot(3,3,2);
plot(ang_deg, squeeze(eps_cont_exp(:,1,2)), 'r', 'LineWidth', lw); hold on; plot(ang_deg, squeeze(eps_cont_sim(:,1,2)), 'b', 'LineWidth', lw); axis(axlim); title('\epsilon_{xy}')
subplot(3,3,3);
plot(ang_deg, squeeze(eps_cont_exp(:,1,3)), 'r', 'LineWidth', lw); hold on; plot(ang_deg, squeeze(eps_cont_sim(:,1,3)), 'b', 'LineWidth', lw); axis(axlim); title('\epsilon_{xz}')
subplot(3,3,4);
plot(ang_deg, squeeze(rot_cont_exp(:,3)), 'r', 'LineWidth', lw); hold on; plot(ang_deg, squeeze(rot_cont_sim(:,3)), 'b', 'LineWidth', lw); axis(axlim); title('\omega_{z}')
subplot(3,3,5);
plot(ang_deg, squeeze(eps_cont_exp(:,2,2)), 'r', 'LineWidth', lw); hold on; plot(ang_deg, squeeze(eps_cont_sim(:,2,2)), 'b', 'LineWidth', lw); axis(axlim); title('\epsilon_{yy}')
subplot(3,3,6);
plot(ang_deg, squeeze(eps_cont_exp(:,2,3)), 'r', 'LineWidth', lw); hold on; plot(ang_deg, squeeze(eps_cont_sim(:,2,3)), 'b', 'LineWidth', lw); axis(axlim); title('\epsilon_{yz}')
subplot(3,3,7);
plot(ang_deg, squeeze(rot_cont_exp(:,2)), 'r', 'LineWidth', lw); hold on; plot(ang_deg, squeeze(rot_cont_sim(:,2)), 'b', 'LineWidth', lw); axis(axlim); title('\omega_{y}')
subplot(3,3,8);
plot(ang_deg, squeeze(rot_cont_exp(:,1)), 'r', 'LineWidth', lw); hold on; plot(ang_deg, squeeze(rot_cont_sim(:,1)), 'b', 'LineWidth', lw); axis(axlim); title('\omega_{x}')
subplot(3,3,9);
plot(ang_deg, squeeze(eps_cont_exp(:,3,3)), 'r', 'LineWidth', lw); hold on; plot(ang_deg, squeeze(eps_cont_sim(:,3,3)), 'b', 'LineWidth', lw); axis(axlim); title('\epsilon_{zz}')
suptitle(['Contour around dislo ', num2str(dislo2plot), ', radius ', num2str(r(1)), ' nm'])

