% This code is associated with the following publication: 
% "Nano-scale imaging of the full strain tensor of specific dislocations extracted from a bulk sample"
% Felix Hofmann, Nicholas W. Phillips, Suchandrima Das, Phani Karamched, Gareth M. Hughes, James O. Douglas, Wonsuk Cha, Wenjun Liu
% Physical Review Materials 2019
% https://arxiv.org/abs/1903.04079

% Felix Hofmann, University of Oxford, 2019


% analysis of multi-reflection data from dislocation 7 sample...

clear all
close all

fact = 1; % this is how to reduce frame size and keep correct scaling...
ss1 = 1:256;
ss2 = 1:256;
ss3 = 1:256;

%% Reflection details & data in detector conjugated frame. 
% the files loaded here are the output of the Phasing codes written by Jesse Clark. 
%%%%%  10-1 reflection of dislo 7 - scans 1064 to 1121: 
R(1,1).ref = [1 0 -1]';
Amp = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_10-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_10-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW-AMP.mat', '-mat')));
Ph  = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_10-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_10-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW-PH.mat', '-mat')));
R(1,1).dat = Amp(ss1, ss2, ss3).* exp(1i * Ph(ss1, ss2, ss3));
% detector angles as given in spec file
R(1,1).gam_spec = 23.3615;
R(1,1).del_spec = 22.64375;
% sample stack angles as given in spec file
R(1,1).the_spec = 78.334297;
R(1,1).chi_spec = 94.8872;
R(1,1).phi_spec = -1.05753;
% scan parameters
R(1,1).inc_the_spec = 0.005; %scanning step increment in theta in degrees
R(1,1).dd = 1.4; % sample to detector distance in m
R(1,1).pix = 55*10^-6; % detector pixel size in m - this is effective, so larger than actual if binnig is used...
R(1,1).lam = 0.12398 * 10^-9;
% twin ??
R(1,1).twin = 1;


%%%%%  -10-1 reflection of dislo 7 - scans 1132 to 1215: 
R(2,1).ref = [-1 0 -1]';
Amp = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_-10-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_-10-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW-AMP.mat', '-mat')));
Ph  = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_-10-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_-10-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW-PH.mat', '-mat')));
R(2,1).dat = Amp(ss1, ss2, ss3).* exp(1i * Ph(ss1, ss2, ss3));
% detector angles as given in spec file
R(2,1).gam_spec = 20.3994;
R(2,1).del_spec = 25.39525;
% sample stack angles as given in spec file
R(2,1).the_spec = -95.847564;
R(2,1).chi_spec = 85.0139;
R(2,1).phi_spec = 0.37276;
% scan parameters
R(2,1).inc_the_spec = 0.005; %scanning step increment in theta in degrees
R(2,1).dd = 1.4; % sample to detector distance in m
R(2,1).pix = 55*10^-6; % detector pixel size in m - this is effective, so larger than actual if binnig is used...
R(2,1).lam = 0.12398 * 10^-9;
% twin??
R(2,1).twin = 1;


%%%%%  1-10 reflection of dislo 7 - scans 1270 to 1327: 
R(3,1).ref = [1 -1 0]';
Amp = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_1-10_fullres_CPU_4000-10-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_1-10_fullres_CPU_4000-10-ERlrHIOlr4000-CVl-SW-AMP.mat', '-mat')));
Ph  = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_1-10_fullres_CPU_4000-10-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_1-10_fullres_CPU_4000-10-ERlrHIOlr4000-CVl-SW-PH.mat', '-mat')));
R(3,1).dat = Amp(ss1, ss2, ss3).* exp(1i * Ph(ss1, ss2, ss3));
% detector angles as given in spec file
R(3,1).gam_spec = 0.1363;
R(3,1).del_spec = 32.3985;
% sample stack angles as given in spec file
R(3,1).the_spec = 35.181896;
R(3,1).chi_spec = 90.2294;
R(3,1).phi_spec = -0.32771;
% scan parameters
R(3,1).inc_the_spec = 0.005; %scanning step increment in theta in degrees
R(3,1).dd = 1.4; % sample to detector distance in m
R(3,1).pix = 55*10^-6; % detector pixel size in m - this is effective, so larger than actual if binnig is used...
R(3,1).lam = 0.12398 * 10^-9;
% twin??
R(3,1).twin = 1;


%%%%%  0-1-1 reflection of dislo 7 - scans 1338 to 1395: 
R(4,1).ref = [0 -1 -1]';
Amp = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_0-1-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_0-1-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW-AMP.mat', '-mat')));
Ph  = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_0-1-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_0-1-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW-PH.mat', '-mat')));
R(4,1).dat = Amp(ss1, ss2, ss3).* exp(1i * Ph(ss1, ss2, ss3));
% detector angles as given in spec file
R(4,1).gam_spec = 20.3612;
R(4,1).del_spec = 25.432;
% sample stack angles as given in spec file
R(4,1).the_spec = -11.771695;
R(4,1).chi_spec = 88.7905;
R(4,1).phi_spec = -4.85186;
% scan parameters
R(4,1).inc_the_spec = 0.005; %scanning step increment in theta in degrees
R(4,1).dd = 1.4; % sample to detector distance in m
R(4,1).pix = 55*10^-6; % detector pixel size in m - this is effective, so larger than actual if binnig is used...
R(4,1).lam = 0.12398 * 10^-9;
% twin??
R(4,1).twin = 1;


%%%%%  110 reflection of dislo 7 - scans 1416 to 1473: 
R(5,1).ref = [1 1 0]';
Amp = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_110_fullres_CPU_4000-110-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_110_fullres_CPU_4000-110-ERlrHIOlr4000-CVl-SW-AMP.mat', '-mat')));
Ph  = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_110_fullres_CPU_4000-110-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_110_fullres_CPU_4000-110-ERlrHIOlr4000-CVl-SW-PH.mat', '-mat')));
R(5,1).dat = Amp(ss1, ss2, ss3).* exp(1i * Ph(ss1, ss2, ss3));
% detector angles as given in spec file
R(5,1).gam_spec = 1.8164;
R(5,1).del_spec = 32.32275;
% sample stack angles as given in spec file
R(5,1).the_spec = 124.926;
R(5,1).chi_spec = 91.6375;
R(5,1).phi_spec = 1.14848;
% scan parameters
R(5,1).inc_the_spec = 0.005; %scanning step increment in theta in degrees
R(5,1).dd = 1.4; % sample to detector distance in m
R(5,1).pix = 55*10^-6; % detector pixel size in m - this is effective, so larger than actual if binnig is used...
R(5,1).lam = 0.12398 * 10^-9;
% twin??
R(5,1).twin = 1;


%%%%%  01-1 reflection of dislo 7 - scans 1484 to 1541:
R(6,1).ref = [0 1 -1]';
Amp = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_01-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_01-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW-AMP.mat', '-mat')));
Ph  = cell2mat(struct2cell(load('Rec-D_gpu_FH_0006_dislo7_01-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW/Rec-D_gpu_FH_0006_dislo7_01-1_fullres_CPU_4000-1-ERlrHIOlr4000-CVl-SW-PH.mat', '-mat')));
R(6,1).dat = Amp(ss1, ss2, ss3).* exp(1i * Ph(ss1, ss2, ss3));
% detector angles as given in spec file
R(6,1).gam_spec = 23.3567;
R(6,1).del_spec = 22.65075;
% sample stack angles as given in spec file
R(6,1).the_spec = 175.04299;
R(6,1).chi_spec = 90.1146;
R(6,1).phi_spec = 4.99869;
% scan parameters
R(6,1).inc_the_spec = 0.005; %scanning step increment in theta in degrees
R(6,1).dd = 1.4; % sample to detector distance in m
R(6,1).pix = 55*10^-6; % detector pixel size in m - this is effective, so larger than actual if binnig is used...
R(6,1).lam = 0.12398 * 10^-9;
% twin??
R(6,1).twin = 1;


% Size of the output image in sample coordinate space. 
Sspace.N1 = 240; % Y direction 
Sspace.N2 = 240; % X direction 
Sspace.N3 = 240; % Z direction
% pixel size in sample coordinates in m
Sspace.pix = 5*10^-9;

%% other parameters
thr = 0.3; %threshold for intensity cutoff...
incl_ref = [1:6]; % numbers of reflections to include in displacement field fitting. at least 3 are needed...
a0 = 316.52*10^-12;% lattice parameter in m...
Sspace.ups = 1; % factor by which the data in the detector plane is up-sampled before mapping to the sample frame

% ramp removal parameters
ups_ph_rr_det = 1; %upsampling for initial phase ramp removal in det coordinates
ups_ph_rr_sam = 1; %upsampling for ramp removal in sample coordinates

% phase offsets applied to each reflection to move phase jumps around. 
ph_off = [0, pi/2, -pi/2]; % these are phase offsets applied to shift the phase ramps around...

% plot titles
crystal_name = 'dislo7'; % titles of figures

%save displacements once computed
save_disp = 0;
disp_save_name = 'dislo7';

%save figure once computed
fig_save = 0;
figure_save_name = 'dislo7';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   No need to change anything below here... %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Treat each reflection in turn...
for iii = 1 : size(R,1)
    R_dat = R(iii,1).dat;
    
    %% correct for twining artefacts:
    if R(iii,1).twin == 1 
        R_dat = conj_reflect(R_dat); 
    end

    %% correct for cropping of the detector frame...
    R(iii,1).inc_the_spec = R(iii,1).inc_the_spec./fact;
    R(iii,1).pix = R(iii,1).pix./fact;


    %% do initial round of phase ramp removal in the detector frame....
    tic
    R_dat_rr = remove_ramp_pn_ups(R_dat, ups_ph_rr_det); %remove phase ramp by centering FFT.
    toc


    %% Move jumps by adding constant phase offsets in detector coordiantes...
    % compute the offsets that should be applied to different reflections: 
    Ph_us = angle(R_dat_rr); % unshifted phase...
    for jjj = 1: size(ph_off,2)
        Ph_temp = Ph_us + ph_off(1, jjj); % add the phase offset...
        Ph_temp(Ph_temp>pi) = Ph_temp(Ph_temp>pi) - 2*pi; % return anything greater than pi into the -pi to pi range
        Ph_temp(Ph_temp<-pi) = Ph_temp(Ph_temp<-pi) + 2*pi; %return anything smaller than pi into the -pi to pi range
        R_ph_rr_oph(:,:,:,jjj) = Ph_temp; % Save offset phases for all reflections...
    end
    clear Ph_us Ph_temp
    
    %% project everything into the sample space:
    % project the full, unshifted complex electron density into the wafer space. 
    [Sam_dat, Xgrid, Ygrid, Zgrid, Sam(iii,1).q_sam] = Support_transfer_ref2sample_ups(R_dat_rr, R(iii,1), Sspace, 'none'); % this is the full, unshifted complex density projected into the sample space. 

    % project the shifted phases into the wafer space.
    for jjj = 1: size(ph_off,2)
        [samcom, ~, ~, ~, ~] = Support_transfer_ref2sample_ups((abs(R_dat_rr).*exp(1i*R_ph_rr_oph(:,:,:,jjj))), R(iii,1), Sspace, 'none'); % project shifted phases into sample space.
        Sam_oph(:,:,:,jjj) = angle(samcom);
        clear samcom
    end

    clear R_dat R_dat_rr R_ph_rr_oph 
    
    %% do another round of phase ramp removal in wafer coordinates:
    for jjj = 1 : size(ph_off,2)
        tic
        amp_ph = abs(Sam_dat).* exp(1i * Sam_oph(:,:,:,jjj));
        AAA = remove_ramp_pn_ups(amp_ph, ups_ph_rr_sam);   % Sam(iii,jjj).oph_rr is complex desnity in lab coordinates with different phase offsets.      
        Sam(iii,jjj).oph_rr = angle(AAA);
        if jjj == 1
            Sam(iii,1).a_rr = abs(AAA);
        end
        toc
    end
    clear AAA Sam_dat Sam_oph
end


% calculate centre of mass of different reflections and shift reflections so CoM is at zero to nearest pixel....
% use mask for each reflection....
for iii =  1: size(R,1)
    mask = Sam(iii,1).a_rr;
    mask(mask>thr) =1;
    mask(mask<1) = 0;
    se = strel('sphere', 3);
    mask = imerode(imdilate(mask, se),se);   
    CoM_rr(iii,:) = (centerOfMass(mask)-(size(mask)+1)/2); %work out offset of CoM from centre of array...
    Sam(iii,1).a_rr = circshift(Sam(iii,1).a_rr, -round(CoM_rr(iii,:)));    %shift the amplitudes...
    for jjj = 1 : size(ph_off,2)
        Sam(iii,jjj).oph_rr = circshift(Sam(iii,jjj).oph_rr, -round(CoM_rr(iii,:)));  %shift the phases...
    end    
end


%% plot crystal phases in the wafer frame...
% 3 figures are plotted, one for each phase offset
plorder = [5, 2, 4, 3, 1, 6]; % order in which to plot reflections...
plot_sam_coords = 1;
if plot_sam_coords == 1;
    for jjj = 1 : size(ph_off,2)
        vp = [0 0 1]; %view point to be used for plots. Magnitude doesn't matter
        ph_caxis = [-pi pi];
        figure; 
        for iii = 1: size(R,1);
            P(plorder(iii)) = subplot(2,3,plorder(iii));
            [faces,verts,colors] = isosurface(Xgrid, Ygrid, Zgrid,Sam(iii,1).a_rr,thr,Sam(iii,jjj).oph_rr);
            patch('Vertices', verts, 'Faces', faces, ...
                'FaceVertexCData', colors, ...
                'FaceColor','interp', ...
                'edgecolor', 'none');
            daspect([1,1,1]);
            view(vp);
            axis equal vis3d xy;
            lighting flat
            caxis(ph_caxis)
            hold on
            %plot coordiantes:
            %quiver3(-750,-750,-750,300,0,0, 'r'); quiver3(-750,-750,-750,0,300,0, 'g'); quiver3(-750,-750,-750,0,0,300, 'b');
            %plot q
            quiver3(0, 0, 0, 10^-17*Sam(iii,1).q_sam(1,1), 10^-17*Sam(iii,1).q_sam(2,1), 10^-17*Sam(iii,1).q_sam(3,1), 'k')
            axis([-5 5 -5 5 -5 5]*10^-7)
            %colorbar
            title(num2str(R(iii,1).ref'))
        end
        suptitle(['Phases in sample coordinates. Phase offset: ', num2str(ph_off(jjj))])
        
        % add a colorbar and improve positioning of everything...
        h=colorbar;
        set(h, 'Position', [.9 .14 .03 .68])
        pos1 = get(P(1), 'Position');
        for i=1:6
            pos=get(P(i), 'Position');
            set(P(i), 'Position', [pos(1)-0.7*pos1(1) 0.9*pos(2)-0.0*pos1(2) pos1(3)*1.2 pos1(4)*1.2]);
        end
    end
end


%% plot crystal amplitudes in the wafer frame...

figure;
for iii = 1: size(R,1);
    P(plorder(iii)) = subplot(2,3,plorder(iii));
    pat = patch(isosurface(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9,Sam(iii,1).a_rr,thr));
    isonormals(Xgrid, Ygrid, Zgrid,Sam(iii,1).a_rr,pat);
    set(pat,'FaceColor','yellow','EdgeColor','none');
    alpha(pat,0.08)
    hold on;
    daspect([1,1,1]);
    view([0 0 1]);
    axis equal vis3d xy off;
    camlight
    %hhh(kkk) = camlight('headlight');
    %set(hhh(kkk), 'Position', [-254.9148 7.2725e+03 4.6953e+03])
    lighting flat
    hold on
    camup([0 1 0])
    view(-180, -55);
    %plot coordiantes:
    quiver3(-490,-490,-490,500./0.9,0,0, 'r', 'LineWidth', 1.5); % oddly Matlab plots vectors too short. With the factor of 1/0.9 axis vectors are correct length, 500 nm. 
    quiver3(-490,-490,-490,0,500./0.9,0, 'g', 'LineWidth', 1.5); % oddly Matlab plots vectors too short. With the factor of 1/0.9 axis vectors are correct length, 500 nm. 
    quiver3(-490,-490,-490,0,0,500./0.9, 'b', 'LineWidth', 1.5); % oddly Matlab plots vectors too short. With the factor of 1/0.9 axis vectors are correct length, 500 nm. 
    %plot q
    qn = Sam(iii,1).q_sam ./ norm(Sam(iii,1).q_sam);
    quiver3(0, 0, 0, qn(1,1)*600, qn(2,1)*600, qn(3,1)*600, 'k', 'LineWidth', 1.5)
    axis([-500 500 -500 500 -500 500])
    t(plorder(iii)) = title([num2str(R(iii,1).ref')], 'position', [0, 700, 0], 'FontSize', 12);
end

    
%% visualising overlap between reconstructions:
fcolours = {'red'; 'green'; 'blue'; 'yellow'; 'magenta'; 'cyan'};
% plot average morphology and overlap of different reflections

% compute intensity averaged over all reflections:
amp_waf_all = zeros(size(Sam(1,1).a_rr)); % intensity averaged over all reflections:
amp_waf_incl = zeros(size(Sam(1,1).a_rr)); % intensity average over reflections included in displacment fitting...

for iii = 1: size(R,1)
    amp_waf_all = amp_waf_all + abs(Sam(iii,1).a_rr);
    if sum(incl_ref == iii)
        amp_waf_incl = amp_waf_incl + abs(Sam(iii,1).a_rr);
    end
end
amp_waf_all = amp_waf_all./size(R,1);
amp_waf_incl = amp_waf_incl./size(incl_ref,2);

% binary mask based on amplitude of all reflections...
mask_all = amp_waf_all; mask_all(mask_all<thr) = 0; mask_all(mask_all>0) = 1;


% Plot average shape recovered from included reflections and overlaps of the morphologies recovered from all reflections...
figure;
subplot(1,2,1)
p_incl = patch(isosurface(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9, amp_waf_incl,thr));
isonormals(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9,amp_waf_incl,p_incl)
set(p_incl,'FaceColor','red','EdgeColor','none');
daspect([1,1,1]);
view([0 1 0]);
axis vis3d;
camlight
lighting flat
axis([-500 500 -500 500 -500 500])
%plot coordiantes:
hold on; quiver3(-490,-490,-490,500./0.9,0,0, 'r', 'LineWidth', 1.5); % oddly Matlab plots vectors too short. With the factor of 1/0.9 axis vectors are correct length, 500 nm.
quiver3(-490,-490,-490,0,500./0.9,0, 'g', 'LineWidth', 1.5); % oddly Matlab plots vectors too short. With the factor of 1/0.9 axis vectors are correct length, 500 nm.
quiver3(-490,-490,-490,0, 0, 500./0.9, 'b', 'LineWidth', 1.5); % oddly Matlab plots vectors too short. With the factor of 1/0.9 axis vectors are correct length, 500 nm.
title('average shape')


%plotting semi-tranparent amplitudes of different included reflections....
fcolours = {'red'; 'green'; 'blue'; 'yellow'; 'magenta'; 'cyan'};
nn = 1;

subplot(1,2,2)
for iii = incl_ref
    Plo(nn) = patch(isosurface(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9, (Sam(iii,1).a_rr),thr));
    isonormals(Xgrid.*10^9, Ygrid.*10^9, Zgrid.*10^9, (Sam(iii,1).a_rr),Plo(nn))
    set(Plo(nn),'FaceColor',fcolours{iii},'EdgeColor','none');
    alpha(Plo(nn),0.05)
    hold on;
    nn = nn+1;
end
daspect([1,1,1]);
view([0 1 0]);
axis vis3d;
camlight
lighting flat
axis([-500 500 -500 500 -500 500])
%plot coordiantes:
hold on; quiver3(-490,-490,-490,500./0.9,0,0, 'r', 'LineWidth', 1.5); % oddly Matlab plots vectors too short. With the factor of 1/0.9 axis vectors are correct length, 500 nm.
quiver3(-490,-490,-490,0,500./0.9,0, 'g', 'LineWidth', 1.5); % oddly Matlab plots vectors too short. With the factor of 1/0.9 axis vectors are correct length, 500 nm.
quiver3(-490,-490,-490,0, 0, 500./0.9, 'b', 'LineWidth', 1.5); % oddly Matlab plots vectors too short. With the factor of 1/0.9 axis vectors are correct length, 500 nm.
title('shapes of included reflections')

suptitle([crystal_name, ' crystal morphology in wafer coordinates'])

if fig_save==1
    saveas(fig6,[figure_save_name, ' crystal morphology in wafer coordinates'],'fig')
    print(fig6,[figure_save_name, ' crystal morphology in wafer coordinates'], '-dpng', '-r600')
end


%% work out strain fields in lab coordinates: 
% Take derivatives of phases in lab coordinates of the different phase shifted versions and keep whichever has the smaller value. 

for iii =  1 : size(R,1)
    % pre-assign space for derivatives
    gr_shx = zeros(size(Sam(iii,1).oph_rr,1), size(Sam(iii,1).oph_rr,2), size(Sam(iii,1).oph_rr,3), size(ph_off,2));
    gr_shy = gr_shx; gr_shz = gr_shx;
    
    %calculate derivatives
    for jjj = 1: size(ph_off,2)
        [gr_shx(:,:,:,jjj), gr_shy(:,:,:,jjj), gr_shz(:,:,:,jjj)] = ...
            gradient(Sam(iii,jjj).oph_rr);
    end
    % now need to filter out smallest gradient:
    % x direction first 
    [~, I] = min(abs(gr_shx),[],4);
    [dim1, dim2, dim3] = ndgrid(1:size(gr_shx,1), 1:size(gr_shx,2), 1:size(gr_shx,3));
    grx_tmp = gr_shx(sub2ind(size(gr_shx),dim1,dim2, dim3, I));
    grx(:,:,:,iii) = reshape(grx_tmp, size(Sam(iii,1).oph_rr,1),size(Sam(iii,1).oph_rr,2), size(Sam(iii,1).oph_rr,3))./Sspace.pix;
    clear I
    clear grx_tmp gr_shx
    
    % y direction
    [~, I] = min(abs(gr_shy),[],4);
    gry_tmp = gr_shy(sub2ind(size(gr_shy),dim1,dim2, dim3, I));
    gry(:,:,:,iii) = reshape(gry_tmp, size(Sam(iii,1).oph_rr,1),size(Sam(iii,1).oph_rr,2), size(Sam(iii,1).oph_rr,3))./Sspace.pix;
    clear I
    clear gry_tmp gr_shy
    
    % z direction 
    [~, I] = min(abs(gr_shz),[],4);
    grz_tmp = gr_shz(sub2ind(size(gr_shz),dim1,dim2, dim3, I));
    grz(:,:,:,iii) = reshape(grz_tmp, size(Sam(iii,1).oph_rr,1),size(Sam(iii,1).oph_rr,2), size(Sam(iii,1).oph_rr,3))./Sspace.pix;
    clear I
    clear grz_tmp gr_shz

end

%% Find displacement gradients in wafer coordiantes:
% d_phi/dx = q. [dux/dx; duy/dx; duz/dx]

% establish q matrix...
for kkk = 1: size(incl_ref,2)
    iii = incl_ref(kkk);
    % matrix of q vectors in wafer coordinates!
    q_waf_mat(kkk,:) = Sam(iii,1).q_sam';
end

% x-derivative terms only: 
for kkk = 1: size(incl_ref,2)
    iii = incl_ref(kkk);
    dphi_dx(kkk,:) = squeeze(reshape(grx(:,:,:,iii),[],1,1)); 
end
% work out x-derivatives - this is the correct way around ...
du_dx_all = q_waf_mat\dphi_dx;
% sort of the terms and get into correct format...
dux_dx = reshape(du_dx_all(1,:), size(Sam(iii,jjj).oph_rr,1), size(Sam(iii,jjj).oph_rr,2), size(Sam(iii,jjj).oph_rr,3));
duy_dx = reshape(du_dx_all(2,:), size(Sam(iii,jjj).oph_rr,1), size(Sam(iii,jjj).oph_rr,2), size(Sam(iii,jjj).oph_rr,3));
duz_dx = reshape(du_dx_all(3,:), size(Sam(iii,jjj).oph_rr,1), size(Sam(iii,jjj).oph_rr,2), size(Sam(iii,jjj).oph_rr,3));
clear dphi_dx du_dx_all

% y-derivative terms only: 
for kkk = 1: size(incl_ref,2)
    iii = incl_ref(kkk);
    dphi_dy(kkk,:) = squeeze(reshape(gry(:,:,:,iii),[],1,1)); 
end
% work out y-derivatives - this is the correct way around ...
du_dy_all = q_waf_mat\dphi_dy;
% sort of the terms and get into correct format...
dux_dy = reshape(du_dy_all(1,:), size(Sam(iii,jjj).oph_rr,1), size(Sam(iii,jjj).oph_rr,2), size(Sam(iii,jjj).oph_rr,3));
duy_dy = reshape(du_dy_all(2,:), size(Sam(iii,jjj).oph_rr,1), size(Sam(iii,jjj).oph_rr,2), size(Sam(iii,jjj).oph_rr,3));
duz_dy = reshape(du_dy_all(3,:), size(Sam(iii,jjj).oph_rr,1), size(Sam(iii,jjj).oph_rr,2), size(Sam(iii,jjj).oph_rr,3));
clear dphi_dy du_dy_all


% z-derivative terms only: 
for kkk = 1: size(incl_ref,2)
    iii = incl_ref(kkk);
    dphi_dz(kkk,:) = squeeze(reshape(grz(:,:,:,iii),[],1,1)); 
end
% work out z-derivatives - this is the correct way around ...
du_dz_all = q_waf_mat\dphi_dz;
% sort of the terms and get into correct format...
dux_dz = reshape(du_dz_all(1,:), size(Sam(iii,jjj).oph_rr,1), size(Sam(iii,jjj).oph_rr,2), size(Sam(iii,jjj).oph_rr,3));
duy_dz = reshape(du_dz_all(2,:), size(Sam(iii,jjj).oph_rr,1), size(Sam(iii,jjj).oph_rr,2), size(Sam(iii,jjj).oph_rr,3));
duz_dz = reshape(du_dz_all(3,:), size(Sam(iii,jjj).oph_rr,1), size(Sam(iii,jjj).oph_rr,2), size(Sam(iii,jjj).oph_rr,3));
clear dphi_dz du_dz_all


%% Assemble strain tensor and lattice rotation tensor...

% format: eps_waf: dims 1,2,3, spatial dimensions. 4th and 5th dimensions are the strain tensor. 
% i.e. :,:,:,1,1 refers to xx, :,:,:,1,2 refers to xy component etc. 
eps_waf(:,:,:,1,1) = dux_dx;                %xx
eps_waf(:,:,:,1,2) = 0.5.*(dux_dy+duy_dx);  %xy
eps_waf(:,:,:,1,3) = 0.5.*(dux_dz+duz_dx);  %xz
eps_waf(:,:,:,2,1) = 0.5.*(duy_dx+dux_dy);  %yx
eps_waf(:,:,:,2,2) = duy_dy;                %yy
eps_waf(:,:,:,2,3) = 0.5.*(duy_dz+duz_dy);  %yz
eps_waf(:,:,:,3,1) = 0.5.*(duz_dx+dux_dz);  %zx
eps_waf(:,:,:,3,2) = 0.5.*(duz_dy+duy_dz);  %zy
eps_waf(:,:,:,3,3) = duz_dz;                %yy

% format: rot_waf: dims 1,2,3, spatial dimensions. 4th dimension is the lattice rotation vector. 
% i.e. :,:,:,1 refers to (rot x), :,:,:,2 refers to (rot y) component etc. 
rot_waf(:,:,:,1) = 0.5.* (duz_dy - duy_dz);
rot_waf(:,:,:,2) = 0.5.* (dux_dz - duz_dx);
rot_waf(:,:,:,3) = 0.5.* (duy_dx - dux_dy);

%% Save outputs for dislocation analysis and plotting
% the file 'recD_sam_str_rot.mat' is needed for the further dislocation analysis and plotting. 
for iii = 1:6; Sam_red(iii,1).q_sam = Sam(iii,1).q_sam; end
Sam_red(1,1).eps_waf = eps_waf;
Sam_red(1,1).rot_waf = rot_waf;
for iii = 1:6; Sam_red(iii,1).a_rr = Sam(iii,1).a_rr; end
save('recD_sam_str_rot.mat', 'Sam_red', 'Xgrid', 'Ygrid', 'Zgrid', '-v7.3')


%% Make a video of reconstructed lattice strain and rotation

v = VideoWriter('strain_rot_tensor.mp4', 'MPEG-4');
v.Quality = 95;
open(v)
calim = [-0.0008 0.0008];

figure;
%colormap(RdYeBlue);
for aaa = 60:184 % size(Sam(iii,1).oph_rr, 2) This gives correct range for 5 nm real-space pixel size... 
P(1) = subplot(3,3,1); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_waf(:,aaa,:,1,1))); view([0 0 1])
shading flat; axis equal; axis([-500 500 -500 500]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{xx}')

P(2) = subplot(3,3,2); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_waf(:,aaa,:,1,2))); view([0 0 1])
shading flat; axis equal; axis([-500 500 -500 500]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{xy}')

P(3) = subplot(3,3,3); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_waf(:,aaa,:,1,3))); view([0 0 1])
shading flat; axis equal; axis([-500 500 -500 500]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{xz}')

P(5) = subplot(3,3,5); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_waf(:,aaa,:,2,2))); view([0 0 1])
shading flat; axis equal; axis([-500 500 -500 500]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{yy}')

P(6) = subplot(3,3,6); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_waf(:,aaa,:,2,3))); view([0 0 1])
shading flat; axis equal; axis([-500 500 -500 500]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{yz}')

P(9) = subplot(3,3,9); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(eps_waf(:,aaa,:,3,3))); view([0 0 1])
shading flat; axis equal; axis([-500 500 -500 500]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\epsilon_{zz}')

P(4) = subplot(3,3,4); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(rot_waf(:,aaa,:,3))); view([0 0 1])
shading flat; axis equal; axis([-500 500 -500 500]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\omega_{z}')

P(7) = subplot(3,3,7); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(rot_waf(:,aaa,:,2))); view([0 0 1])
shading flat; axis equal; axis([-500 500 -500 500]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\omega_{y}')

P(8) = subplot(3,3,8); 
str_plo = surf(squeeze(Zgrid(:,aaa,:)*10^9), squeeze(Ygrid(:,aaa,:)*10^9), squeeze(rot_waf(:,aaa,:,1))); view([0 0 1])
shading flat; axis equal; axis([-500 500 -500 500]); caxis(calim); alpha(str_plo, squeeze(mask_all(:,aaa,:)));
title('\omega_{x}')


h=colorbar;
set(h, 'Position', [.9 .18 .03 .68])
pos1 = get(P(1), 'Position');
for i=1:9
    pos=get(P(i), 'Position');
    set(P(i), 'Position', [pos(1)-0.5*pos1(1) 1.05*pos(2)-0.08*pos1(2) 1.1*pos(3) 1.1*pos(4)]);
end

suptitle(['Lattice strain and rotation. Xpos: ', num2str(Xgrid(1,aaa,1)*10^9), ' nm'])

% make video 
frame = print('-RGBImage', '-r300');
writeVideo(v,frame) ;
end
close(v)


