close all

Fs = 1000;
allow_save = true;

N_comps = 2;%2 %ncomps for pca
Nhigh = 3;%3 % upper freq limit (including) (DC considered as freq with index 0)
Nlow = 1; %1%lower limit (including)
DoWhitening = false;
pre_time = 20;%20
post_time = 30;%30
T = pre_time+post_time+1; %num of samples
Nsim = 200; % Nsim for wave and Nsim for static
SNR = 3;


params = struct();
params.duration = T/Fs;
params.sampling_rate = Fs;
params.nspoints = 51;%for 0.5 mm/ms - 20 mm
params.speed = 0.8;
params.draw_paths = false;
params.draw_wave = false;



subj_name = '0026';
protocol_dir = '/home/ksasha/anaconda3/Documents/brainstorm_db/F0026/';
file_name = 'block001';
channel_type = 'grad'; % channels you want to analyse ('grad' or 'mag')
  
cortex = load('C:/Users/Fedosov/Documents/projects/waves/tess_cortex_pial_high.mat');


% Initial parameters
if strcmp(channel_type, 'grad') == 1
    channel_idx = setdiff(1:306, 3:3:306);
elseif strcmp(channel_type, 'mag') == 1
    channel_idx = 3:3:306;
end

G3_dir=  'C:/Users/Fedosov/Downloads/K_waves/data/headmodel_surf_os_meg.mat';
G = gain_orient(G3_dir,channel_idx);








rng(2)
F_storage = zeros(1,Nsim+Nsim);

waves = zeros(Nsim+Nsim,length(channel_idx),T);
blobs = zeros(Nsim,length(channel_idx),T);
for i = 1:Nsim
    i
 
    [wave_raw,strt,nn]= wave_on_sensor(cortex,params, G);
    blob_raw = blob_on_sensor(cortex,params, G,strt,nn);

    wave_norm = wave_raw/norm(wave_raw);
    blob_norm = blob_raw/norm(blob_raw);
    
    noise = generate_brain_noise(G, 1000, T, Fs);
    noise_norm = noise/norm(noise);
    waves(i,:,:) = SNR*wave_norm +noise_norm;
    waves(i+Nsim,:,:) = SNR*blob_norm +noise_norm;
end



labels = zeros(1,Nsim+Nsim);
labels(1:Nsim) = 1;


save('C:/Users/Fedosov/Documents/projects/waves/Wave_Catcher/waves_speed_8.mat','waves','labels','N_comps',...
'Nhigh','Nlow','DoWhitening','pre_time','post_time','T','Nsim','SNR','params','subj_name','channel_type');
  


for k = 1:Nsim*2
k
y  = squeeze(waves(k,:,:))';
[u s v] = svd(y,'econ');
        
pcs = y*v(:,1:N_comps);%*v(:,1:N_comps)';
if DoWhitening
        pcs = pcs./diag(s(1:N_comps,1:N_comps))';
end
Nch = N_comps;
X = pcs;    
        
% make H matrix for symmetric model containing  6 free variables (cols)
Hsym = zeros(Nch^2,(Nch^2-Nch)/2+Nch);
Pos = reshape([1:Nch^2],Nch,Nch);
counter = 1;
for i = 1:Nch
     for j = i:Nch
       Hsym([Pos(i,j),Pos(j,i)],counter) = 1;
       counter = counter+1;
     end 
 end

 Mask = zeros(size(X,1),1);
 Mask(1) = 0; % note we also get rid of the DC
 Mask(Nlow+1:Nhigh+1) = 1;
 Mask(end-Nhigh+1:end-Nlow+1) = 1; % make frequency domain window symmetric 
    
 XF = fft(X).*(Mask*ones(1,Nch));
 j = sqrt(-1);
 N = size(X,1);
 kk = (0:N-1)';

        XFW = XF.*(exp(-j*2*pi*kk/N)*ones(1,Nch));
        dXF = XF - XFW;
        XFWblk = [];
        for i = 1:Nch
        XFWblk =  blkdiag(XFWblk,XFW); % 300 x 9
        end

        mopt_ful_fft = pinv(XFWblk'*XFWblk)*XFWblk'*dXF(:);
        mopt_sym_fft = pinv((XFWblk*Hsym)'*(XFWblk*Hsym))*(XFWblk*Hsym)'*dXF(:);
        
        
        %mopt_tra_fft = travelling_matrix(XFW,dXF);  %%%%
        

        Chi2_ful_fft = norm(XFWblk*mopt_ful_fft - dXF(:));
        Chi2_sym_fft = norm(XFWblk*Hsym*mopt_sym_fft - dXF(:));
        
        %Chi2_sym_fft = norm(XFW*mopt_tra_fft-dXF);%%%%

        F_fft = Chi2_sym_fft/Chi2_ful_fft;
        F_fft = (Chi2_sym_fft^2-Chi2_ful_fft^2)/Chi2_ful_fft^2;
        
        %eigvals = eig(mopt_tra_fft);
        %F_fft = norm(real(eigvals))/norm(imag(eigvals));
        
        
        %F_fft = Chi2_skew_fft/Chi2_sym_fft;
            
       
      
        F_storage(k) = F_fft;
        
        
        
        
      
  
        
        


%     edges = [1:0.2:6];
%     figure
%     histogram(Ff_storage{n_clust},edges,'FaceColor',c(n_clust,:))
%     ylim([0,50])
% 
%     title(['cluster ',num2str(n_clust)]);
%     if allow_save
%         saveas(gcf,['/home/ksasha/projects/jPCA/spikes_0026_new/_', num2str(n_clust),'_distrib.png']);
%     end    
%     hold on
%     histogram(Ffake_storage{n_clust},edges,'FaceColor','k')



end

close all

[X,Y,Thr,AUC] = perfcurve(labels,F_storage,1);

figure
plot(X,Y)
AUC

save('C:/Users/Fedosov/Documents/projects/waves/Wave_Catcher/speed_8_xyAUC.mat','X','Y','Thr','AUC');
  



