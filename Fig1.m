%%% code matlab of Fig1 BD-RPCA%%%%%                                              ;
clear  all;
close all 
%% Add Path
running_folder = 'C:\Users\dpham\ownCloud\Working\Atempo\';
addpath(genpath(fullfile(running_folder,'BD-RPCA-GitHub')));
addpath(genpath(fullfile(running_folder,'Data')));

%%
load psf_simu;
load pht_data.mat;
load mask.mat;

f0=3e6;                  %  Transducer center frequency [Hz]
fs=9e6;                %  Sampling frequency [Hz]
c=1540;       
no_lines=300;              %  Number of lines in image
image_width=100/1000;      %  Size of image sector
d_x=image_width/no_lines; %  Increment for image
d_z = (c/fs)/2;

x = min(phantom_positions(:,1)):d_x:max(phantom_positions(:,1));
z = min(phantom_positions(:,3)):d_z:max(phantom_positions(:,3));

[X Z] = meshgrid(x,z);
refl = griddata(phantom_positions(:,1),phantom_positions(:,3),phantom_amplitudes,X,Z);
refl(find(isnan(refl)==1)) = 0;
refl = refl(:,:)/max(refl(:));
refl=refl(700:1150,40:200);
Nt= 400; %number of layers

region1=refl(200:270,45:57);
region2=refl(145:180,90:100);

for i=1:Nt
   refl_prov_3D(:,:,i)=refl;
   region1= circshift(region1,-1,1);
   region1= circshift(region1,1,2);
   region2= circshift(region2,1,1);
   region2= circshift(region2,1,2); 
   refl_prov_3D(200:270,45:57,i)=region1;
   refl_prov_3D(145:180,90:100,i)=region2;
end

% figure,
% for i=1:Nt
% imagesc((refl_prov_3D(:,:,i)));colormap('gray')
% drawnow
% end

for i=1:Nt
   rf_image_3D(:,:,i) = conv2(refl_prov_3D(:,:,i),psf,'same');
   rf_image=rf_image_3D(:,:,i);
%    rf_image=refl_prov_3D(:,:,i);
   rf_image_3D(:,:,i) = rf_image/max(rf_image(:));
end


% figure,
% for i=1:Nt
% imagesc(rf2bmode(rf_image_3D(:,:,i),1));colormap('gray')
% drawnow
% end

Fig = figure;
rf_image_3D(200:270,45,1)=0;
rf_image_3D(200:270,57,1)=0;
rf_image_3D(200,45:57,1)=0;
rf_image_3D(270,45:57,1)=0;
rf_image_3D(145:180,90,1)=0;
rf_image_3D(145:180,100,1)=0;
rf_image_3D(145,90:100,1)=0;
rf_image_3D(180,90:100,1)=0;
imagesc(rf2bmode(rf_image_3D(:,:,1),1)); colormap('gray')
[Nz,Nx,Nt] = size(rf_image_3D);    
espace_xx=1:Nx;
espace_zz=1:Nz;
xlabel('N_X [nb]')  ; 
ylabel('N_Z [nb]')                                  ; 
set(gca,'XTick', 1:40:espace_xx(end));
set(gca,'YTick', 1:100:espace_zz(end)); 
hold on
rectangle('Position',[45 200 12 70],'EdgeColor','r','LineWidth',2)
rectangle('Position',[90 145 10 35],'EdgeColor','r','LineWidth',2)

result_folder = fullfile(running_folder,'BD-RPCA-GitHub','Results');

%% Change if 1 if need to print this figure. 

if 0
export_fig(Fig, ... % figure handle
        sprintf('%s/simuConv', result_folder),... % name of output file without extension
        '-painters', ...      % renderer
        '-transparent', ...   % renderer
        '-pdf', ...         % file format
        '-r500',... 
        '-nocrop' );             % resolution in dpi  
end
% figure;subplot(121), imagesc(rf2bmode(refl,0.1)); colormap gray
% title('TRF')
% subplot(122),imagesc(rf2bmode(rf_image(:,:,1),0.1)); colormap gray
% title('Bmode image')



