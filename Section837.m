%% Section 8.3.7
%Using Shepp-Logan filter and others to optimize SNR of low dose image data

%% Call Image Data

data_struct=load('C:\Users\Brandon\Documents\MATLAB\CTLab\Raw\Section75.mat');

%Create individual variables
data=data_struct.data;
data_low=data_struct.data_lowdose;
theta=data_struct.theta;
p=data_struct.p;

%Show the original high and low dose images using a Ramp filter for
%comparison

recon=iradon(data, theta, 'linear', 'Ram-Lak', 1, 512);
recon_low=iradon(data_low, theta, 'linear', 'Ram-Lak', 1, 512);

% imshow(recon, [])
% 
% figure
% imshow(recon_low, [])

%Define which filter we use

filter='Ram-Lak';



%Vary frequency scaling for recontructions

% recon_low_1=iradon(data_low, theta, 'linear', filter, 0.001, 512);
recon_low_2=iradon(data_low, theta, 'linear', filter, 0.2, 512);
recon_low_3=iradon(data_low, theta, 'linear', filter, 0.4, 512);
recon_low_4=iradon(data_low, theta, 'linear', filter, 0.6, 512);
recon_low_5=iradon(data_low, theta, 'linear', filter, 0.8, 512);
recon_low_6=iradon(data_low, theta, 'linear', filter, 1, 512);

% figure
% imshow(recon_low_1, []);
% axis on
% 
% figure
% imshow(recon_low_2, []);
% axis on
% 
% figure
% imshow(recon_low_3, []);
% axis on
% 
% figure
% imshow(recon_low_4, []);
% axis on
% 
% figure
% imshow(recon_low_5, []);
% axis on
% 
% figure
% imshow(recon_low_6, []);
% axis on
%Outputs an 3D analyze image stack. Images 1 and 2 are the high and low
%dose data reconned with a regular Ram-Lak filter. Images 3 onwards are
%images reconned with filters other than Ram-Lak (e.g. Cosine Shepp-Logan,
%etc.) Images 3 and onwards represent declining levels of high pass
%filtering.

image3d=cat(3, recon,recon_low,recon_low_1,recon_low_2,recon_low_3,recon_low_4,recon_low_5,recon_low_6);

datatype=16;

nii=make_ana(image3d,[],[],datatype,[]);

filename=input('Name of file?', 's');

save_untouch_nii(nii, filename);


