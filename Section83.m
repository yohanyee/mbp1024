%% Radon and Sinogram 
clear

%% Read in the dicom files

original=load('C:\Users\Brandon\Documents\MATLAB\CTLab\Raw\Section75.mat');

data=original.data;
data_low=original.data_lowdose;
theta=original.theta;
p=original.p;

%Determine size and resolution of sinograms

N=length(data);
l=1:N;
l=p*(l-(N/2));
r=-N/2:N/2;

%Plot sinogram

imshow(data, [], 'XData', theta, 'YData', l)
axis square
axis on
figure
 
%% Create various reconstructions of the sinogram data


res=1:512;
res=p*(res-(512/2));

%laminogram

filter1='none';

recon1=iradon(data, theta, 'linear', filter1, 1, 512)*16;
% recon1=recon1*20;
recon1_low=iradon(data_low, theta, 'linear', filter1, 1, 512);

imshow(recon1, [], 'xdata', res, 'ydata', res)
% axis on
figure


%%Filtered Backprojection

filter2='Ram-Lak';


%Using projections from 0 and 90 degrees
angle1=[0,90];
data1=[data(:,1), data(:,181)];
recon2=iradon(data1, angle1, 'linear', filter2, 1, 512)*16;


imshow(recon2, [],  'xdata', res, 'ydata', res)
% axis on
figure


%Using every 10th projection
m=1:10:360;

for i=1:length(m)
    angle2(i)=theta(m(i));
    data2(:,i)=data(:,m(i));
end

recon3=iradon(data2, angle2, 'linear', filter2, 1, 512)*16;


imshow(recon3, [], 'xdata', res, 'ydata', res)
% axis on
figure


%Using all prjections

recon4=iradon(data, theta, 'linear', filter2, 1, 512)*16;


imshow(recon4, [], 'xdata', res, 'ydata', res)
% axis on



