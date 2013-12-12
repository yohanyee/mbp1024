%% Section 8.2 %%
%%
%Load XrayObject.mat

clear
object=load('XrayObject.mat');

%create proper res

dist=256*object.p;
res=1:256;
res=object.p*(res-(256/2));


A=object.Object;
imshow(A, [], 'xdata', res, 'ydata', res)
axis square
% axis on



%create Radon transform of object at theta=0 degrees

[radon_0, xp]=radon(A,0);
radon_0=radon_0/20;

%create proper position along detector

p_new=dist/length(radon_0);
l=xp*p_new;

%Plot radon transform at 0 degrees as a function of dectector position. 0mm
%representes central axis of detector.
figure
plot(l,radon_0)

% Create a sinogam, i.e. a radon tranform at all angles, in step size of 1
% degree.

theta=0:179;

sinogram=radon(A, theta);
sinogram=sinogram/20;

figure
imshow(sinogram,[], 'xdata',theta, 'ydata', l)
axis square
axis on
