%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bicubic interpolation and extrapolation iteration method for high %
%resolution digital holographic reconstruction

%Zhengzhong Huang, Liangcai Cao*
%State Key Laboratory of Precision Measurement Technology and Instruments, %Department of Precision Instruments, Tsinghua University, Beijing 100084, %China
%hzz19@mails.tsinghua.edu.cn
%clc@tsinghua.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Zhengzhong Huang 2019
% The software is Matlab

a_simulation_hologram_plane_wave.m simulates the hologram
b_reconstruction_iterative_plane_wave.m iteratively reconstructs the hologram

subroutines:
FT2Dc.m
IFT2Dc.m
Propagator.m

should be placed into C/...MATLAB folder, as for example 
into "C:/Program Files/MATLAB/R2010b/myfiles"


clc
clear all;
close all;
%% read hologram
hologram = fopen('hologram.bin', 'r');
hologram = fread(hologram, [400,400], 'double');% Object hologram

background = fopen('background.bin', 'r');
back = fread(background, [400,400], 'double');% Background hologram without object
holo=hologram./back;% Normalized hologram

[M,N]=size(holo);
%% Experiment Parameters
wavelength=532*10^(-9);% wavelength
pitch = 3.8*10^(-6);% pitch size of detector
z=0.237;% recording distance
prop = Propagator_function(M, N, wavelength, pitch, z);% Transfer Function

%% Directly Reconstruction
measured=sqrt(holo);
recons = IFT((FT(holo)).*(prop));

%% BIP
num=2;% Interpolation rate
holo0=BIPfunction(holo,num);
[M0,N0]=size(holo0);

%% Define loop parameters
NN=2000;
A = ones(NN,NN);
phase = zeros(NN,NN);
mask = fopen('mask.bin', 'r');
support = fread(mask, [NN,NN], 'double');% Support Function
N1 = (NN - M0)/2;
N2 = (NN + N0)/2;
prop1 = Propagator_function(NN, NN, wavelength, (pitch/num), z);% Transfer Function
Loops=50;% number of iterations

%% EPI 
for tt = 1:Loops
fprintf(': %d\n', tt)

A(N1+1:N2,N1+1:N2)=holo0(1:M0,1:N0);
holo_field = A.*exp(1i.*phase);

recons1 = IFT((FT(holo_field)).*prop1);
object=1-abs(recons1);

for ii=1:N0
    for jj=1:N0
        if (object(ii,jj)<0)
            object(ii,jj)=0;
        end
    end
end

object=object.*support;
recons1 = 1-object;


holo_field_updated = IFT((FT(recons1)).*conj(prop1));
A = abs(holo_field_updated);
phase = angle(holo_field_updated);
end

figure,imshow(rot90(rot90(abs(recons1))),[]);



