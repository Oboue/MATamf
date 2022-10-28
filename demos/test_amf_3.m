% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets
% 
%  Copyright (C) Oboue et al., 2022
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/

%  References

%  Oboue et al., 2022
%  Huang, G., M. Bai, Q. Zhao, W. Chen, and Y. Chen, 2021, Erratic noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, 69, 101-121.
%  Chen, Y., S. Zu, Y. Wang, and X. Chen, 2020, Deblending of simultaneous-source data using a structure-oriented space varying median filter, Geophysical Journal International, 222, 1805�1�723.
%  Zhao, Q., Q. Du, X. Gong, and Y. Chen, 2018, Signal-preserving erratic noise attenuation via iterative robust sparsity-promoting filter, IEEE Transactions on Geoscience and Remote Sensing, 56, 1558-0644.
%%
clc;clear;close all;
%% addpath
addpath /Users/oboue/Desktop/MATamf/amf_data/
addpath /Users/oboue/Desktop/MATamf/amfsrcs/
addpath /Users/oboue/Desktop/MATamf/seistr/
%% load data clean data

load micro_sf_3001_334_3.mat
d=amf_scale(data(:,:,1),2);
[n1,n2]=size(d);
%% load horizontal noise 
load d_noise_micro.mat

randn('state',211111);
nnn=amf_seisdither(nn,round(amf_meanf(20*randn(1,nx),20,1,2)));

dn1=d+0.1*nnn;
%%  Denosing using the BP method
% Parameter tuning for the BP method

dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
tic
d1=amf_bandpass(dn1,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%
%% Denosing using the BP+SOSVMF method 
% Parameter tuning：add the key parameters of the SOSVMF method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=21;                   % rect:  smoothing radius (ndim*1)
rect(2)=7;                    % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=9;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
%
tic
d2=amf_bandpasssosvmf(dn1,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0.2;        % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
% 
tic
d3=amf_bandpasssosvmffk(dn1,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% Denoising data using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the key parameters of the curvelet method
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=10;               % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
tic
d4=amf_bandpasssosvmffkcurvelet(dn1,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
%% Denoising using the AMF method 
%  Parameter tuning: add the key parameters for the local orthogonalization operation

rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
rec(1) = 30;
rec(2) = 30;
rec(3) = 1;
eps1=0;               % regularization parameter, default 0.0
niter2=20;            % number of CG iterations
verb=1;               % verbosity flag (default: 0) 

%
tic
d_amf=amf(dn1,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1,rec,eps1,niter2,verb);
d5=d_amf;
toc
%
%% Vertical noise
%% load vertical noise 
load d_noise_micro.mat

randn('state',2111112);
nn2=amf_seisdither(nn2,round(amf_meanfs(10*randn(1,nt),20,1,2,100)));

dn2=d+0.2*nn2';
%%  Denosing using the BP method
% Parameter tuning for the BP method

dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
tic
d21=amf_bandpass(dn2,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%
%% Denosing using the BP+SOSVMF method 
% Parameter tuning：add the key parameters of the SOSVMF method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                    % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=8;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
%                                                                                                             
tic
d22=amf_bandpasssosvmf(dn2,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0; % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
tic
d23=amf_bandpasssosvmffk(dn2,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% Denoising data using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the key parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=2.5;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
%
tic
d24=amf_bandpasssosvmffkcurvelet(dn2,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
%% Denoising using the AMF method 
% % Parameter tuning: add the key parameters for the local orthogonalization operation

rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
rec(1) = 30;
rec(2) = 30;
rec(3) = 1;
eps1=0;               % regularization parameter, default 0.0
niter2=20;            % number of CG iterations
verb=1;               % verbosity flag (default: 0) 

%
tic
dout=amf(dn2,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1,rec,eps1,niter2,verb);
d25=dout;
toc
%
%% plot figures

dt=0.0005;
t=[0:nt-1]*dt;
h=[0:nx-1]*dx+x0;

figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
subplot(2,4,1);
imagesc(h,t,d);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('Clean','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);
%
subplot(2,4,2);
imagesc(h,t,dn1);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('Raw','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,3);
imagesc(h,t,d1);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('BP','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,4);
imagesc(h,t,d2);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('BP+SOSVMF','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,5);
imagesc(h,t,d3);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('BP+SOSVMF+FK','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,6);
imagesc(h,t,d4);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('BP+SOSVMF+FK+Curvelet','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,7);
imagesc(h,t,d5);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('AMF','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);
print(gcf,'-depsc','-r300','fig5.pdf');
print(gcf,'-depsc','-r300','fig5.eps');
%
figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
subplot(2,4,1);
imagesc(h,t,d);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('Clean','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,2);
imagesc(h,t,dn2);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('Raw','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,3);
imagesc(h,t,d21);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('BP','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,4);
imagesc(h,t,d22);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('BP+SOSVMF','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,5);
imagesc(h,t,d23);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('BP+SOSVMF+FK','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,6);
imagesc(h,t,d24);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('BP+SOSVMF+FK+Curvelet','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);

subplot(2,4,7);
imagesc(h,t,d25);caxis([-0.5,0.5]);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Lateral distance (km)','Fontsize',13,'fontweight','bold');
title('AMF','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
colormap(amf_seis);
print(gcf,'-depsc','-r300','fig6.pdf');
print(gcf,'-depsc','-r300','fig6.eps');

%% Numerical test using signal-to-noise ratio (S/N)
% S/N --- horizontal noise

amf_snr(d,dn1) % Noisy
amf_snr(d,d1) % BP
amf_snr(d,d2) % BP+SOSVMF
amf_snr(d,d3) %Bp+SOSVMF+Curvelet
amf_snr(d,d4) % BP+SOSVMF+Curvelet+FK
amf_snr(d,d5) % AMF

% S/N --- horizontal noise

amf_snr(d,dn2) % Noisy
amf_snr(d,d21) % BP
amf_snr(d,d22) % BP+SOSVMF
amf_snr(d,d23) %Bp+SOSVMF+Curvelet
amf_snr(d,d24) % BP+SOSVMF+Curvelet+FK
amf_snr(d,d25) % AMF



