% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

%  Script to plot Figures 6 and 7

%  Copyright (C) Oboue et al., 2022

%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.

%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/

%  References

%  Oboue et al., 2022
%  Huang, G., M. Bai, Q. Zhao, W. Chen, and Y. Chen, 2021, Erratic noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, 69, 101-121.
%  Chen, Y., S. Zu, Y. Wang, and X. Chen, 2020, Deblending of simultaneous-source data using a structure-oriented space varying median filter, Geophysical Journal International, 222, 1805�1�723.
%  Zhao, Q., Q. Du, X. Gong, and Y. Chen, 2018, Signal-preserving erratic noise attenuation via iterative robust sparsity-promoting filter, IEEE Transactions on Geoscience and Remote Sensing, 56, 1558-0644.
%  -------------------------------------------------------------------------
clc;clear;close all;
%% addpath
addpath('../amf_data/');
addpath('../amfsrcs/');
addpath('../seistr/');
%% load data clean data

load micro_sf_3001_334_3.mat
d=amf_scale(data(:,:,1),2);      % Clean data
[n1,n2]=size(d);
%% load horizontal noise 
load d_noise_micro.mat

randn('state',211111);
nnn=amf_seisdither(nn,round(amf_meanf(20*randn(1,nx),20,1,2)));

dn1=d+0.1*nnn;                    % Clean data corrupted by horizontal noise
%%  Denosing using the BP method
%   Parameter tuning for the BP method

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
% save d1.mat d1 
% load d1.mat
%
%% Denosing using the BP+SOSVMF method 
%  Parameter tuning：add the key parameters of the SOSVMF method

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
% save d2.mat d2
% load d2.mat
%
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0.2;                        % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
% 
tic
d3=amf_bandpasssosvmffk(dn1,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
% save d3.mat d3
% load d3.mat
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
% save d4.mat d4
% load d4.mat
% 
%% Denoising using the AMF method 
%  Parameter tuning: add the key parameters for the local orthogonalization operation

par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=200;   % High frequency in band, default is Nyquist
par.nplo=6;    % number of poles for low cutoff
par.nphi=6;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag

par.niter=niter;                      % number of nonlinear iterations
par.liter=liter;                      % liter: number of linear iterations (in divn)
par.order1=order1;                    % order: accuracy order
par.eps_dv=eps_dv;                    % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=eps_cg;                    % eps_cg: eps for CG    (default: 1)
par.tol_cg=tol_cg;                    % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=rect(1);                  % rect:  smoothing radius (ndim*1)
par.rect(2)=rect(2);                  % "      "        "
par.rect(3)=rect(3);                  % "      "        "
par.verb1=verb1;                      % verbosity flag

par.adj=adj;                          % adjoint flag
par.add=add;                          % adding flag
par.ns=ns;                            % spray radius
par.order2=order2;                    % PWD order
par.eps=eps;                          % regularization (default:0.01);
par.ndn=ndn;                          % size of dn (n1*n2)
par.nds=nds;                          % size of ds (n1*n2)
par.type_mf=type_mf;                  % 0 (MF) or 1 (SVMF)
par.ifsmooth=ifsmooth;                % 1 (if smooth) or 0 (only MF)

par.w=w;                              % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

par.c1=c1;                            % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c2;                            % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=c3;                            % Thresholding parameter (alpha)
par.niter1=niter1;                    % Number of iteration

rec = zeros(3, 1);                    % 3-D vector denoting smooth radius 
par.rec(1) = 30;
par.rec(2) = 30;
par.rec(3) = 30;
par.eps1=0;                           % regularization parameter, default 0.0
par.niter2=20;                        % number of CG iterations
par.verb=1;                           % verbosity flag (default: 0) 

%
% tic
% d_amf1=amf(dn1,par);
% d5=d_amf1;
% toc
% save d_amf1.mat d_amf1
% load d_amf1.mat
%%
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
% tic
% d21=amf_bandpass(dn2,dt,flo,fhi,nplo,nphi,phase,verb0);
% save d21.mat d21
load d21.mat
% toc
% %
% %% Denosing using the BP+SOSVMF method 
% % Parameter tuning：add the key parameters of the SOSVMF method

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
% save d22.mat d22
% load d22.mat
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0.2;        % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
tic
d23=amf_bandpasssosvmffk(dn2,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
% save d23.mat d23
% load d23.mat
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
% save d24.mat d24
% load d24.mat
%
%% Denoising using the AMF method 
%  Parameter tuning: add the key parameters for the local orthogonalization operation
% 
par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=200;   % High frequency in band, default is Nyquist
par.nplo=6;% number of poles for low cutoff
par.nphi=6;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag

par.niter=niter;                    % number of nonlinear iterations
par.liter=liter;                    % liter: number of linear iterations (in divn)
par.order1=order1;                  % order: accuracy order
par.eps_dv=eps_dv;                  % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=eps_cg;                  % eps_cg: eps for CG    (default: 1)
par.tol_cg=tol_cg;                  % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=rect(1);                % rect:  smoothing radius (ndim*1)
par.rect(2)=rect(2);                % "      "        "
par.rect(3)=rect(3);                % "      "        "
par.verb1=verb1;                    % verbosity flag

par.adj=adj;                        % adjoint flag
par.add=add;                        % adding flag
par.ns=ns;                          % spray radius
par.order2=order2;                  % PWD order
par.eps=eps;                        % regularization (default:0.01);
par.ndn=ndn;                        % size of dn (n1*n2)
par.nds=nds;                        % size of ds (n1*n2)
par.type_mf=type_mf;                % 0 (MF) or 1 (SVMF)
par.ifsmooth=ifsmooth;              % 1 (if smooth) or 0 (only MF)

par.w=0.00;                         % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

par.c1=c1;                          % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c2;                          % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=c3;                          % Thresholding parameter (alpha)
par.niter1=niter1;                  % Number of iteration
rec = zeros(3, 1);                  % 3-D vector denoting smooth radius 
par.rec(1) = 30;
par.rec(2) = 30;
par.rec(3) = 30;
par.eps1=0;                         % regularization parameter, default 0.0
par.niter2=20;                      % number of CG iterations
par.verb=1;                         % verbosity flag (default: 0) 
%
tic
d_amf2=amf(dn2,par);
d25=d_amf2;
toc
% load d_amf2.mat
% save d_amf2.mat d_amf2
% save d25.mat d25
%% Plot figures

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

print(gcf,'-depsc','-r300','fig6.eps');
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

print(gcf,'-depsc','-r300','fig7.eps');
%% Numerical test using signal-to-noise ratio (S/N)

%S/N --- horizontal noise

amf_snr(d,dn1) % Noisy
amf_snr(d,d1)  % BP
amf_snr(d,d2)  % BP+SOSVMF
amf_snr(d,d3)  % Bp+SOSVMF+Curvelet
amf_snr(d,d4)  % BP+SOSVMF+Curvelet+FK
amf_snr(d,d_amf)  % AMF

% S/N --- Vertical noise

amf_snr(d,dn2) % Noisy
amf_snr(d,d21) % BP
amf_snr(d,d22) % BP+SOSVMF
amf_snr(d,d23) %Bp+SOSVMF+Curvelet
amf_snr(d,d24) % BP+SOSVMF+Curvelet+FK
amf_snr(d,d_amf2) % AMF

