% Demo for the advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

% This script to plot Figure 16

%  Copyright (C) Oboue et al., 2022

%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Genenral Public License as published
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
%-------------------------------------------------------------------------
clc;clear;close all;
%% addpath
addpath('../amf_data/');
addpath('../amfsrcs/');
addpath('../seistr/');
%% load ss precursor data

load d2dssp_nmo.mat
[n1,n2]=size(d2dssp_nmo);
%%  Denosing using the BP method
% Parameter tuning for the BP method
%
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
tic
d2_bp=amf_bp(d2dssp_nmo,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%
%% Denosing using the SOSVMF method 
% Parameter tuning：add the key parameters of the SOSVMF method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=2;                   % rect:  smoothing radius (ndim*1)
rect(2)=2;                   % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=2;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
%                                                                                                             
tic
d_sosvmf=amf_sosvmf(d2dssp_nmo,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
% 
%% Denoising using curvelet method

d_est=d2dssp_nmo;
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=2.5;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
tic
d_ct=amf_ct(d2dssp_nmo,d_est,n1,n2,c1,c2,c3,niter1);
toc
%
%% AMF
% dt=0.0005; % sampling
par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=200;   % High frequency in band, default is Nyquist
par.nplo=6;    % number of poles for low cutoff
par.nphi=6;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag

par.niter=5;                      % number of nonlinear iterations
par.liter=20;                     % liter: number of linear iterations (in divn)
par.order1=2;                     % order: accuracy order
par.eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
par.tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=2;                   % rect:  smoothing radius (ndim*1)
par.rect(2)=2;                    % "      "        "
par.rect(3)=1;                    % "      "        "
par.verb1=1;                      % verbosity flag

par.adj=0;                        % adjoint flag
par.add=0;                        % adding flag
par.ns=1;                         % spray radius
par.order2=2;                     % PWD order
par.eps=0.01;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % size of ds (n1*n2)
par.type_mf=1;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)

par.w=0; 

par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=2.5;             % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration

rec = zeros(3, 1);       % 3-D vector denoting smooth radius 
par.rec(1) = 100;
par.rec(2) = 100;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 
%
tic
d_amf=amf(d2dssp_nmo,par);
d5l=d_amf;
toc
%%
x=1:nh;
t0=-300;
t1=-100;
%t0=-500;
%t1=100;
keep = t>=t0 & t<=t1;
%%
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(3,2,1); imagesc(1:nh,t,reshape(d2dssp_nmo,n1,nh)); hold on;
text(-3.5,-75,'(a)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
plot(1:nh,t410_ref*ones(1,nh),'--r','linewidth',4)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',4)
title('Raw')
colormap(amf_seis);caxis([-0.02 0.02]); 
axis xy
xlim([1 35])
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.386555555555555 0.869918699186999 0.0723333333333335 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.382111111111111 0.788617886178868 0.0723333333333337 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

subplot(3,2,2); imagesc(1:nh,t,reshape(d2_bp,n1,nh)); hold on;
text(-3.5,-75,'(b)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
plot(1:nh,t410_ref*ones(1,nh),'--r','linewidth',4)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',4)
title('BP')
colormap(amf_seis);caxis([-0.02 0.02]); 
axis xy
xlim([1 35])
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.828777777777778 0.870934959349596 0.0723333333333336 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.825444444444444 0.790650406504075 0.0723333333333337 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

subplot(3,2,3); imagesc(1:nh,t,reshape(d_sosvmf,n1,nh)); hold on;
text(-3.5,-75,'(c)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
plot(1:nh,t410_ref*ones(1,nh),'--r','linewidth',4)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',4)
title('SOSVMF')
colormap(amf_seis);caxis([-0.02 0.02]); 
axis xy
xlim([1 35])
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.384333333333333 0.574186991869923 0.0723333333333336 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.383222222222222 0.490853658536594 0.0723333333333337 0.0270961896366003],...
    'Color',[1 0 0],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
%
subplot(3,2,4); imagesc(1:nh,t,reshape(d_ct,n1,nh)); hold on;
text(-3.5,-75,'(d)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
plot(1:nh,t410_ref*ones(1,nh),'--r','linewidth',4)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',4)
title('Curvelet')
colormap(amf_seis);caxis([-0.02 0.02]); 
axis xy
xlim([1 35])
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.829888888888889 0.575203252032525 0.0723333333333336 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.828777777777777 0.488821138211394 0.0723333333333337 0.0270961896366007],...
    'Color',[1 0 0],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

subplot(3,2,5); imagesc(1:nh,t,reshape(d_amf,n1,nh)); hold on;
text(-3.5,-75,'(e)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
plot(1:nh,t410_ref*ones(1,nh),'--r','linewidth',4)
plot(1:nh,t660_ref*ones(1,nh),'--r','linewidth',4)
title('AMF')
colormap(amf_seis);caxis([-0.02 0.02]);
axis xy
xlim([1 35])
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textbox',...
    [0.386555555555555 0.269308943089435 0.0723333333333332 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.384333333333333 0.188008130081308 0.0723333333333337 0.0270961896366003],...
    'Color',[1 0 0],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);


print(gcf,'-depsc','-r300','fig17.eps');

