% Demo for the advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

% This script to plot Figure 11

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
d_bp=amf_bp(d2dssp_nmo,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%
dn_bp=d2dssp_nmo-d_bp;
%% Denosing using the BP+SOSVMF method 
% Parameter tuning：add the key parameters of the SOSVMF method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=2;                   % rect:  smoothing radius (ndim*1)
rect(2)=2;                    % "      "        "
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
d_bpsosvmf=amf_bpsosvmf(d2dssp_nmo,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%
dn_bpsosvmf=d2dssp_nmo-d_bpsosvmf;
%% Denoising data using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the key parameters of the curvelet method

w=0.0;
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=1;               % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
tic
d_bpsosvmffkct=amf_bpsosvmffkct(d2dssp_nmo,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
dn_bpsosvmffkct=d2dssp_nmo-d_bpsosvmffkct;
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
par.rect(1)=200;                   % rect:  smoothing radius (ndim*1)
par.rect(2)=200;                    % "      "        "
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
par.c3=1;             % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration

rec = zeros(3, 1);       % 3-D vector denoting smooth radius 
par.rec(1) = 1000;
par.rec(2) = 1000;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 
%
tic
d_amf=amf(d2dssp_nmo,par);
d5l=d_amf;
toc
dn_amf=d2dssp_nmo-d_amf;
%% Noise attenuation  using drr method 

flow=1/75;
fhigh=1/15;
dt=1; 
N=6;% rank
K=2;% damping factor
verb=0;

tic
d_drr=amf_fxydrr(d2dssp_nmo,flow,fhigh,dt,N,K,verb);    
toc
dn_drr=d2dssp_nmo-d_drr;
% amf_snr(dcl,d_drr1) % drr
%%
x=1:nh;
t0=-300;
t1=-100;
%t0=-500;
%t1=100;
keep = t>=t0 & t<=t1;
%%
figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(2,2,1); imagesc(1:nh,t,reshape(d2dssp_nmo,n1,nh)); hold on;
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
annotation(gcf,'textarrow',[0.171111111111111 0.152222222222222],...
    [0.676683143725999 0.699316888581966],'Color',[0 0 1],'TextColor',[0 0 1],...
    'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',20);
annotation(gcf,'textbox',...
    [0.389888888888888 0.829683174162464 0.0723333333333335 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.388777777777778 0.702258710516452 0.0723333333333337 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
subplot(2,2,2); imagesc(1:nh,t,reshape(d_bpsosvmf,n1,nh)); hold on;
text(-3.5,-75,'(b)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
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
annotation(gcf,'textarrow',[0.616666666666667 0.597777777777777],...
    [0.670795018112652 0.693428762968619],'Color',[0 0 1],'TextColor',[0 0 1],...
    'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',20);
annotation(gcf,'textbox',...
    [0.827666666666666 0.829683174162461 0.0723333333333336 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.827666666666667 0.702258710516448 0.0723333333333337 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
subplot(2,2,3); imagesc(1:nh,t,reshape(d_bpsosvmffkct,n1,nh)); hold on;
text(-3.5,-75,'(c)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
plot(1:nh,t410_ref*ones(1,nh),'--r','linewidth',4)
plot(1:nx*ny,t660_ref*ones(1,nx*ny),'--r','linewidth',4)
title('SOSVMF+Curvelet')
colormap(amf_seis);caxis([-0.02 0.02]); 
axis xy
xlim([1 35])
ylim([-300 -100])
ylabel('Time to SS (sec)')
xlabel('Trace');
set(gca,'fontsize',16)
annotation(gcf,'textarrow',[0.176666666666667 0.157777777777778],...
    [0.199558188654254 0.222191933510221],'Color',[0 0 1],'TextColor',[0 0 1],...
    'String',{'Weaker signal energy'},...
    'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',20);
annotation(gcf,'textbox',...
    [0.389888888888888 0.355668118751846 0.0723333333333331 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.388777777777778 0.230206363643623 0.0723333333333337 0.0270961896365999],...
    'Color',[1 0 0],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
subplot(2,2,4); imagesc(1:nh,t,reshape(d_amf,n1,nh)); hold on;
text(-3.5,-75,'(d)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
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
annotation(gcf,'textarrow',[0.616666666666667 0.597777777777778],...
    [0.200561754645278 0.223195499501246],'Color',[0 0 1],'TextColor',[0 0 1],...
    'String',{'Improved signal energy'},...
    'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',20);
annotation(gcf,'textbox',...
    [0.828777777777778 0.354398541532034 0.0723333333333336 0.0270961896366002],...
    'Color',[1 0 0],...
    'String','S410S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
annotation(gcf,'textbox',...
    [0.828777777777777 0.228762256955259 0.0723333333333337 0.0270961896366001],...
    'Color',[1 0 0],...
    'String','S660S',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);
print(gcf,'-depsc','-r300','fig11.eps');
















