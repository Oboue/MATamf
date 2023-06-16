% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

% Script to plot Figures 15 and 16

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
%-------------------------------------------------------------------------
clc;clear;close all;
%% addpath
addpath('../amf_data/');
addpath('../amfsrcs/');
addpath('../seistr/');
%% read field reflection seismic data 1

fid=fopen('real2.bin','r');
d2d=fread(fid,[256,256],'float');

dn=amf_scale(d2d,2);
[n1,n2]=size(dn);
din=dn;
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
d_bp=amf_bp(din,dt,flo,fhi,nplo,nphi,phase,verb0);
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
rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                   % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=5;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
%                                                                                                             
tic
d_sosvmf=amf_sosvmf(din,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
% 
%% Denoising using curvelet method

din=dn;
d_est=din;
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=2.0;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
tic
d_ct=amf_ct(din,d_est,n1,n2,c1,c2,c3,niter1);
toc
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

par.niter=niter;                 % number of nonlinear iterations
par.liter=liter;                 % liter: number of linear iterations (in divn)
par.order1=order1;               % order: accuracy order
par.eps_dv=eps_dv;               % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=eps_cg;               % eps_cg: eps for CG    (default: 1)
par.tol_cg=tol_cg;               % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=50;             % rect:  smoothing radius (ndim*1)
par.rect(2)=50;             % "      "        "
par.rect(3)=1;             % "      "        "
par.verb1=verb1;                 % verbosity flag

par.adj=adj;                     % adjoint flag
par.add=add;                     % adding flag
par.ns=5;                       % spray radius
par.order2=order2;               % PWD order
par.eps=eps;                     % regularization (default:0.01);
par.ndn=ndn;                     % size of dn (n1*n2)
par.nds=nds;                     % size of ds (n1*n2)
par.type_mf=type_mf;             % 0 (MF) or 1 (SVMF)
par.ifsmooth=ifsmooth;                   % 1 (if smooth) or 0 (only MF)

par.w=0.0;                        % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

par.c1=c1;                      % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c2;                      % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=0.5;                      % Thresholding parameter (alpha)
par.niter1=niter1;              % Number of iteration
%
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
par.rec(1) = 10;
par.rec(2) = 10;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0)
%
tic
d_amf=amf(din,par);
toc 
%
%% Noise attenuation using drr method
flow=0;
fhigh=250;
dt=0.004; 
N=10;% rank
K=2;% damping factor
verb=0;

tic
d_drr=amf_fxydrr(din,flow,fhigh,dt,N,K,verb);    
toc
dn_drr=din-d_drr;

% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% subplot(1,2,1);amf_imagesc(d_drr,100,2,x,t);caxis([-0.4,0.4]);
% subplot(1,2,2);amf_imagesc(dn_drr,100,2,x,t);caxis([-0.4,0.4]);
% % amf_snr(dcl,d_drr1) % drr
%%
%% Plot figures

t=[0:n1]*0.002;
x=1:n2;
ngap=5;
d_d1=[din,zeros(n1,ngap),d_bp,din-d_bp]; % BP
d_d2=[din,zeros(n1,ngap),d_sosvmf,din-d_sosvmf]; % SOSVMF
% d_d3=[din,zeros(n1,ngap),d_fk,din-d_fk]; % FK
d_d4=[din,zeros(n1,ngap),d_ct,din-d_ct]; % Curvelet
d_d5=[din,zeros(n1,ngap),d_amf,din-d_amf]; % AMF
d_d6=[din,zeros(n1,ngap),d_drr,din-d_drr]; % drr

%%Plot figures

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(3,2,1); imagesc(x,t,din);colormap(amf_seis);hold on
text(n1/-10,-0.05,'(a)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(n2/1.50,-0.02,'BP','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%
subplot(3,2,2); imagesc(x,t,d_bp);colormap(amf_seis);hold on
text(n1/-10,-0.05,'(b)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'BP','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%
subplot(3,2,3); imagesc(x,t,d_sosvmf);colormap(amf_seis);hold on
text(n1/-10,-0.05,'(c)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'SOSVMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% 
subplot(3,2,4); imagesc(x,t,d_ct);colormap(amf_seis);hold on
text(n1/-10,-0.05,'(d)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'Curvelet','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(n2/1.50,-0.02,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% 
subplot(3,2,5); imagesc(x,t,d_amf);colormap(amf_seis);hold on
text(n1/-10,-0.05,'(e)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% 
subplot(3,2,6); imagesc(x,t,d_drr);colormap(amf_seis);hold on
text(n1/-10,-0.05,'(f)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'DRR','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%
annotation(gcf,'textarrow',[0.157777777777778 0.182222222222222],...
    [0.769381746810597 0.720314033366043],'Color',[1 0 0],'String','',...
    'LineWidth',2,...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',25);
annotation(gcf,'textarrow',[0.258888888888889 0.283333333333333],...
    [0.811579980372914 0.76251226692836],'Color',[1 0 0],'String','',...
    'LineWidth',2,...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',25);
annotation(gcf,'textarrow',[0.384444444444444 0.408888888888889],...
    [0.875368007850829 0.826300294406275],'Color',[1 0 0],'String','',...
    'LineWidth',2,...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',25);
annotation(gcf,'textarrow',[0.6 0.624444444444444],...
    [0.773307163886159 0.724239450441605],'Color',[1 0 0],'String','',...
    'LineWidth',2,...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',25);
annotation(gcf,'textarrow',[0.7 0.724444444444444],...
    [0.816486751717366 0.767419038272812],'Color',[1 0 0],'String','',...
    'LineWidth',2,...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',25);
annotation(gcf,'textarrow',[0.825555555555556 0.85],...
    [0.922473012757599 0.873405299313045],'Color',[1 0 0],'String','',...
    'LineWidth',2,...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',25);
annotation(gcf,'textarrow',[0.36 0.384444444444444],...
    [0.61727183513248 0.568204121687926],'Color',[1 0 0],'String','',...
    'LineWidth',2,...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',25);
annotation(gcf,'textarrow',[0.808888888888889 0.833333333333333],...
    [0.613346418056916 0.564278704612362],'Color',[1 0 0],'String','',...
    'LineWidth',2,...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',25);
annotation(gcf,'textarrow',[0.356666666666667 0.381111111111111],...
    [0.318940137389594 0.269872423945041],'Color',[1 0 0],'String','',...
    'LineWidth',2,...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',25);
annotation(gcf,'textarrow',[0.807777777777778 0.832222222222222],...
    [0.317958783120704 0.26889106967615],'Color',[1 0 0],'String','',...
    'LineWidth',2,...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',25);

print(gcf,'-depsc','-r300','fig15.eps');
%%

d_dn1=din-d_bp; % BP
d_dn2=din-d_sosvmf; % SOSVMF
d_dn3=din-d_ct; % Curvelet
d_dn4=din-d_amf; % AMF
d_dn5=din-d_drr; % drr

figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
ax(1)=subplot(2,5,1); imagesc(x,t,d_dn1);colormap(ax(1),amf_seis);hold on
text(n1/-10,-0.05,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'BP','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%
ax(1)=subplot(2,5,2); imagesc(x,t,d_dn2);colormap(ax(1),amf_seis);hold on
text(n1/-10,-0.05,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'SOSVMF','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');

ax(1)=subplot(2,5,3); imagesc(x,t,d_dn3);colormap(ax(1),amf_seis);hold on
text(n1/-10,-0.05,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'Curvelet','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(n2/1.50,-0.02,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');

ax(1)=subplot(2,5,4); imagesc(x,t,d_dn4);colormap(ax(1),amf_seis);hold on
text(n1/-10,-0.05,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'AMF','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');

ax(1)=subplot(2,5,5); imagesc(x,t,d_dn5);colormap(ax(1),amf_seis);hold on
text(n1/-10,-0.05,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'DRR','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
print(gcf,'-depsc','-r300','fig15.eps');
%% local similarity metrics
rect=[5,5,1];niter=20;eps=0;verb=0;
% [simid1]=amf_localsimi(d2dssp_nmo-d_bp,d_bp,rect,niter,eps,verb);
[simid1]=amf_localsimi(din-d_bp,d_bp,rect,niter,eps,verb);
[simid2]=amf_localsimi(din-d_sosvmf,d_sosvmf,rect,niter,eps,verb);
[simid3]=amf_localsimi(din-d_ct,d_ct,rect,niter,eps,verb);
[simid4]=amf_localsimi(din-d_amf,d_amf,rect,niter,eps,verb);
[simid5]=amf_localsimi(din-d_drr,d_drr,rect,niter,eps,verb);

% figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
ax(2)=subplot(2,5,6); imagesc(x,t,simid1);colormap(ax(2),jet); hold on;
text(n1/-10,-0.05,'(f)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% title('BP')
% text(n1/-10,-0.05,'(d)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'BP','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 8; % c.Label.FontWeight = bold;
 caxis([0 0.75]);
ylabel('Time (s)','Fontsize',8,'fontweight','bold');
xlabel('Trace','Fontsize',8,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',8,'Fontweight','bold');

ax(2)=subplot(2,5,7); imagesc(x,t,simid2); colormap(ax(2),jet);hold on;
text(n1/-10,-0.05,'(g)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% title('DRR')
% text(n1/-10,-0.05,'(d)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'SOSVMF','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 8; % c.Label.FontWeight = bold;
 caxis([0 0.5]);
ylabel('Time (s)','Fontsize',8,'fontweight','bold');
xlabel('Trace','Fontsize',8,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',8,'Fontweight','bold');

ax(2)=subplot(2,5,8); imagesc(x,t,simid3); colormap(ax(2),jet); hold on;
text(n1/-10,-0.05,'(h)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% title('Curvelet')
% text(n1/-10,-0.05,'(d)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'Curvelet','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 8; % c.Label.FontWeight = bold;
 caxis([0 0.75]);
ylabel('Time (s)','Fontsize',8,'fontweight','bold');
xlabel('Trace','Fontsize',8,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',8,'Fontweight','bold');

ax(2)=subplot(2,5,9); imagesc(x,t,simid4); colormap(ax(2),jet);hold on;
text(-10,-0.05,'(i)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% title('AMF')
% text(n1/-10,-0.05,'(d)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'AMF','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 8; % c.Label.FontWeight = bold;
caxis([0 0.75]);
ylabel('Time (s)','Fontsize',8,'fontweight','bold');
xlabel('Trace','Fontsize',8,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',8,'Fontweight','bold');


ax(2)=subplot(2,5,10); imagesc(x,t,simid5); colormap(ax(2),jet); hold on;
% text(-3.5,-75,'(d)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% title('DRR')
text(n1/-10,-0.05,'(j)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(n2/5,-0.02,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'DRR','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% colormap(amf_seis);caxis([-0.4,0.4]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Trace','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');


c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 8; % c.Label.FontWeight = bold;
caxis([0 0.75]);
ylabel('Time (s)','Fontsize',8,'fontweight','bold');
xlabel('Trace','Fontsize',8,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',8,'Fontweight','bold');



% set(gca,'XTick',0:200:1200)
% set(gca,'XTickLabel',0:2:12)
