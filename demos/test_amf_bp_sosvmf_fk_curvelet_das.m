% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

% Script to plot Figure 26

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
%% Load the raw DAS seismic data corrupted by strong hogh-frequency noise, high-amplitude 

NOs=[1,20,10,25,11,2];
labels={...                                             %P-arrival sample NO from the SEGY file
    'FORGE\_78-32\_iDASv3-P11\_UTC190423150554.sgy',... %24169
    'FORGE\_78-32\_iDASv3-P11\_UTC190426070723.sgy',... %24811
    'FORGE\_78-32\_iDASv3-P11\_UTC190426062208.sgy',... %26090
    'FORGE\_78-32\_iDASv3-P11\_UTC190426110008.sgy',... %4921
    'FORGE\_78-32\_iDASv3-P11\_UTC190426062553.sgy',... %8934
    'FORGE\_78-32\_iDASv3-P11\_UTC190423182409.sgy'};   %4210
eq=zeros(2000,960);
[n1,n2]=size(eq);
t=[0:n1]*0.0005;

for ii=4
    
if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('amf_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
din=d1;
%% Denosing using the BP method
%  Parameter tuning for the BP method
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
d_bp=amf_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%
%% 
%%  Denosing using the SOSVMF method
%   Parameter tuning for the BP method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=500;                   % rect:  smoothing radius (ndim*1)
rect(2)=500;                   % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=20;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
tic
d_sosvmf=amf_sosvmf(din,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%
%% Denoising using the dip filter in FK domain method

w=0.2;                       % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
tic
d_fk=din-amf_fk_dip(din,w);
toc
%
%% Denoising using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the key parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=3;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
d_est=din;
% 
tic
d_curvelet=amf_curvelet(din,d_est,n1,n2,c1,c2,c3,niter1);
toc
%
%% Denoising using the AMF method 
%  Parameter tuning: add the key parameters for the local orthogonalization operation
par.dt=dt;                           % sampling
par.flo=flo;                         % Low frequency in band, default is 0
par.fhi=fhi;                         % High frequency in band, default is Nyquist
par.nplo=nplo;                       % number of poles for low cutoff
par.nphi=nphi;                       % number of poles for high cutoff
par.phase=phase;                     % y: minimum phase, n: zero phase
par.verb0=verb0;                     % verbosity flag

par.niter=niter;                     % number of nonlinear iterations
par.liter=liter;                     % liter: number of linear iterations (in divn)
par.order1=order1;                   % order: accuracy order
par.eps_dv=eps_dv;                   % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=eps_cg;                   % eps_cg: eps for CG    (default: 1)
par.tol_cg=tol_cg;                   % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=50;                      % rect:  smoothing radius (ndim*1)
par.rect(2)=50;                      % "      "        "
par.rect(3)=1;                       % "      "        "
par.verb1=verb1;                     % verbosity flag

par.adj=adj;                         % adjoint flag
par.add=add;                         % adding flag
par.ns=15;                           % spray radius
par.order2=order2;                   % PWD order
par.eps=eps;                         % regularization (default:0.01);
par.ndn=ndn;                         % size of dn (n1*n2)
par.nds=nds;                         % size of ds (n1*n2)
par.type_mf=type_mf;                 % 0 (MF) or 1 (SVMF)
par.ifsmooth=ifsmooth;               % 1 (if smooth) or 0 (only MF)

par.w=0.05;                             % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

par.c1=c1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=0.25;                % Thresholding parameter (alpha)
par.niter1=niter1;        % Number of iteration
%
rec = zeros(3, 1);        % 3-D vector denoting smooth radius 
par.rec(1) = 150;
par.rec(2) = 150;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 
tic
d_amf=amf(din,par);
toc
%
end
%% Plot figures

t=[0:n1]*0.004;
x=1:n2;
ngap=5;
d_d1=[din,zeros(n1,ngap),d_bp,din-d_bp]; % BP
d_d2=[din,zeros(n1,ngap),d_sosvmf,zeros(n1,ngap),din-d_sosvmf]; % SOSVMF
d_d3=[din,zeros(n1,ngap),d_fk,zeros(n1,ngap),din-d_fk]; % FK
d_d4=[din,zeros(n1,ngap),d_curvelet,zeros(n1,ngap),din-d_curvelet]; % Curvelet
d_d5=[din,d_amf,zeros(n1,ngap),din-d_amf]; % AMF

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(3,2,1);amf_imagesc(d_d1,100,2,x,t);caxis([-25,25]);

text(n1/-20,-0.5,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/5,-0.5,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.5,-0.5,'BP','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);;caxis([-25,25]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
%
subplot(3,2,2);amf_imagesc(d_d2,100,2,x,t);caxis([-25,25]);;hold on
text(n1/-20,-0.5,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/5,-0.5,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.5,-0.5,'SOSVMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-25,25]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');

subplot(3,2,3);amf_imagesc(d_d3,100,2,x,t);caxis([-25,25]); hold on
text(n1/-20,-0.5,'(c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/5,-0.5,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.5,-0.5,'FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-25,25]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');

subplot(3,2,4);amf_imagesc(d_d4,100,2,x,t);caxis([-25,25]); hold on
text(n1/-20,-0.5,'(d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/5,-0.5,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.5,-0.5,'Curvelet','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-25,25]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');

subplot(3,2,5);amf_imagesc(d_d5,100,2,x,t);caxis([-25,25]); hold on
text(n1/-20,-0.5,'(e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/5,-0.5,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.5,-0.5,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-25,25]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');

print(gcf,'-depsc','-r300','fig26.eps');
