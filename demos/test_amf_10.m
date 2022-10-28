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
%% AMF method 
%% Load the clean synthetic data (dc) including linear events

load Curveddata
dc=data;
[n1,n2]=size(dc);
%% Load the noisy data corrupted by both random and erratic noise

load dncurved 
% dn=dn;
[n1,n2]=size(dn);

din=dn;
%%  Denosing using the BP method
% Parameter tuning for the BP method

dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=251;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=7;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb=0;    % verbosity flag

%
tic
d_bp=amf_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb);
toc
%
%% Denoising using the SOSVMF method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=27;                   % rect:  smoothing radius (ndim*1)
rect(2)=7;                    % "      "        "
rect(3)=1;                    % "      "        "
verb=1;                       % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=9;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)

d_sosvmf=amf_sosvmf(din,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
%% Denoising using the dip filter in FK domain method

w=0.0001;                    % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
d_fk=din-amf_fk_dip(din,w);
%% Denoising using curvelet method

din=dn;
d_est=din;

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=2.9;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

dout=amf_curvelet(din,d_est,n1,n2,c1,c2,c3,niter1);
d_curvelet=dout;

%% AMF method
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=251;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=7;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=11;                   % rect:  smoothing radius (ndim*1)
rect(2)=5;                    % "      "        "
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
w=0.0;               % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.4;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

% 
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
rec(1) = 8;
rec(2) = 8;
rec(3) = 1;
eps1=0;               % regularization parameter, default 0.0
niter2=20;            % number of CG iterations
verb=1;               % verbosity flag (default: 0) 

dout=amf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1,rec,eps1,niter2,verb);
d_amf=dout;
%% Plot figure

t=[0:n1]*0.004;
x=1:n2;
ngap=5;
d_d1=[dc,zeros(n1,ngap),dn]; %Noisy
d_d2=[d_bp,zeros(n1,ngap),dn-d_bp]; % BP
d_d3=[d_sosvmf,zeros(n1,ngap),dn-d_sosvmf]; % SOSVMF
d_d4=[d_fk,zeros(n1,ngap),dn-d_fk]; % FK
d_d5=[d_curvelet,zeros(n1,ngap),dn-d_curvelet]; % Curvelet
d_d6=[d_amf,zeros(n1,ngap),dn-d_amf]; % AMF

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(3,2,1); imagesc(x,t,d_d1);colormap(amf_seis);hold on
text(n2/-500,-0.3,'(a)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'Clean','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.25,-0.08,'Raw','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');

colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
%
subplot(3,2,2); imagesc(x,t,d_d2);colormap(amf_seis);hold on
text(n2/-500,-0.3,'(b)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'BP','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
%
subplot(3,2,3); imagesc(x,t,d_d3);colormap(amf_seis);hold on
text(n2/-500,-0.3,'(c)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'SOSVMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');

subplot(3,2,4); imagesc(x,t,d_d4);colormap(amf_seis);hold on
text(n2/-500,-0.3,'(d)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');

subplot(3,2,5); imagesc(x,t,d_d5);colormap(amf_seis);hold on
text(n2/-500,-0.3,'(e)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');

subplot(3,2,6); imagesc(x,t,d_d6);colormap(amf_seis);hold on
text(n2/-500,-0.3,'(f)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'AMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');

print(gcf,'-depsc','-r300','fig18-amf.pdf');
print(gcf,'-depsc','-r300','fig18-amf.eps');
%% Numerical test using signal-to-noise ratio (S/N)

amf_snr(dc,dn) %Noisy
amf_snr(dc,d_bp) % BP
amf_snr(dc,d_sosvmf) % SOSVMF
amf_snr(dc,d_fk) % FK
amf_snr(dc,d_curvelet) % Curvelet 
amf_snr(dc,d_amf) % AMF


