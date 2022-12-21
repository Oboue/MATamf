% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

% Script to plot Figures 13 and 14

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
%% Generate the clean synthetic data (dc) including linear events

dcl=amf_levents(200);dcl=amf_scale(dcl);

[n1,n2]=size(dcl);

mask=rand(1,n2);
mask(logical(mask<0.9))=0;
mask(logical(mask>=0.9))=1;

err_n=zeros(size(dcl));
for i=1:n1
    randn('state',123456+i);
    err_n(i,:)=0.5*randn(1,n2).*mask;
end

randn('state',201920);
ran_n=0.1*randn(n1,n2);

dnl=dcl+err_n+ran_n;
dinl=dnl;

dt=0.004;
t=[0:n1-1]*dt; 
x=[1:n2];
%%  Denosing using the BP method
% Parameter tuning for the BP method

dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=251;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=7;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb=0;    % verbosity flag

tic
d_bp1=amf_bandpass(dinl,dt,flo,fhi,nplo,nphi,phase,verb);
toc
%% Denoising using the SOSVMF method

niter=5;                      % number of nonlinear iterations
liter=20;                     % liter: number of linear iterations (in divn)
order1=2;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=11;                   % rect:  smoothing radius (ndim*1)
rect(2)=5;                    % "      "        "
rect(3)=1;                    % "      "        "
verb=1;                       % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=4;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)

d_sosvmf1=amf_sosvmf(dinl,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
%% Denoising using the dip filter in FK domain method

w=0.0001;                     % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

d_fk1=dinl-amf_fk_dip(dinl,w);
%% Denoising using curvelet method

dinl=dnl;
d_est=dinl;

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=3;                % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

d_curvelet1=amf_curvelet(dinl,d_est,n1,n2,c1,c2,c3,niter1);
%% AMF method
par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=251;   % High frequency in band, default is Nyquist
par.nplo=6;    % number of poles for low cutoff
par.nphi=7;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag
%
par.niter=5;                      % number of nonlinear iterations
par.liter=20;                     % liter: number of linear iterations (in divn)
par.order1=2;                     % order: accuracy order
par.eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
par.tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=11;                   % rect:  smoothing radius (ndim*1)
par.rect(2)=5;                    % "      "        "
par.rect(3)=1;                    % "      "        "
par.verb1=1;                      % verbosity flag

par.adj=0;                        % adjoint flag
par.add=0;                        % adding flag
par.ns=4;                         % spray radius
par.order2=2;                     % PWD order
par.eps=0.01;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % size of ds (n1*n2)
par.type_mf=1;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
% 
par.w=0.0;               % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=0.9;              % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration

% 
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
par.rec(1) = 5;
par.rec(2) = 5;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 

d_amf1=amf(dinl,par);
%% Plot figures/linear events

t=[0:n1]*0.004;
x=1:n2;
ngap=5;
d1_d1=[dcl,zeros(n1,ngap),dnl]; % Noisy
d1_d2=[d_bp1,zeros(n1,ngap),dnl-d_bp1]; % BP
d1_d3=[d_sosvmf1,zeros(n1,ngap),dnl-d_sosvmf1]; % SOSVMF
d1_d4=[d_fk1,zeros(n1,ngap),dnl-d_fk1]; % FK
d1_d5=[d_curvelet1,zeros(n1,ngap),dnl-d_curvelet1]; % Curvelet
d1_d6=[d_amf1,zeros(n1,ngap),dnl-d_amf1];% AMF

%% Numerical test using signal-to-noise ratio (S/N)----linear events
amf_snr(dcl,dinl) %Noisy
amf_snr(dcl,d_bp1) % BP
amf_snr(dcl,d_sosvmf1) % SOSVMF
amf_snr(dcl,d_fk1) % FK
amf_snr(dcl,d_curvelet1) % Curvelet 
amf_snr(dcl,d_amf1) % AMF
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
par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=251;   % High frequency in band, default is Nyquist
par.nplo=6;    % number of poles for low cutoff
par.nphi=7;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag
%
par.niter=2;                      % number of nonlinear iterations
par.liter=10;                     % liter: number of linear iterations (in divn)
par.order1=3;                     % order: accuracy order
par.eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
par.tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=11;                   % rect:  smoothing radius (ndim*1)
par.rect(2)=5;                    % "      "        "
par.rect(3)=1;                    % "      "        "
par.verb1=1;                      % verbosity flag

par.adj=0;                        % adjoint flag
par.add=0;                        % adding flag
par.ns=9;                         % spray radius
par.order2=2;                     % PWD order
par.eps=0.01;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % size of ds (n1*n2)
par.type_mf=1;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
% 
par.w=0.0;               % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=0.4;              % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration

% 
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
par.rec(1) = 8;
par.rec(2) = 8;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 

dout=amf(din,par);
d_amf=dout;
%%
t=[0:n1]*0.004;
x=1:n2;
ngap=5;
d_d1=[dc,zeros(n1,ngap),dn]; %Noisy
d_d2=[d_bp,zeros(n1,ngap),dn-d_bp]; % BP
d_d3=[d_sosvmf,zeros(n1,ngap),dn-d_sosvmf]; % SOSVMF
d_d4=[d_fk,zeros(n1,ngap),dn-d_fk]; % FK
d_d5=[d_curvelet,zeros(n1,ngap),dn-d_curvelet]; % Curvelet
d_d6=[d_amf,zeros(n1,ngap),dn-d_amf]; % AMF

% Numerical test using signal-to-noise ratio (S/N)----Curved events

amf_snr(dc,dn) %Noisy
amf_snr(dc,d_bp) % BP
amf_snr(dc,d_sosvmf) % SOSVMF
amf_snr(dc,d_fk) % FK
amf_snr(dc,d_curvelet) % Curvelet 
amf_snr(dc,d_amf) % AMF
%% Plot figures
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(3,2,1); imagesc(x,t,d1_d1);colormap(amf_seis);hold on
text(n1/-40,-0.2,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'Clean','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.25,-0.08,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
%
subplot(3,2,2); imagesc(x,t,d1_d2);colormap(amf_seis);hold on
text(n1/-40,-0.2,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'BP','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
annotation(gcf,'textarrow',[0.654444444444444 0.632222222222222],...
    [0.898477157360408 0.911675126903556],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.701111111111111 0.678888888888889],...
    [0.758375634517767 0.771573604060915],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
%
subplot(3,2,3); imagesc(x,t,d1_d3);colormap(amf_seis);hold on
text(n1/-40,-0.2,'(c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'SOSVMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
annotation(gcf,'textbox',...
    [0.149888888888888 0.453807106598985 0.133444444444444 0.0426395939086297],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textarrow',[0.38 0.376666666666666],...
    [0.555329949238579 0.584261113095579],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.714444444444444 0.683333333333333],...
    [0.608121827411168 0.607611366902685],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

subplot(3,2,4); imagesc(x,t,d1_d5);colormap(amf_seis);hold on
text(n1/-40,-0.2,'(d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'Curvelet','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
annotation(gcf,'textarrow',[0.714444444444444 0.683333333333333],...
    [0.608121827411168 0.607611366902685],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.634444444444444 0.611111111111111],...
    [0.441624365482235 0.453296646090502],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.771111111111111 0.752222222222222],...
    [0.580710659898478 0.601015228426396],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.875555555555556 0.863333333333333],...
    [0.56548223350254 0.587306798374768],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
subplot(3,2,5); imagesc(x,t,d1_d6);colormap(amf_seis);hold on
text(n1/-40,-0.2,'(e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
annotation(gcf,'textarrow',[0.272222222222222 0.242222222222222],...
    [0.226395939086294 0.206596138476287],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Improved','signal','energy'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.431111111111111 0.418888888888888],...
    [0.265989847715736 0.287814412587964],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

print(gcf,'-depsc','-r300','fig13.eps');

%% crueved events

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(3,2,1); imagesc(x,t,d_d1);colormap(amf_seis);hold on
text(n1/-40,-0.2,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'Clean','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.25,-0.08,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
%
subplot(3,2,2); imagesc(x,t,d_d2);colormap(amf_seis);hold on
text(n1/-40,-0.2,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'BP','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
annotation(gcf,'textarrow',[0.728888888888889 0.694444444444444],...
    [0.88832487309645 0.877662128324002],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.638888888888889 0.604444444444444],...
    [0.748223350253811 0.737560605481363],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
%
subplot(3,2,3); imagesc(x,t,d_d3);colormap(amf_seis);hold on
text(n1/-40,-0.2,'(c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'SOSVMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
annotation(gcf,'textarrow',[0.268888888888889 0.252222222222222],...
    [0.551269035532996 0.532484463349385],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textbox',...
    [0.165444444444444 0.455837563451777 0.133444444444444 0.0426395939086297],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textarrow',[0.391111111111111 0.416666666666667],...
    [0.502538071065992 0.52131695065903],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

subplot(3,2,4); imagesc(x,t,d_d5);colormap(amf_seis);hold on
text(n1/-40,-0.2,'(d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'Curvelet','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
annotation(gcf,'textarrow',[0.622222222222222 0.605555555555555],...
    [0.569543147208123 0.550758575024513],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.708888888888889 0.692222222222222],...
    [0.55532994923858 0.53654537705497],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.663333333333333 0.646666666666667],...
    [0.473096446700508 0.501522842639594],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.827777777777778 0.853333333333333],...
    [0.508629441624367 0.527408321217405],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

subplot(3,2,5); imagesc(x,t,d_d6);colormap(amf_seis);hold on
text(n1/-40,-0.2,'(e)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,-0.08,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.3,-0.08,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
annotation(gcf,'textarrow',[0.283333333333333 0.253333333333333],...
    [0.250761421319797 0.23096162070979],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Improved','signal','energy'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
print(gcf,'-depsc','-r300','fig14.eps');