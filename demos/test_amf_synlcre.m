% Demo for the advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

% This script to plot Figures 2, 3, 4 and 5

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
%% Generate the clean synthetic data including linear events

dcl=amf_levents(200);
dcl=amf_scale(dcl); % Clean data
[n1,n2]=size(dcl);
%% %Generate the noisy synthetic data having linear events

% Create the mask operator
mask=rand(1,n2);
mask(logical(mask<0.9))=0;
mask(logical(mask>=0.9))=1;

% Adding erratic noise 

err_n=zeros(size(dcl));
for i=1:n1
    randn('state',123456+i);
    err_n(i,:)=0.5*randn(1,n2).*mask;
end

% Adding random noise 

randn('state',201920);
ran_n=0.1*randn(n1,n2);

dnl=dcl+err_n+ran_n;
dinl=dnl;            % Noisy data

%%  Denosing using the bandpass method
% Parameter tuning for the BP method
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=251;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=7;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
tic
d1l=amf_bandpass(dinl,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%
%% Denosing using the BP+SOSVMF method 
% Parameter tuning：add the key parameters of the SOSVMF method

niter=5;                      % number of nonlinear iterations
liter=20;                     % liter: number of linear iterations (in divn)
order1=2;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=11;                   % rect:  smoothing radius (ndim*1)
rect(2)=5;                    % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=4;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)

%                                                                                                            
tic
d2l=amf_bandpasssosvmf(dinl,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
% 
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method

w=0; % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

% 
tic
d3l = amf_bandpasssosvmffk(dinl,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% Denoising data using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the key parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.9;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
tic
d4l=amf_bandpasssosvmffkcurvelet(dinl,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
%% AMF
% dt=0.0005; % sampling
par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=251;   % High frequency in band, default is Nyquist
par.nplo=6;    % number of poles for low cutoff
par.nphi=7;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag

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

par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=0.9;              % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration

rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
par.rec(1) = 5;
par.rec(2) = 5;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 
%
tic
d_amf=amf(dinl,par);
d5l=d_amf;
toc
%%
ngap=5;
d_1l=[dcl,dinl,zeros(n1,ngap),d1l,dinl-d1l,zeros(n1,ngap),d2l,dinl-d2l];
d_2l=[d3l,dinl-d3l,zeros(n1,ngap),d4l,d4l-dinl,zeros(n1,ngap),d5l,dinl-d5l];
%% Numerical test using signal-to-noise ratio (S/N) 
amf_snr(dcl,dnl) % Noisy
amf_snr(dcl,d1l) % BP
amf_snr(dcl,d2l) % BP+SOSVMF
amf_snr(dcl,d3l) % BP+SOSVMF+FK
amf_snr(dcl,d4l) % BP+SOSVMF+FK+curvelet
amf_snr(dcl,d5l) % AMF
%% local similarity metrics
rect=[10,10,1];niter=20;eps=0;verb=0;
[simid1]=amf_localsimi(dnl-d1l,d1l,rect,niter,eps,verb);
[simid2]=amf_localsimi(dnl-d2l,d2l,rect,niter,eps,verb);
[simid3]=amf_localsimi(dnl-d3l,d3l,rect,niter,eps,verb);
[simid4]=amf_localsimi(dnl-d4l,d4l,rect,niter,eps,verb);
[simid5]=amf_localsimi(dnl-d5l,d5l,rect,niter,eps,verb);
d_siml=[simid1,zeros(n1,ngap),simid2,zeros(n1,ngap),simid3,zeros(n1,ngap),simid4,zeros(n1,ngap),simid5];
%% Generate the clean synthetic data (dc) including curved events

load Curveddata
dc=data;         % Clean data
% [n1,n2]=size(dc);
%% Load the noisy data corrupted by both random and erratic noise

load dncurved 

[n1,n2]=size(dn);
din=dn;          % Noisy data
%%  Denosing using the BP method
% Parameter tuning for the BP method

dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=251;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=7;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag

%
tic
d1=amf_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);
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
d2=amf_bandpasssosvmf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method

w=0; % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

%
tic
d3=amf_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% Denoising data using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the key parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.4;              % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

% 
tic
d4=amf_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
%% Denoising using the AMF method 
% Parameter tuning: add the key parameters for the local orthogonalization operation

par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=251;   % High frequency in band, default is Nyquist
par.nplo=6;    % number of poles for low cutoff
par.nphi=7;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag

par.niter=2;                      % number of nonlinear iterations
par.liter=10;                     % liter: number of linear iterations (in divn)
par.order1=3;                     % order: accuracy order
par.eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
par.tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=21;                   % rect:  smoothing radius (ndim*1)
par.rect(2)=7;                    % "      "        "
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

par.w=0;                          % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=0.4;              % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration

rec = zeros(3, 1);       % 3-D vector denoting smooth radius 
par.rec(1) = 8;
par.rec(2) = 8;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 
%
tic
d_amf=amf(din,par);
d5=d_amf;
toc
%%
ngap=5;
d_1=[dc,din,zeros(n1,ngap),d1,din-d1,zeros(n1,ngap),d2,din-d2];
d_2=[d3,din-d3,zeros(n1,ngap),d4,d4-din,zeros(n1,ngap),d5,din-d5];
%% Numerical test using signal-to-noise ratio (S/N) 
amf_snr(dc,dn) % Noisy
amf_snr(dc,d1) % BP
amf_snr(dc,d2) % BP+SOSVMF
amf_snr(dc,d3) % BP+SOSVMF+FK
amf_snr(dc,d4) % BP+SOSVMF+FK+curvelet
amf_snr(dc,d5) % AMF
%% local similarity metrics
rect=[10,10,1];niter=20;eps=0;verb=0;
[simid1]=amf_localsimi(din-d1,d1,rect,niter,eps,verb);
[simid2]=amf_localsimi(din-d2,d2,rect,niter,eps,verb);
[simid3]=amf_localsimi(din-d3,d3,rect,niter,eps,verb);
[simid4]=amf_localsimi(din-d4,d4,rect,niter,eps,verb);
[simid5]=amf_localsimi(din-d5,d5,rect,niter,eps,verb);
d_sim=[simid1,zeros(n1,ngap),simid2,zeros(n1,ngap),simid3,zeros(n1,ngap),simid4,zeros(n1,ngap),simid5];
%% Plot figures linear events

dt=0.004;
t=[0:n1-1]*dt; 
x=[1:n2];

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1); imagesc(x,t,d_1l);colormap(amf_seis);hold on
text(n1/-70,-0.1,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/9.5,-0.07,'Clean','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/4,-0.07,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.07,'BP','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.17,-0.07,'BP+SOSVMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
annotation(gcf,'textarrow',[0.445555555555556 0.47],...
    [0.658251769464105 0.682507583417594],'Color',[1 0 0],...
    'String',{'High-amplitude','erratic noise'},...
    'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textbox',...
    [0.649888888888889 0.601617795753287 0.133444444444444 0.0758341759352916],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textarrow',[0.836666666666667 0.854444444444444],...
    [0.828120323559153 0.864509605662287],'Color',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.817777777777778 0.836666666666667],...
    [0.910021233569268 0.866531850353895],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
subplot(2,1,2); imagesc(x,t,d_2l);colormap(amf_seis);hold on
text(n1/-70,-0.1,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/6.0,-0.07,'BP+SOSVMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.07,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,-0.07,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
annotation(gcf,'textbox',...
    [0.131 0.128412537917089 0.133444444444445 0.0758341759352915],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textarrow',[0.321111111111111 0.331111111111111],...
    [0.355915065722954 0.389282103134484],'Color',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.287777777777778 0.306666666666667],...
    [0.435804853387265 0.392315470171892],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.551111111111111 0.57],...
    [0.436815975733068 0.393326592517695],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.574444444444445 0.592222222222223],...
    [0.350870576339742 0.387259858442875],'Color',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.812222222222223 0.831111111111111],...
    [0.440860465116284 0.397371081900911],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.841111111111111 0.858888888888889],...
    [0.35289282103135 0.389282103134483],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');


print(gcf,'-depsc','-r300','fig2.eps');
%
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
imagesc(x,t,d_siml);colormap(jet);hold on
text(n2/9,-0.02,'BP','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.25,-0.02,'BP+SOSVMF','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'BP+SOSVMF+FK','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.40,-0.02,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.10,-0.02,'AMF','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 16;%c.Label.FontWeight = bold;
caxis([0,0.5]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
annotation(gcf,'textarrow',[0.471111111111111 0.617777777777778],...
    [0.640040444893833 0.763397371081902],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',30,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',25);
annotation(gcf,'textarrow',[0.47 0.505555555555555],...
    [0.640051567239637 0.724974721941355],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',30,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',25);
annotation(gcf,'textarrow',[0.47 0.372222222222222],...
    [0.641051567239636 0.715874620829121],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',30,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',25);

print(gcf,'-depsc','-r300','fig3.eps');
%% Plot figures curved events

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1); imagesc(x,t,d_1);colormap(amf_seis);hold on
text(n1/-70,-0.1,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/9.5,-0.05,'Clean','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/4,-0.05,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.05,'BP','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.17,-0.05,'BP+SOSVMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
annotation(gcf,'textarrow',[0.446666666666667 0.475555555555556],...
    [0.64106268958544 0.699696663296259],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'High-amplitude','erratic noise'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textbox',...
    [0.651 0.587462082912033 0.125666666666667 0.0768452982810964],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textarrow',[0.813333333333333 0.814444444444444],...
    [0.850353892821031 0.893832153690603],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.871111111111111 0.9],...
    [0.662296258847325 0.720930232558143],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'VIsible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

subplot(2,1,2); imagesc(x,t,d_2);colormap(amf_seis);hold on
text(n1/-70,-0.1,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/6.0,-0.05,'BP+SOSVMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.05,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,-0.05,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');

annotation(gcf,'textbox',...
    [0.134333333333333 0.124368048533874 0.125666666666667 0.0768452982810964],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textarrow',[0.323333333333333 0.324444444444444],...
    [0.353892821031345 0.404448938321543],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.348888888888889 0.378888888888889],...
    [0.194146612740146 0.248736097067746],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.558888888888889 0.573333333333333],...
    [0.347837209302332 0.409504550050557],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'VIsible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.821111111111111 0.832222222222222],...
    [0.351870576339737 0.405460060667348],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

print(gcf,'-depsc','-r300','fig4.eps');
%%
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
imagesc(x,t,d_sim);colormap(jet);hold on
text(n2/9,-0.02,'BP','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.25,-0.02,'BP+SOSVMF','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/2,-0.02,'BP+SOSVMF+FK','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.40,-0.02,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.10,-0.02,'AMF','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 16; % c.Label.FontWeight = bold;
caxis([0,0.5]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');

annotation(gcf,'textarrow',[0.551111111111114 0.553333333333336],...
    [0.304358948432761 0.378159757330637],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',25);
annotation(gcf,'arrow',[0.55 0.408888888888889],...
    [0.306370070778564 0.380182002022245],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.552222222222222 0.684444444444444],...
    [0.306381193124368 0.389282103134479],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15);

print(gcf,'-depsc','-r300','fig5');
%%






