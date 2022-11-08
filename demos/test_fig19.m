% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

%  Copyright (C) Oboue et al., 2022

%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

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
%% %Generate the clean synthetic data including linear events

dc=amf_levents(200);
dc=amf_scale(dc);
[n1,n2]=size(dc);
%% %Generate the noisy synthetic data having linear events

% Create the mask operator
mask=rand(1,n2);
mask(logical(mask<0.9))=0;
mask(logical(mask>=0.9))=1;

% Adding erratic noise 

err_n=zeros(size(dc));
for i=1:n1
    randn('state',123456+i);
    err_n(i,:)=0.5*randn(1,n2).*mask;
end

% Adding random noise 

randn('state',201920);
ran_n=0.1*randn(n1,n2);

dn=dc+err_n+ran_n;
din=dn;
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

%                                                                                                             0,0,dipn,[],n1,n2,ns,2,0.01,n1*n2,n1*n2,type_mf,ifsmooth,d1,[]
tic
d2=amf_bandpasssosvmf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
% 
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method

w=0.0001; % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

% 
tic
d3 = amf_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
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
d4=amf_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
%% Denoising using the AMF method 
% % Parameter tuning: add the key parameters for the local orthogonalization operation




%% AMF method
par.dt=dt; % sampling
par.flo=flo;     % Low frequency in band, default is 0
par.fhi=fhi;   % High frequency in band, default is Nyquist
par.nplo=nplo;    % number of poles for low cutoff
par.nphi=nphi;    % number of poles for high cutoff
par.phase=phase;   % y: minimum phase, n: zero phase
par.verb0=verb0;   % verbosity flag
%
par.niter=niter;                      % number of nonlinear iterations
par.liter=liter;                     % liter: number of linear iterations (in divn)
par.order1=order1;                     % order: accuracy order
par.eps_dv=eps_dv;                  % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=eps_cg;                     % eps_cg: eps for CG    (default: 1)
par.tol_cg=tol_cg;              % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=11;                   % rect:  smoothing radius (ndim*1)
par.rect(2)=5;                    % "      "        "
par.rect(3)=1;                    % "      "        "
par.verb1=verb1;                      % verbosity flag

par.adj=adj;                        % adjoint flag
par.add=add;                        % adding flag
par.ns=ns;                         % spray radius
par.order2=order2;                     % PWD order
par.eps=eps;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % size of ds (n1*n2)
par.type_mf=type_mf;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=ifsmooth;                   % 1 (if smooth) or 0 (only MF)
% 
par.w=w;               % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=c1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=c3;              % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration

rec = zeros(3, 1);       % 3-D vector denoting smooth radius 
par.rec(1) = 5;
par.rec(2) = 5;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 
%
tic
d_amf=amf(din,par);
d5=d_amf;
toc
%
%% Plot figures

t=[0:n1]*0.004;
x=1:n2;
ngap=5;
d_1=[d3,dn-d3,zeros(n1,ngap),d4,d4-dn,zeros(n1,ngap),d5,dn-d5];

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1); imagesc(x,t,d_1);colormap(amf_seis);hold on
text(n2/200,-0.2,'(a)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/6.0,-0.05,'BP+SOSVMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.05,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,-0.05,'AMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
colormap(amf_seis);caxis([-0.6,0.6]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');

%%
% local similarity maps
rect=[10,10,1];niter=20;eps=0;verb=0;
[simid3]=amf_localsimi(din-d3,d3,rect,niter,eps,verb);
[simid4]=amf_localsimi(din-d4,d4,rect,niter,eps,verb);
[simid5]=amf_localsimi(din-d5,d5,rect,niter,eps,verb);

dsim=[simid3,zeros(n1,ngap),simid4,zeros(n1,ngap),simid5];

subplot(2,1,2); imagesc(x,t,dsim);colormap(jet);hold on
text(n2/200,-0.2,'(b)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/6.0,-0.05,'BP+SOSVMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,-0.05,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.19,-0.05,'AMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 16;%c.Label.FontWeight = bold;
caxis([0,0.5]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Trace','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');

print(gcf,'-depsc','-r300','fig19-amf.eps');
%% Numerical test using signal-to-noise ratio (S/N)

amf_snr(dc,dn) %Noisy
amf_snr(dc,d1) %BP
amf_snr(dc,d2) %BP+SOSVMF
amf_snr(dc,d3) %BP+SOSVMF+FK
amf_snr(dc,d4) %BP+SOSVMF+FK+Curvelet
amf_snr(dc,d5) %AMF
