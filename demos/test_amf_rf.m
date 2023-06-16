% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

% Script to plot Figures 9 and 10

%  Copyright (C) Oboue et al., 2022

% Oboue Nov. 2022
%%
% Oboue Nov. 2022 Global Seismology Group, Zhejiang University
% Denoise receiver functions using AMF method
%
% Oct. 19, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University
% Denoise receiver functions using Radon transform
%
% clear; clc; close all;
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
addpath('../amfsrcs/');
addpath('../seistr/');
addpath('../amf_data/');
%% load the vertical component of the receiver function 

load d_z.mat % vertical component 

load d_r.mat % radial component

[n1,n2]=size(d_z);
%%  Denosing using the BP method
% Parameter tuning for the BP method

% dt=t; % sampling
flo=0;  % Low frequency in band, default is 0
fhi=200;  % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
tic
d_bp_z=amf_bp(d_z,t,flo,fhi,nplo,nphi,phase,verb0);
toc

tic
d_bp_r=amf_bp(d_r,t,flo,fhi,nplo,nphi,phase,verb0);
toc
%% Denosing using the BP+SOSVMF method 
% Parameter tuning：add the key parameters of the SOSVMF method

niter=5;                      % number of nonlinear iterations
liter=20;                     % liter: number of linear iterations (in divn)
order1=2;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=200;                   % rect:  smoothing radius (ndim*1)
rect(2)=200;                   % "      "        "
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
d_bpsosvmf_z=amf_bpsosvmf(d_z,t,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%

niter=5;                      % number of nonlinear iterations
liter=20;                     % liter: number of linear iterations (in divn)
order1=2;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=1000;                   % rect:  smoothing radius (ndim*1)
rect(2)=1000;                   % "      "        "
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

tic
d_bpsosvmf_r=amf_bpsosvmf(d_r,t,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
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
rect(1)=200;                   % rect:  smoothing radius (ndim*1)
rect(2)=200;                   % "      "        "
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

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=500;               % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

w=0.0;
% 
tic
d_bpsosvmffkct_z=amf_bpsosvmffkct(d_z,t,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc

niter=5;                      % number of nonlinear iterations
liter=20;                     % liter: number of linear iterations (in divn)
order1=2;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=1000;                   % rect:  smoothing radius (ndim*1)
rect(2)=1000;                   % "      "        "
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

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=500;               % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
tic
d_bpsosvmffkct_r=amf_bpsosvmffkct(d_r,t,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%% AMF
% dt=0.0005; % sampling
par.dt=t; % sampling
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
par.rect(1)=200;                  % rect:  smoothing radius (ndim*1)
par.rect(2)=200;                  % "      "        "
par.rect(3)=1;                    % "      "        "
par.verb1=1;                      % verbosity flag

par.adj=0;                        % adjoint flag
par.add=0;                        % adding flag
par.ns=2;                         % spray radius
par.order2=2;                     % PWD order
par.eps=0.01;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % size of ds (n1*n2)
par.type_mf=1;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)

par.w=0;

par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=500;               % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration

rec = zeros(3, 1);       % 3-D vector denoting smooth radius 
par.rec(1) = 200;
par.rec(2) = 200;
par.rec(3) = 1;
par.eps1=0;              % regularization parameter, default 0.0
par.niter2=20;           % number of CG iterations
par.verb=1;              % verbosity flag (default: 0) 
%
tic
d_amf_z=amf(d_z,par);
toc
%

%% load the radial component of the receiver function
 
load d_r.mat
[n1,n2]=size(d_r);
%% AMF
% dt=0.0005; % sampling
par.dt=t; % sampling
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
par.rect(1)=1000;                  % rect:  smoothing radius (ndim*1)
par.rect(2)=1000;                  % "      "        "
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

par.w=0;

par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=2;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=500;               % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration

rec = zeros(3, 1);       % 3-D vector denoting smooth radius 
par.rec(1) = 500;
par.rec(2) = 500;
par.rec(3) = 1;
par.eps1=0;              % regularization parameter, default 0.0
par.niter2=20;           % number of CG iterations
par.verb=1;              % verbosity flag (default: 0) 
%
tic
d_amf_r=amf(d_r,par);
toc
%%
% deconvolution
% raw data

    itmax = 400;
    minderr = 0.001;
    ph = 5; % phase delay
    VB = 0; % Allow verbose output
    gauss=2.5;

    for n=1:size(d_r,2)
        R=d_r(:,n);
        Z=d_bp_z(:,n);
        [itr,itrms] = amf_makeRFitdecon_la_norm(R,Z,dt,nt,ph,gauss,itmax,minderr);
        sub(n).itr=itr';
    end
    
    
 % BP 

    for n=1:size(d_bp_r,2)
        R=d_bp_r(:,n);
        Z=d_bp_z(:,n);
        [itr,itrms] = amf_makeRFitdecon_la_norm(R,Z,dt,nt,ph,gauss,itmax,minderr);
        sub(n).itrbp=itr';
    end
    
% BP+SOSVMF

    for n=1:size(d_bpsosvmf_r,2)
        R=d_bpsosvmf_r(:,n);
        Z=d_bpsosvmf_z(:,n);
        [itr,itrms] = amf_makeRFitdecon_la_norm(R,Z,dt,nt,ph,gauss,itmax,minderr);
        sub(n).itrbpsosvmf=itr';
    end

% BP+SOSVMF+Curvelet

    for n=1:size(d_bpsosvmffkct_r,2)
        R=d_bpsosvmffkct_r(:,n);
        Z=d_bpsosvmffkct_z(:,n);
        [itr,itrms] = amf_makeRFitdecon_la_norm(R,Z,dt,nt,ph,gauss,itmax,minderr);
        sub(n).itrbpsosvmfct=itr';
    end
    
% AMF

    for n=1:size(d_amf_r,2)
        R=d_amf_r(:,n);
        Z=d_amf_z(:,n);
        [itr,itrms] = amf_makeRFitdecon_la_norm(R,Z,dt,nt,ph,gauss,itmax,minderr);
        sub(n).itramf=itr';
    end
    %% Plot figures
    
    % BP

    figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');

    subplot(3,2,1)
    amf_wigb(d_z,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    text(-0.15,0.98,'(a)','Units','normalized','FontSize',16)
    title('Vertical (raw)')
    subplot(3,2,2)
    amf_wigb(d_bp_z,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    title('Vertical (BP)')
    text(-0.15,0.98,'(b)','Units','normalized','FontSize',18)
    annotation(gcf,'arrow',[0.228888888888889 0.265555555555556],...
    [0.921052631578948 0.903846153846154],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    annotation(gcf,'arrow',[0.431111111111111 0.385555555555556],...
    [0.918016194331984 0.897773279352227],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    % 
    subplot(3,2,3)
    amf_wigb(d_r,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    text(-0.15,0.98,'(c)','Units','normalized','FontSize',18)
    title('Radial (raw)')
    subplot(3,2,4)
    amf_wigb(d_bp_r,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    title('Radial (BP)')
    text(-0.15,0.98,'(d)','Units','normalized','FontSize',18)     
    annotation(gcf,'arrow',[0.677777777777778 0.714444444444445],...
    [0.919028340080972 0.901821862348179],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    annotation(gcf,'arrow',[0.886666666666667 0.841111111111111],...
    [0.915991902834008 0.895748987854251],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.18 0.216666666666667],...
    [0.607287449392713 0.590080971659919],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.232222222222222 0.268888888888889],...
    [0.615384615384615 0.598178137651822],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.232222222222222 0.268888888888889],...
    [0.615384615384615 0.598178137651822],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.425555555555556 0.392222222222222],...
    [0.612348178137652 0.597165991902834],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.62 0.656666666666667],...
    [0.61336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.62 0.656666666666667],...
    [0.61336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.694444444444444 0.731111111111111],...
    [0.614372469635628 0.597165991902834],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.884444444444444 0.851111111111111],...
    [0.611336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);

    subplot(3,2,5)
    amf_wigb([sub.itr],4,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    ylim([-5 30])
    set(gca,'fontsize',16)
    text(-0.15,0.98,'(e)','Units','normalized','FontSize',18)
    title('RF (raw)')
   subplot(3,2,6)
    amf_wigb([sub.itrbp],4,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    ylim([-5 30])
    set(gca,'fontsize',16)
    title('RF (BP)')
    text(-0.15,0.98,'(f)','Units','normalized','FontSize',18)
    annotation(gcf,'textarrow',[0.357777777777778 0.318888888888889],...
    [0.204465587044535 0.232793522267206],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Moho'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.29 0.274444444444444],...
    [0.163029264884195 0.185786802030457],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Strong','   erratic noise'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');

annotation(gcf,'textarrow',[0.795555555555555 0.764444444444444],...
    [0.203436983641093 0.228658536585366],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Moho'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');   
  %% 
    % BP+SOSVMF
    
 figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');

    subplot(3,2,1)
    amf_wigb(d_z,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    text(-0.15,0.98,'(a)','Units','normalized','FontSize',16)
    title('Vertical (raw)')
    subplot(3,2,2)
    amf_wigb(d_bpsosvmf_z,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    title('Vertical (BP+SOSVMF)')
    text(-0.15,0.98,'(b)','Units','normalized','FontSize',18)
    annotation(gcf,'arrow',[0.228888888888889 0.265555555555556],...
    [0.921052631578948 0.903846153846154],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    annotation(gcf,'arrow',[0.431111111111111 0.385555555555556],...
    [0.918016194331984 0.897773279352227],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    % 
    subplot(3,2,3)
    amf_wigb(d_r,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    text(-0.15,0.98,'(c)','Units','normalized','FontSize',18)
    title('Radial (raw)')
    subplot(3,2,4)
    amf_wigb(d_bpsosvmf_r,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    title('Radial (BP+SOSVMF)')
    text(-0.15,0.98,'(d)','Units','normalized','FontSize',18)     
    annotation(gcf,'arrow',[0.677777777777778 0.714444444444445],...
    [0.919028340080972 0.901821862348179],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    annotation(gcf,'arrow',[0.886666666666667 0.841111111111111],...
    [0.915991902834008 0.895748987854251],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.18 0.216666666666667],...
    [0.607287449392713 0.590080971659919],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.232222222222222 0.268888888888889],...
    [0.615384615384615 0.598178137651822],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.232222222222222 0.268888888888889],...
    [0.615384615384615 0.598178137651822],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.425555555555556 0.392222222222222],...
    [0.612348178137652 0.597165991902834],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.62 0.656666666666667],...
    [0.61336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.62 0.656666666666667],...
    [0.61336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.694444444444444 0.731111111111111],...
    [0.614372469635628 0.597165991902834],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.884444444444444 0.851111111111111],...
    [0.611336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);

    subplot(3,2,5)
    amf_wigb([sub.itr],4,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    ylim([-5 30])
    set(gca,'fontsize',16)
    text(-0.15,0.98,'(e)','Units','normalized','FontSize',18)
    title('RF (raw)')
   subplot(3,2,6)
    amf_wigb([sub.itrbpsosvmf],4,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    ylim([-5 30])
    set(gca,'fontsize',16)
    title('RF (BP+SOSVMF)')
    text(-0.15,0.98,'(f)','Units','normalized','FontSize',18)
    annotation(gcf,'textarrow',[0.357777777777778 0.318888888888889],...
    [0.204465587044535 0.232793522267206],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Moho'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.29 0.274444444444444],...
    [0.163029264884195 0.185786802030457],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Strong','   erratic noise'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');

annotation(gcf,'textarrow',[0.795555555555555 0.764444444444444],...
    [0.203436983641093 0.228658536585366],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Moho'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');   
    
  %%  
    % BP+SOSVMF+Curvelet
    
    figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');

    subplot(3,2,1)
    amf_wigb(d_z,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    text(-0.15,0.98,'(a)','Units','normalized','FontSize',20)
    title('Vertical (raw)')
    subplot(3,2,2)
    amf_wigb(d_bpsosvmffkct_z,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    title('Vertical (BP+SOSVMF+Curvelet)')
    text(-0.15,0.98,'(b)','Units','normalized','FontSize',20)
    annotation(gcf,'arrow',[0.228888888888889 0.265555555555556],...
    [0.921052631578948 0.903846153846154],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    annotation(gcf,'arrow',[0.431111111111111 0.385555555555556],...
    [0.918016194331984 0.897773279352227],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    % 
    subplot(3,2,3)
    amf_wigb(d_r,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    text(-0.15,0.98,'(c)','Units','normalized','FontSize',20)
    title('Radial (raw)')
    subplot(3,2,4)
    amf_wigb(d_bpsosvmffkct_r,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    title('Radial (BP+SOSVMF+Curvelet)')
    text(-0.15,0.98,'(d)','Units','normalized','FontSize',18)     
    annotation(gcf,'arrow',[0.677777777777778 0.714444444444445],...
    [0.919028340080972 0.901821862348179],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    annotation(gcf,'arrow',[0.886666666666667 0.841111111111111],...
    [0.915991902834008 0.895748987854251],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.18 0.216666666666667],...
    [0.607287449392713 0.590080971659919],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.232222222222222 0.268888888888889],...
    [0.615384615384615 0.598178137651822],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.232222222222222 0.268888888888889],...
    [0.615384615384615 0.598178137651822],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.425555555555556 0.392222222222222],...
    [0.612348178137652 0.597165991902834],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.62 0.656666666666667],...
    [0.61336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.62 0.656666666666667],...
    [0.61336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.694444444444444 0.731111111111111],...
    [0.614372469635628 0.597165991902834],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.884444444444444 0.851111111111111],...
    [0.611336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);

    subplot(3,2,5)
    amf_wigb([sub.itr],4,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    ylim([-5 30])
    set(gca,'fontsize',16)
    text(-0.15,0.98,'(e)','Units','normalized','FontSize',20)
    title('RF (raw)')
   subplot(3,2,6)
    amf_wigb([sub.itrbpsosvmfct],4,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    ylim([-5 30])
    set(gca,'fontsize',16)
    title('RF (BP+SOSVMF+Curvelet)')
    text(-0.15,0.98,'(f)','Units','normalized','FontSize',20)
    annotation(gcf,'textarrow',[0.357777777777778 0.318888888888889],...
    [0.204465587044535 0.232793522267206],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Moho'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.29 0.274444444444444],...
    [0.163029264884195 0.185786802030457],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Strong','   erratic noise'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.795555555555555 0.764444444444444],...
    [0.203436983641093 0.228658536585366],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Moho'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');   
  %%  
% AMF
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');

    subplot(3,2,1)
    amf_wigb(d_z,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    text(-0.2,1.05,'(a)','Units','normalized','FontSize',20)
    title('Vertical (Raw)')
    subplot(3,2,2)
    amf_wigb(d_amf_z,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    title('Vertical (AMF)')
    text(-0.2,1.05,'(b)','Units','normalized','FontSize',20)
    annotation(gcf,'arrow',[0.228888888888889 0.265555555555556],...
    [0.921052631578948 0.903846153846154],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    annotation(gcf,'arrow',[0.431111111111111 0.385555555555556],...
    [0.918016194331984 0.897773279352227],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    % 
    subplot(3,2,3)
    amf_wigb(d_r,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    text(-0.2,1.05,'(c)','Units','normalized','FontSize',20)
    title('Radial raw')
    subplot(3,2,4)
    amf_wigb(d_amf_r,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    title('Radial (AMF)')
    text(-0.2,1.05,'(d)','Units','normalized','FontSize',20)     
    annotation(gcf,'arrow',[0.677777777777778 0.714444444444445],...
    [0.919028340080972 0.901821862348179],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
    annotation(gcf,'arrow',[0.886666666666667 0.841111111111111],...
    [0.915991902834008 0.895748987854251],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.18 0.216666666666667],...
    [0.607287449392713 0.590080971659919],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.232222222222222 0.268888888888889],...
    [0.615384615384615 0.598178137651822],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.232222222222222 0.268888888888889],...
    [0.615384615384615 0.598178137651822],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.425555555555556 0.392222222222222],...
    [0.612348178137652 0.597165991902834],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.62 0.656666666666667],...
    [0.61336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.62 0.656666666666667],...
    [0.61336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.694444444444444 0.731111111111111],...
    [0.614372469635628 0.597165991902834],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.884444444444444 0.851111111111111],...
    [0.611336032388664 0.596153846153846],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',15,...
    'HeadLength',15);

    subplot(3,2,5)
    amf_wigb([sub.itr],4,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    ylim([-5 30])
    set(gca,'fontsize',16)
    text(-0.2,1.05,'(e)','Units','normalized','FontSize',20)
    title('RF (Raw)')
   subplot(3,2,6)
    amf_wigb([sub.itramf],4,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    ylim([-5 30])
    set(gca,'fontsize',16)
    title('RF (AMF)')
    text(-0.2,1.05,'(f)','Units','normalized','FontSize',20)
    annotation(gcf,'textarrow',[0.357777777777778 0.318888888888889],...
    [0.204465587044535 0.232793522267206],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Moho'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.29 0.274444444444444],...
    [0.163029264884195 0.185786802030457],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Strong','   erratic noise'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.771111111111111 0.732222222222222],...
    [0.207502024291499 0.23582995951417],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Moho'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');

print(gcf,'-depsc','-r300','fig9.eps');
%% plot receiver functions data together
    
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');

    subplot(3,2,1)
    imagesc(dist0,t,[sub.itr]);
    colormap(amf_seis);
    caxis([-0.4 0.4])
    xlabel('Distance (km)');
    ylabel('Time (sec)');
    set(gca,'fontsize',16)
    title('RF (Raw)')
    ylim([-5 30])
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
    ylim([-5 30])
    text(-0.2,1.05,'(a)','Units','normalized','FontSize',25);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Distance (km)','Fontsize',16,'fontweight','bold');
annotation(gcf,'textarrow',[0.241111111111111 0.221111111111111],...
    [0.788832487309646 0.808418792001484],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Strong','erratic','noise'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.313333333333333 0.295555555555555],...
    [0.806091370558376 0.83379950266138],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Moho'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textbox',...
    [0.368777777777778 0.719796954314721 0.0823333333333334 0.0497461928934038],...
    'Color',[1 0 0],...
    'String',{'Random','   noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');

 subplot(3,2,2)
    imagesc(dist0,t,[sub.itrbp]);
    colormap(amf_seis);
    caxis([-0.4 0.4])
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
    title('RF (BP)')
    ylim([-5 30])
    text(-0.2,1.05,'(b)','Units','normalized','FontSize',25)
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Distance (km)','Fontsize',16,'fontweight','bold');
annotation(gcf,'textarrow',[0.678888888888889 0.661111111111111],...
    [0.768527918781728 0.800296964590317],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.757777777777778 0.737777777777778],...
    [0.80507614213198 0.833799502661382],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');

 subplot(3,2,3)
    imagesc(dist0,t,[sub.itrbpsosvmf]);
    colormap(amf_seis);
    caxis([-0.4 0.4])
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
    title('RF (BP+SOSVMF)')
    ylim([-5 30])
    text(-0.2,1.05,'(c)','Units','normalized','FontSize',25)
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Distance (km)','Fontsize',16,'fontweight','bold');
annotation(gcf,'textarrow',[0.235555555555556 0.217777777777778],...
    [0.479187817258884 0.510956863067474],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.304444444444444 0.286666666666666],...
    [0.50253807106599 0.53430711687458],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');

    subplot(3,2,4)
    imagesc(dist0,t,[sub.itrbpsosvmfct]);
    colormap(amf_seis);
    caxis([-0.4 0.4])
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
    title('RF (BP+SOSVMF+Curvelet)')
    ylim([-5 30])
    text(-0.2,1.05,'(d)','Units','normalized','FontSize',25)
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Distance (km)','Fontsize',16,'fontweight','bold');
   annotation(gcf,'textarrow',[0.737777777777778 0.723333333333333],...
    [0.505583756345178 0.535322345300976],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.882222222222222 0.867777777777778],...
    [0.507614213197971 0.537352802153768],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
    
 subplot(3,2,5)
    imagesc(dist0,t,[sub.itramf]);
    colormap(amf_seis);
    caxis([-0.4 0.4])
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
    title('RF (AMF)')
    ylim([-5 30])
    text(-0.2,1.05,'(e)','Units','normalized','FontSize',25)
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Distance (km)','Fontsize',16,'fontweight','bold');
annotation(gcf,'textarrow',[0.29 0.272222222222222],...
    [0.201015228426396 0.236845187940569],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Improved','signal','energy  of Moho'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.437777777777778 0.42],...
    [0.2 0.235829959514173],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');

print(gcf,'-depsc','-r300','fig10.eps');