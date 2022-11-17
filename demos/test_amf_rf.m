% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

% Script to plot Figure 17

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

load d_z.mat 

[n1,n2]=size(d_z);
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
%
%% deconvolution
    itmax = 400;
    minderr = 0.001;
    ph = 5; % phase delay
    VB = 0; % Allow verbose output
    gauss=2.5;
    for n=1:size(d_amf_r,2)
        R=d_amf_r(:,n);
        Z=d_amf_z(:,n);
        [itr,itrms] = amf_makeRFitdecon_la_norm(R,Z,dt,nt,ph,gauss,itmax,minderr);
        sub(n).itr1=itr';
    end

    %% Plot figures

    figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');

    subplot(3,2,1)
    amf_wigb(d_z,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    text(-0.15,0.98,'(a)','Units','normalized','FontSize',16)
    title('Vertical raw')
    subplot(3,2,2)
    amf_wigb(d_amf_z,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    title('AMF')
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
    title('Radial raw')
    subplot(3,2,4)
    amf_wigb(d_amf_r,2,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    set(gca,'fontsize',16)
    title('AMF')
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
    title('RF raw')
   subplot(3,2,6)
    amf_wigb([sub.itr1],4,h,t)
    xlabel('Distance (km)')
    ylabel('Time (s)')
    ylim([-5 30])
    set(gca,'fontsize',16)
    title('AMF')
    text(-0.15,0.98,'(f)','Units','normalized','FontSize',18)
    annotation(gcf,'textarrow',[0.357777777777778 0.318888888888889],...
    [0.204465587044535 0.232793522267206],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'Moho'},...
    'LineWidth',3,...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.294444444444444 0.221111111111111],...
    [0.142724696356276 0.172064777327935],'Color',[1 0 0],'TextColor',[1 0 0],...
    'TextBackgroundColor',[1 1 1],...
    'String',{'High-amplitude','   erratic noise'},...
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


print(gcf,'-depsc','-r300','fig16]7.eps');
% end