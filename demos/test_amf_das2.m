% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets

% Script to plot Figure 8 

%  Copyright (C) Oboue et al., 2022

%  This program is free software: you canz redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.

%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
%  GNU General Public License for more details: http://www.gnu.org/licenses/

%  References

%  Oboue et al., 2022
%  Huang, G., M. Bai, Q. Zhao, W. Chen, and Y. Chen, 2021, Erratic noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, 69, 101-121.
%  Chen, Y., S. Zu, Y. Wang, and X. Chen, 2020, Deblending of simultaneous-source data using a structure-oriented space varying median filter, Geophysical Journal International, 222, 1805�1�723.
%  Zhao, Q., Q. Du, X. Gong, and Y. Chen, 2018, Signal-preserving erratic noise attenuation via iterative robust sparsity-promoting filter, IEEE Transactions on Geoscience and Remote Sensing, 56, 1558-0644.
%
clc;clear;close all;
%% addpath
addpath('../amf_data/');
addpath('../amfsrcs/');
addpath('../seistr/');
%% Load the DAS data

NOs=[1,3,20,10,25,11,2];
labels={...                                           %P-arrival sample NO from the SEGY file
    'FORGE\_78-32\_iDASv3-P11\_UTC190423150554.sgy',... %24169
    'FORGE\_78-32\_iDASv3-P11\_UTC190423213209.sgy'};   %4210
eq=zeros(2000,960);
[n1,n2]=size(eq);
t=[0:n1]*0.0005;
ngap=50;
x=1:n2*5+4*ngap;

%% Process DAS seismic data with weak signal energy corrupted by strong noise using the AMF method
%%
for ii=1
    
if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('amf_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
din=d1;
%% Parameter tuning 
%
par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=200;   % High frequency in band, default is Nyquist
par.nplo=6;    % number of poles for low cutoff
par.nphi=6;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag
%
par.niter=2;                      % number of nonlinear iterations
par.liter=10;                     % liter: number of linear iterations (in divn)
par.order1=3;                     % order: accuracy order
par.eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
par.tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=50;                   % rect:  smoothing radius (ndim*1)
par.rect(2)=50;                   % "      "        "
par.rect(3)=1;                    % "      "        "
par.verb1=1;                      % verbosity flag

par.adj=0;                        % adjoint flag
par.add=0;                        % adding flag
par.ns=8;                         % spray radius
par.order2=2;                     % PWD order
par.eps=0.01;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % size of ds (n1*n2)
par.type_mf=1;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
par.w=0.05;              % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=0.5;               % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration
%
rec = zeros(3, 1);   % 3-D vector denoting smooth radius 
par.rec(1) = 1000;
par.rec(2) = 1000;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 

tic
d_amf=amf(din,par);
d_amf1=d_amf;
toc
%
end

comp1=[din,d_amf1,din-d_amf1]; 
%% Process the DAS seismic data with strong signal energy corrupted by strong noise using the AMF method
%%

eq=zeros(2000,960);
[n1,n2]=size(eq);

for ii=3
    if ~ismember(ii,[14,16,17,27,47,52])
        strcat('amf_data/eq-',num2str(ii),'.mat')
        load(strcat('amf_data/eq-',num2str(ii),'.mat'));
    end
    d1=d1;
    eq=d1;
din=d1;
%% Parameter tuning 
%
par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=200;   % High frequency in band, default is Nyquist
par.nplo=6;    % number of poles for low cutoff
par.nphi=6;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag
%
par.niter=2;                      % number of nonlinear iterations
par.liter=10;                     % liter: number of linear iterations (in divn)
par.order1=3;                     % order: accuracy order
par.eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
par.tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=100;                   % rect:  smoothing radius (ndim*1)
par.rect(2)=100;                   % "      "        "
par.rect(3)=1;                    % "      "        "
par.verb1=1;                      % verbosity flag

par.adj=0;                        % adjoint flag
par.add=0;                        % adding flag
par.ns=8;                        % spray radius
par.order2=2;                     % PWD order
par.eps=0.01;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % size of ds (n1*n2)
par.type_mf=1;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
par.w=0.05;              % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=5;                 % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration
%
rec = zeros(3, 1);       % 3-D vector denoting smooth radius 
par.rec(1) = 30;
par.rec(2) = 30;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 
%
tic
d_amf=amf(din,par);
d_amf2=d_amf;
toc
%
end
%% Plot figures 

comp2=[din,d_amf2,din-d_amf2]; 

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1);amf_imagesc(comp1,95,1,x,t);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Panel','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
set(gca,'XTick',0)
text(n2/1.2,-0.03,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.38,-0.03,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.23,-0.03,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{1},'color','b','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'textarrow',[0.511111111111111 0.486666666666667],...
    [0.827456104944501 0.870837537840565],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible signal'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.523333333333333 0.495555555555555],...
    [0.732602421796166 0.688193743693239],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
subplot(2,1,2);amf_imagesc(comp2,95,1,x,t);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Panel','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
set(gca,'XTick',0)

text(n2/1.2,-0.03,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.38,-0.03,'AMF','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.23,-0.03,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(b)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{2},'color','b','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'textarrow',[0.318888888888889 0.358888888888889],...
    [0.386487386478305 0.420787083753784],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.217777777777778 0.217777777777778],...
    [0.25732492431887 0.214934409687185],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.487777777777778 0.487777777777778],...
    [0.259343087790111 0.216952573158426],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.604444444444444 0.622222222222222],...
    [0.37840565085772 0.420787083753784],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Improved','signal'},...
    'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

annotation(gcf,'textbox',...
    [0.365444444444444 0.900883218841989 0.0312222222222223 0.030421982335618],...
    'String','1',...
    'FontWeight','bold',...
    'FontSize',25,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textbox',...
    [0.621 0.899901864573108 0.0256666666666667 0.0314033366045119],...
    'String',{'2'},...
    'FontWeight','bold',...
    'FontSize',25,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textbox',...
    [0.876555555555556 0.897939156035327 0.029 0.0343473994111859],...
    'String',{'3'},...
    'FontWeight','bold',...
    'FontSize',25,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'line',[0.387777777777778 0.387777777777778],...
    [0.925417075564279 0.581943081452404],'LineWidth',3);
annotation(gcf,'line',[0.646666666666667 0.646666666666667],...
    [0.924435721295387 0.584887144259078],'LineWidth',3);

annotation(gcf,'textbox',...
    [0.361 0.426889106967615 0.0334444444444445 0.0304219823356219],...
    'String','1',...
    'FontWeight','bold',...
    'FontSize',25,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'line',[0.388888888888889 0.387777777777778],...
    [0.449441609421001 0.107948969578018],'LineWidth',3);
annotation(gcf,'line',[0.648888888888889 0.648888888888889],...
    [0.450441609421001 0.108930323846909],'LineWidth',3);

annotation(gcf,'textbox',...
    [0.623222222222222 0.428851815505396 0.027888888888889 0.0294406280667329],...
    'String',{'2'},...
    'FontWeight','bold',...
    'FontSize',25,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textbox',...
    [0.876555555555556 0.433758586849851 0.027888888888889 0.0245338567222756],...
    'String',{'3'},...
    'FontWeight','bold',...
    'FontSize',25,...
    'FitBoxToText','off',...
    'EdgeColor','none');
print(gcf,'-depsc','-r300','fig8.eps');