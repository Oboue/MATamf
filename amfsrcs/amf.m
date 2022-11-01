function [dout] = amf(din,par)
%  AMF: Advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets
%
% By Oboue et al., 2022
%
% INPUT
% din:      input data
% par:      parameter file
% IN par:
% dt:       sampling
% flo:      Low frequency in band, default is 0
% fhi:      High frequency in band, default is Nyquist
% nplo:     number of poles for low cutoff
% nphi:     number of poles for high cutoff
% phase:    y: minimum phase, n: zero phase
% verb0:    verbosity flag
%
% niter:    number of nonlinear iterations
% liter:    number of linear iterations (in divn)
% order1:   accuracy order
% eps_dv:   eps for divn  (default: 0.01)
% eps_cg:   eps for CG    (default: 1)
% tol_cg:   tolerence for CG (default: 0.000001)
% rect:     smoothing radius (ndim*1)
% verb1:    verbosity flag
%
% adj:      adjoint flag
% add:      adding flag
% ns:       spray radius
% order2:   PWD order
% eps:      regularization (default:0.01);
% ndn:      size of dn (n1*n2)
% nds:      size of ds (n1*n2)
% type_mf:  0 (MF) or 1 (SVMF)
% ifsmooth: 1 (if smooth) or 0 (only MF)
%
% w:        half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
% c1:       Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
% c2:       Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
% c3:       Thresholding parameter (alpha)
% niter1:   Number of iteration
%
% rec:      3-D vector denoting smooth radius 
% eps1:     regularization parameter, default 0.0
% niter2:   number of CG iterations
% verb:     verbosity flag (default: 0) 
%
% OUTPUT
% dout:     output data
%%
if isfield(par,'dt')
    dt=par.dt;
else
    dt=0.004;
end

if isfield(par,'flo')
    flo=par.flo;
else
    flo=0;
end

if isfield(par,'fhi')
    fhi=par.fhi;
else
    fhi=1/dt/2;
end

nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1,rec,eps1,niter2,verb




[d0]=amf_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
% amf_bandpasssosvmffkcurvelet: robust denoising method composed of the bandpass filter, SOSVMF, FK, and the curvelet methods. 
% This function estimates the initially denoised signal d0

%% Local orthogonalization operation 
 nois_0=din-d0; % compute the initial noise section
[dout,nois2,low]=amf_localortho(d0,nois_0,rec,niter2,eps1,verb);

% amf_localortho: amf_localortho compensates for the signal leakage energy
% causes by ''amf_bandpasssosvmffkcurvelet'' method.
% dout is the final result from the AMF method.
return
