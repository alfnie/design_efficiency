
function [E,params]=design_efficiency(varargin)
% DESIGN_EFFICIENCY Efficiency computation for fMRI block designs
%
% E=design_efficiency(argname1,argvalue1,argname2,argvalue2,...)
% possible argument names are:
%
%  mandatory arguments:
%   'conditions'    : Condition indexes for each block (vector of M numbers representing indexes ranging from 1 to N)
%   'blocklengths'  : Block lengths for each condition (vector of N numbers representing length in seconds of each condition block)
%   'TR'			: Repetition time (scalar representing acquisition period in seconds)
%	'contrast'	    : Contrast vector or matrix (vector of KxN numbers representing contrast weights)
%
%  optional arguments:
%	'hparam'		: High-pass filter value (in seconds). It defaults to 120
%	'nscans'		: Total number of scans. It defaults to last scan modeled by design
%   'blockdelays'   : Block onset delay for each condition (fixed block-onset delay) (vector of N numbers representing length in seconds of each condition fixed delay). It defaults to 0
%   'blockjitters'  : Block onset jitter for each condition (random block-onset delay) (vector of N numbers representing maximum length in seconds of each condition random delay). It defaults to 0
%   'model'         : 1: hrf model; 2: hrf+deriv model; 3: FIR model. It defaults to 1
%   'modelFIRorder' : FIR order. It defaults to 5
%   'modelFIRmaxT'  : FIR window size (in seconds). It defaults to 12s
%
% Example:
% E=design_efficiency('conditions',repmat([1,2],[1,10]),'blocklengths',[16,16],'contrast',[1,-1]);
%

% alfnie@gmail.com

fields={'conditions',[],...    
        'blocklengths',[],...
        'blockjitters',[],...
        'blockdelays',[],...
        'tr',[],...
        'model',1,...
        'modelFIRorder',5,...
        'modelFIRmaxT',12,...
        'contrast',[],...
		'hparam',120,...
		'nparam',120,...
		'nscans',[]};
params=[]; for n1=1:2:nargin, params=setfield(params,lower(varargin{n1}),varargin{n1+1}); end
for n1=1:2:length(fields), if ~isfield(params,fields{n1}) | isempty(getfield(params,fields{n1})), if isstr(fields{n1+1}), params=setfield(params,fields{n1},eval(fields{n1+1})); else, params=setfield(params,fields{n1},fields{n1+1}); end; end; end

nconditions=size(params.contrast,2);
params.conditions=double(params.conditions(:)'-min(params.conditions)+1);
params.blocklengths=params.blocklengths(:)';
if isempty(params.blockdelays), params.blockdelays=zeros(1,nconditions); end
if isempty(params.blockjitters), params.blockjitters=zeros(1,nconditions); end
params.blockdelays=params.blockdelays(:)';
params.blockjitters=params.blockjitters(:)';
blocklengths_percondition=true; %numel(params.blocklengths)==nconditions;

params.tr0=min(1,params.tr/16);
h=spm_hrf(params.tr0);
dh=diff(h);
mdh=sqrt(sum(dh.^2)/sum(h.^2));
maxT=ceil(params.modelFIRmaxT/params.tr0);
maxi=ceil(params.nscans*params.tr/params.tr0);

if blocklengths_percondition&&any(params.blockjitters)>0, nrepeat=10;
else nrepeat=1;
end

E=zeros(size(params.contrast,1),numel(params.model));
for inrepeat=1:nrepeat
    
% Creates design matrix
if blocklengths_percondition, times=[0,cumsum(params.blocklengths(params.conditions))]+[cumsum(params.blockdelays(params.conditions)+rand(size(params.conditions)).*params.blockjitters(params.conditions)),0];
else times=[0,cumsum(params.blocklengths)]; 
end
Ltotal=times(end);
%params.tr0=params.tr/ceil(params.tr*10);
if isempty(params.nscans), params.nscans=ceil(Ltotal/params.tr); 
elseif Ltotal<params.nscans*params.tr;
    while Ltotal<params.nscans*params.tr,
        params.conditions=[params.conditions,params.conditions];
        Ltotal=2*Ltotal;
    end
    if blocklengths_percondition, times=[0,cumsum(params.blocklengths(params.conditions))]+[cumsum(params.blockdelays(params.conditions)+rand(size(params.conditions)).*params.blockjitters(params.conditions)),0];
    else times=[0,cumsum(params.blocklengths)];
    end
    maskin=times<=params.nscans*params.tr;
    maskin(end)=false;
    params.conditions=params.conditions(maskin);
    times=times(maskin);
    Ltotal=times(end);
end
%t=(0:params.nscans-1)*params.tr0;

for nmodel=1:numel(params.model)
    model=params.model(nmodel);
    switch(model)
        case 1, nh=1; % hrf response
        case 2, nh=2; % hrf+hrf_derivative response
        case 3, nh=params.modelFIRorder; % FIR order
    end
    %basis=@(x,n)(1+cos(pi*(2*x-1))).*((2*x-1).^n)/2;
    basis=@(x,h)max(0,(nh-1)*min(x-(h-1)/(nh-1)+1,(h-1)/(nh-1)+1-x)-(nh-2));
    X=zeros(ceil(params.nscans*params.tr/params.tr0),nconditions*nh);
    for n1=1:length(params.conditions),
        if times(n1)/params.tr0+1>maxi, break; end
        switch(model)
            case 1
                idx=round(1+times(n1)/params.tr0):round(1+(times(n1)+params.blocklengths(params.conditions(n1)))/params.tr0);
                X(max(1,idx),params.conditions(n1))=1;
            case 2
                t0=round(times(n1)/params.tr0)+1;
                idx=t0:round(1+(times(n1)+params.blocklengths(params.conditions(n1)))/params.tr0);
                X(max(1,idx),1+nh*(params.conditions(n1)-1))=1;
                if t0>1, X(max(1,t0+(-1:0)),2+nh*(params.conditions(n1)-1))=params.tr/params.tr0*[-1;1]/mdh; end
            case 3,
                t0=round(times(n1)/params.tr0)+1;
                idx=t0+(0:maxT-1);
                pidx=find(idx>=1&idx<=size(X,1),maxT);
                a=(pidx-1)/(maxT-1);
                for inh=1:nh, X(idx(pidx),inh+nh*(params.conditions(n1)-1))=basis(a,inh); end
        end
    end
    
    % Convolves design matrix with hrf
    switch(model)
        case {1,2}, X=convn(X,h,'full');
        case 3,
    end
    X=X(1:maxi,:);
    params.X=X;
    
    % Convolves design matrix with high-pass filter & inverse noise
    n=2*size(X,1); %pow2(ceil(log2(size(X,1))));
    f=[0:n/2,n/2-1:-1:1]'/params.tr0/n;
    hpfilt=(f>1/params.hparam);
    noiseinv=hpfilt./(1+1./max(eps,f*params.nparam));
    Xf=fft([X;flipud(X)],n);
    Xf=repmat(noiseinv,[1,size(Xf,2)]).*Xf; % whitens design
    X=real(ifft(Xf,n));
    X=X(1:n/2,:);
    k=1/mean(1./noiseinv(noiseinv>0).^2); % 1/mse
    params.X2=X;
    params.t=(0:size(X,1)-1)*params.tr0;
    
    % Resamples design matrix and computes efficiency
    n0=round(params.tr/params.tr0);
    x=0; for n1=1:n0 x=x+X(min(size(X,1),n1+(0:params.nscans-1)*n0),:); end; X=x/n0;
    iXX=pinv(X'*X);
    switch(model)
        case 1,     C=permute(params.contrast,[3,2,1]);
        case {2,3}, C=permute(reshape(kron(params.contrast,eye(nh)),nh,size(params.contrast,1),size(params.contrast,2)*nh),[1,3,2]);
    end
    for n1=1:size(params.contrast,1),
        r=trace(C(:,:,n1)*iXX*C(:,:,n1)')/size(C,1);
        E(n1,nmodel)=E(n1,nmodel)+max(0,k./r)/nrepeat;
        %E(n1,nmodel)=sqrt(max(0,k./r));
    end
end
end

end




