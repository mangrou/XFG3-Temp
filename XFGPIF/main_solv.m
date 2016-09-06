function main_solv

load('test')
t=0;
mixsig=msig;W=B;p=pa;lambda=0.04;ll=1;
p=solvepa(mixsig,W,p);
[W,p]=v(p,W);
err=1;
while(err>0.01)
    Wt=solveB(mixsig,W,p,lambda,ll);
    pt=solvepa(mixsig,Wt,p);
    [Wt,pt]=v(pt,Wt);
    err=max(max(max(abs(W-Wt))),max(max(abs(p-pt))));
    W=Wt;
    p=pt;
    t=t+1;
    if mod(t, 10)==0
        ll=ll*0.5;
    end
end
for i=1:3
    for j=1:3
        if abs(W(i,j))<0.0001
            W(i,j)=0;
        end
    end
end

PW=W;
pa=p;
sig=PW*msig;
save('result')
ff=figure
subplot(3,1,1)
plot(sig(1,:))
subplot(3,1,2)
plot(sig(2,:))
subplot(3,1,3)
plot(sig(3,:))
print('Sig','-dpng')

close(ff)

ff=imagesc(zeros(3));
colormap(flipud(gray));
textStrings = num2str(W(:),'%0.2f');
textStrings = strtrim(cellstr(textStrings));
[x,y] = meshgrid(1:3);
hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
set(gca,'visible','off')
print('Loading Matrix','-dpng')
close

end


function [Wt,pt]=v(pa,wa)
% scale the loading matrix and NIG parameters such that the IFs are of unit
% variance. 
[n,~]=size(wa);
x=pa(:,4).*pa(:,1).^2./sqrt(pa(:,1).^2-pa(:,2).^2).^3;
x=sqrt(x);
pa(:,1)=pa(:,1).*x;
pa(:,2)=pa(:,2).*x;
pa(:,3)=pa(:,3)./x;
pa(:,4)=pa(:,4)./x;
pt=pa;
xx=zeros(n);
for i=1:n
    xx(:,i)=x;
end
Wt=wa./xx;
end

function pa=solvepa(mixsig,W,PA)
[n,~]=size(mixsig);
p=zeros(n,4);
x=(W*mixsig)';
for i=1:n
    [p(i,1),p(i,2),p(i,3),p(i,4)]=nigem(x(:,i)',0.01,100,PA(i,:));
end
pa=p;
end

function y=solveB(mixsig,W,pa,lambda,l)
t=0;
alpha=l;
error=1;
while(error>0.0001 && t<200) 
    d=dldbM(mixsig,W,lambda,pa);
    dw=W+alpha*d;
    % [d; dw]  uncommend to watch the change of matrix and derevative
    t=t+1;
    if mod(t,50)==0
         alpha=alpha*0.3;
    end


    error=max(max(abs(d)));
    W=dw;
end
y=W;
end

function [alpha, beta, mu, delta] = nigpar(m, v, s, k)
%NIGPAR Parameters for the Normal-Inverse-Gaussian distribution.
%   [ALPHA, BETA, MU, DELTA] = NIGPAR(M, V, S, K) returns the scale 
%   paramter ALPHA, BETA which determines the skewness, the location 
%   parameter MU and the scale parameter DELTA of the Normal-Inverse-Gaussian 
%   distribution with mean M, variance V, skewness S and kurtosis K.
%
%   See also NIGPDF, NIGCDF, NIGINV, NIGRND, NIGSTATS.
%
%   References:
%   [1] Prause, K. (1999). The Generalized Hyperbolic Model

% -------------------------------------------------------------------------
% 
% Allianz, Group Risk Controlling
% Risk Methodology
% Koeniginstr. 28
% D-80802 Muenchen
% Germany
% Internet:    www.allianz.de
% email:       ralf.werner@allianz.de
% 
% Implementation Date:  2006 - 05 - 01
% Author:               Dr. Ralf Werner
%
% -------------------------------------------------------------------------

Tol = 1.0e-007;

%% Default values
if nargin < 1
    m = 0;
end
if nargin < 2
    v = 1;
end
if nargin < 4
    if nargin < 3
        s = 0;
        k = 3 + Tol;
    else
        k = 3 + 5/3*s^2 + Tol;
    end
end

%% Constraints for the parameters

if v <= 0
    error('The variance V must be positive.');
end

if (k - 5/3*s^2 - 3 <= 0)
    error('K must be greater than 3 + 5/3 S^2.');
end
    
alpha = sqrt((3*k - 4*s^2 - 9) / (v*(k-5/3*s^2 - 3)^2));
beta = s/(sqrt(v)*(k - 5/3*s^2 - 3));
mu = m - 3*s*sqrt(v)/(3*k - 4*s^2 - 9);
delta = 3^(3/2)*sqrt(v*(k - 5/3*s^2 - 3))/(3*k - 4*s^2 - 9);
end


function y=dldbM(x,B,lambda,pa)  
[n,T]=size(x);
de=zeros(n);
for i=1:T
    p=zeros(n,1);
    s=B*x(:,i);  %a vector of signals
    temp=pa(:,4).^2+(s-pa(:,3)).^2;
    p=pa(:,2)-(besselk(0,pa(:,1).*sqrt(temp))+besselk(2,pa(:,1).*sqrt(temp)))./besselk(1,pa(:,1).*sqrt(temp))/2.*(s-pa(:,3))./sqrt(temp).*pa(:,1)-(s-pa(:,3))./temp;
    de=de+p*(x(:,i)');
end
y=de/T+inv(B')-penalty(B,lambda); 
end


function y=penalty(W,lambda) %the penalty matrix used in gradient method
[~,n]=size(W);
p=zeros(n);
absw=abs(W);
for i=1:n
    for j=1:n
        if absw(i,j)>lambda
            p(i,j)=sign(W(i,j))*max(3.7*lambda-absw(i,j),0)/(3.7-1);
        else p(i,j)=sign(W(i,j))*lambda;
        end
    end
end
y=p;
end

function [nig_alpha,nig_beta,nig_mu,nig_delta]=nigem(sampoints,tol,maxiter,pa)
  
  % NIGEM A Maximum Likelihood parameter estimation of the Normal-Inverse
  %  Gaussian using an EM type algorithm.
  %  
  %   [NIG_ALPHA, NIG_BETA, NIG_MU, NIG_DELTA] = NIGEM(SAMPOINTS, TOL, MAXITER) returns the scale 
  %   parameter NIG_ALPHA, the skewness parameter NIG_BETA, the location 
  %   parameter NIG_MU and the scale parameter NIG_DELTA, given the one
  %   dimensional set of points SAMPOINTS. The algorithm terminates if
  %   the number of iterations exceed MAXITER or 
  %   $max_{i=1..4} |\theta_i(k+1)-\theta_i(k)/theta_i(k)| < TOL$, where $i$
  %   is the iteration number and $\theta$ is the set of all parameters.
  %    
  %   SAMPOINTS = row vector, TOL = scalar, MAXITER = scalar
  %
  % 
  %   References: 
  %      [1] Karlis, D. (2002) An EM type algorithm for maximum
  %      likelihood estimation of the Normal-Inverse Gaussian distribution
  %
  % -----------------------------------------------------------------------
  % Saket Sathe
  %  email: saket@ee.iitb.ac.in
  % -----------------------------------------------------------------------
    
    
  %
  % General statistics of the input points.
  %
  samp_stddev = std(sampoints); % Standerd Deviation
  samp_var = var(sampoints);    % Variance, second moment
  samp_mean = mean(sampoints);  % Mean
  samp_kurt = kurtosis(sampoints); % Kurtosis, fourth moment
  samp_skew = skewness(sampoints); % Skewness, third moment

  %
  % Initial estimates of alpha, beta, mu, delta and gamma.
  %
  nig_gamma1 = samp_skew/((samp_var)^(3/2)); % gamma_1
  nig_gamma2 = samp_kurt/((samp_var^2)-3); % gamma_2
  nig_gamma = sqrt(pa(1)^2-pa(2)^2); % gamma  
  nig_alpha = pa(1);
  nig_beta = pa(2);
  nig_mu = pa(3);
  nig_delta = pa(4);

  for iter=1:maxiter
    %disp(sprintf('Iteration: %f',iter));

    %
    % E-step
    %
    
    s_i = zeros(size(sampoints));
    w_i = zeros(size(sampoints));

    for i=1:size(sampoints,2)
      temp_phi = (1 + ((sampoints(i)-nig_mu)/nig_delta)^2)^0.5;
            
      if besselk(1,(nig_delta*nig_alpha*temp_phi))~=0 && nig_delta*nig_alpha*temp_phi>0
        s_i(1,i) = (nig_delta*temp_phi*besselk(0,(nig_delta*nig_alpha*temp_phi)))/(nig_alpha*besselk(1,(nig_delta*nig_alpha*temp_phi)));    
      else
        disp('zero')
        nig_delta*nig_alpha*temp_phi
      end
  
      if besselk(-1,(nig_delta*nig_alpha*temp_phi))~=0 && (nig_delta* ...
                                                      nig_alpha*temp_phi) >0
        w_i(1,i) = (nig_alpha*besselk(-2,(nig_delta*nig_alpha*temp_phi)))/(nig_delta*temp_phi*besselk(-1,(nig_delta*nig_alpha*temp_phi)));
      else
        disp('zero')
        nig_delta*nig_alpha*temp_phi
      end
      
    end
    % nig_gamma
    old_params = [nig_beta nig_delta nig_mu nig_alpha];

    % M-step

    M_hat = sum(s_i,2)/size(s_i,2);
    %lamda_cap_hat = size(s_i,2)*(sum(((w_i-(M_hat^-1)),2))^-1;

    tempvar = w_i-(1/M_hat);

    lamda_cap_hat = size(s_i,2)/(sum(tempvar,2));

    nig_delta = lamda_cap_hat^0.5;

    nig_gamma = nig_delta/M_hat;

    nig_beta = (dot(sampoints,w_i) - samp_mean*sum(w_i,2))/(size(s_i,2)- ...
                                                  M_hat*sum(w_i,2));

    nig_mu = samp_mean - nig_beta*M_hat;

    nig_alpha = (nig_gamma^2+nig_beta^2 )^0.5;

    % nig_gamma
    new_params=[nig_beta nig_delta nig_mu nig_alpha];

    paramdiff = new_params-old_params;
    perdiff = paramdiff./old_params;

    if max(abs(perdiff))<tol
      nig_alpha;
      nig_beta;
      nig_mu;
      nig_delta;
      break
    end

  end

end



