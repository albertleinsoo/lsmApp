function [c0Out, c1Out, GOut, b0eOut, b1eOut, aeOut, sigOut, sigB0eOut, sigB1eOut, b0eCIOut, b1eCIOut, sigAeOut, aeCIOut, aAveOut, sigAOut, nOut] = lsmFunc(b0,b1,th,a,ph0,ph1,phi0Str,phi1Str,Th0,Th1,dt)
%% Least Square Method for calculating model parameters. Two parameter case
%
%Input params
%
% b0        - param initial value
% b1        - param initial value
% th        - vector of θ (θῆτα) 
% a         - vector of a (άλφα) parameters
% ph0       - ɸ0(θ) (φεῖ) function handler
% ph1       - ɸ1(θ) (φεῖ) function handler
% phi0Str   - string name for ɸ0(θ)
% phi1Str   - string name for ɸ1(θ)
% Th0       - lower theta boundary
% Th1       - upper theta boundary
% dt        - theta step
%
%Output params
% 
% Model parameters:
%   c0Out - parameterfor calculatioin estimated b0   
%   c1Out - parameterfor calculatioin estimated b1 
%   GOut - matrix for calculation b0, b1 params
%   b0eOut - estimated b0 model parameter
%   b1eOut - estimated b1 model parameter
%   aeOut - estimated a (άλφα) parameters
%   sigOut - standard deviation (Root mean square - RMS) of measurement errors
%   sigB0eOut - standard deviation error of estimated b0 parameter
%   sigB1eOut - standard deviation error of estimated b1 parameter
%   b0eCIOut - confidence interval for b0
%   b1eCIOut - confidence interval for b1
%   sigAeOut - array of σ  a^ (άλφα) estimated
%   aeCIOut - confidence intervals for a (άλφα) estimated parameters
%   aAveOut -  average a (άλφα)
%   sigAOut -  sigma for measure model uncertainty
%   nOut - model uncertainty
%
% Throws error if ɸ(θ) is infinity or Nan
%

%   Albert Leinsoo 2023.

%% 1. Calculation of model parameter estimates (b0e b1e)
% N - 6 for two parameter case
N = 6;

% Two parameter case allows to use 2 c coefficients.
c0 = 0;
c1 = 0;

% While calculating c0, c1, check functions on infinity and Nan
for i = 1:length(th)
    if (isinf(ph0(th(i))) || isnan(ph0(th(i))))
        error('\nɸ0(θ) is Infinity or NaN . Check ɸ0(θ) at θ = %g',phi0Str, th(i))
    end
    if (isinf(ph1(th(i))) || isnan(ph1(th(i))))
        error('\nɸ1(θ) is Infinity or NaN. Check ɸ1(θ) at θ = %g',phi1Str, th(i))
    end

    c0 = c0 + a(i)*ph0(th(i));
    c1 = c1 + a(i)*ph1(th(i));
end

% print с0, с1 to terminal
c0
c1

% initialize G matrix
G = zeros(2,2);

% calculate G matrix
for i = 1:length(th)
    G(1,1) = G(1,1) + power(ph0(th(i)),2);
    G(2,2) = G(2,2) + power(ph1(th(i)),2);
    G(1,2) = G(1,2) + ph0(th(i))*ph1(th(i));
end
G(2,1) = G(1,2);

% print G matrix
G

% b estimated params
b0e = (G(2,2)*c0 - G(1,2)*c1) / (G(1,1)*G(2,2) - G(2,1)*G(1,2))
b1e = (-G(2,1)*c0 + G(1,1)*c1) / (G(1,1)*G(2,2) - G(2,1)*G(1,2))

%% 2. Calculation of estimates of standard deviation (Root mean square - RMS) of measurement errors
% ae - a (άλφα) estimated
ae = b0e*ph0(th) + b1e*ph1(th)

sig = 0;
for i = 1:length(th)
    sig = sig + (a(i) - ae(i))^2;
end
sig = sqrt(sig / N)

%% 3. Calculation of estimates of the standard deviation of errors in estimating model parameters
sigB0e = sig * sqrt(G(2,2) / (G(1,1)*G(2,2)-G(2,1)*G(1,2)))
sigB1e = sig * sqrt(G(1,1) / (G(1,1)*G(2,2)-G(2,1)*G(1,2)))

%% 4. Determining the boundaries of confidence intervals for parameter estimates using the “two sigma rule”
% confidence interval for b estimated b0e b1e
b0eCI = [b0e - 2*sigB0e, b0e + 2*sigB0e]
b1eCI = [b1e - 2*sigB1e, b1e + 2*sigB1e]

%% 5. Significance check of b0e b1e params

% If b0 confidence interval contains 0, warn user that it may not be
% significant
if discretize(0, b0eCI) == 1
    warning('Confidence interval of b0e contains Zero. This coef should be excluded from the model.')
end

% If b1 confidence interval contains 0, warn user that it may not be
% significant
if discretize(0, b1eCI) == 1
    warning('Confidence interval of b1e contains Zero. This coef should be excluded from the model.')
end

% If bouth b0 and b1 confidence intervals contain 0, warn user that they may not be
% significant
if ( (discretize(0, b0eCI) == 1) && (discretize(0, b1eCI) == 1) )
    warning('Confidence interval of b0e and b1e contains Zero')
end

%% 6. Determining the boundaries of confidence intervals for a (άλφα) estimated
sigAe = @(thi) sig * sqrt( ( (ph0(thi).^2)*G(2,2) + (ph1(thi).^2)*G(1,1) - 2*ph0(thi).*ph1(thi)*G(1,2) ) / ...
    ( G(1,1)*G(2,2)-G(2,1)*G(1,2) ) )

% Confidence intervals matrix of ae
aeCI = [(ae - 2*sigAe(th))', (ae + 2*sigAe(th))']

%% 7. Determining the measure of model uncertainty
% n (ἦτα) - model uncertainty
% aAve - a average
aAve = sum(a)/N
sigA = sqrt((sum((a-aAve).^2)) / 5)

n = sig / sigA

%% 8. Output parameters

c0Out = c0;
c1Out = c1;
GOut = G;
b0eOut = b0e;
b1eOut = b1e;
aeOut = ae;
sigOut = sig;
sigB0eOut = sigB0e;
sigB1eOut = sigB1e;
b0eCIOut= b0eCI;
b1eCIOut= b1eCI;
sigAeOut = sigAe(th);
aeCIOut = aeCI;
aAveOut = aAve;
sigAOut = sigA;
nOut = n;

%% 9. Plot

% α(θ) - initial function
aTh = @(thi) b0*ph0(thi) + b1*ph1(thi);

% α^(θ) - estimated function
aeTh = @(thi) b0e*ph0(thi) + b1e*ph1(thi);

% α^(θ) + 2 σ
aeP2sig = @(thi) aeTh(thi) + 2*sigAe(thi);

% α^(θ) - 2 σ
aeM2sig = @(thi) aeTh(thi) - 2*sigAe(thi);

% Last param sets position and size proportional to screen
f = figure('Name','Least Square Method', 'units','normalized','outerposition',[0.4 0.4 0.6 0.6]);

hold on;
 
% functions names  names on plot
phi0FuncName = strcat('\phi_{0}(\theta) = ', phi0Str);
phi1FuncName = strcat('\phi_{1}(\theta) = ', phi1Str);

% '$' added for latex format to show text symbols on plot
plotTitle = strcat('$',phi0FuncName, '\hspace{1cm}', phi1FuncName, '$');
title(plotTitle,'interpreter','latex', 'FontSize', 25)
xlabel('θ (θήτα)', 'FontSize',25) 
ylabel('α (άλφα)', 'FontSize',25) 

%initial values points
plot(th,a,'o','LineWidth',2,'Color',[0.5,0.5,0.5], 'DisplayName', '$\alpha$ values');

% 2 σα^ error in α point
errorbar(th,a, repmat(2*sig,length(th),1),'Color',[0, 0.6, 1], 'LineStyle','none', 'DisplayName','2$\hat{\sigma}_{\epsilon}$');

% Xi Yi for smooth plots
Xi = Th0:dt:Th1;
Yi = aTh(Xi);

% α(θ)
plot(Xi,Yi, 'DisplayName','$\alpha(\theta)$');

% α^(θ)
Yi = aeTh(Xi);
plot(Xi,Yi,'Color',[0, 0, 0], 'DisplayName','$\hat{\alpha}({\theta})$');

% α^(θ) + 2 σ
Yi = aeP2sig(Xi);
plot(Xi,Yi,'Color',[0, 0.5, 1], 'DisplayName','$\hat{\alpha}({\theta}) + 2\hat{\sigma}_{\hat\alpha}$');

% α^(θ) - 2 σ
Yi = aeM2sig(Xi);
plot(Xi,Yi,'Color',[1, 0, 0], 'DisplayName','$\hat{\alpha}({\theta}) - 2\hat{\sigma}_{\hat\alpha}$');

% legend Interpreter of LaTeX markup 
legend('Interpreter','latex','FontSize',20);
grid on;
hold off;
end