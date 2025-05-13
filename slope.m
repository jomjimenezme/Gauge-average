 close all; clc;
tau = 0:max_tau;   

% 
logC = log10(C_tau);          %  log(C);

% Least-squares straight-line fit in log space ----------
p = polyfit(tau, logC, 1);  % p(1) = slope, p(2) = intercept
slopeLog = p(1);

%  Display 
fprintf('Slope in log10 space \t %f\n', slopeLog);
fprintf('value of slope: \t %f\n', slopeLog*log(10));

% -display fit -------------------------------------
figure
semilogy(tau, C_tau, 'o', 'MarkerFaceColor', 'b', 'DisplayName', 'data'); hold on
Cfit = 10.^(polyval(p, tau));   % back-transform to original scale
semilogy(tau, Cfit, '-r', 'DisplayName', ...
         sprintf('fit (slope = %.3g)', slopeLog));
xlabel('\tau'); ylabel('C(\tau)'); grid on
legend('Location','best')