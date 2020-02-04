var A Y K L r w C I G T i gamma tau psi pi;
varexo epsgamma epstau;

parameters alpha delta rho phi eta gammastar taustar;
alpha = 0.5;
delta = 0.01;
rho = 0.01;
phi = 0.9;
eta = 1;
gammastar = 0.5;
taustar = 0.49;

parameters Astar rstar K_L Y_L C_L wstar Lstar Kstar Ystar Cstar;
Astar = 1;
rstar = rho + delta;
K_L = (rstar / alpha / Astar)^(1 / (alpha - 1));
Y_L = Astar * K_L^alpha;
C_L = (1 - gammastar) * Y_L - delta * K_L;
wstar = (1 - alpha) * Astar * K_L^alpha;
Lstar = (wstar / (eta + 1))^(1 / (eta + 1)) * C_L^(-1 / (eta + 1));
Kstar = K_L * Lstar;
Ystar = Y_L * Lstar;
Cstar = C_L * Lstar;

parameters Istar Gstar pistar psistar istar;
Istar = delta * Kstar;
Gstar = gammastar * Ystar;
a = taustar;
b = (rstar^2 - rstar + 2) * taustar + (1 - 2 * rstar) * gammastar;
c = rstar * (gammastar * rstar - taustar);
pistar = (-b + (b^2 - 4 * a * c)^(1/2)) / (2 * a);
psistar = (1 + pistar) * rstar * taustar / ((rstar + pistar) * ((2 + pistar) * taustar - gammastar));
istar = (1 + rstar) * (1 + pistar) - 1;

model(linear);
A = 0;
Y = A + alpha * K + (1 - alpha) * L;
r = A + (alpha - 1) * K + (1-alpha) * L;
w = A + alpha * K - alpha * L;
Y = (Cstar * C + Istar * I + Gstar * G) / Ystar;
K = (1 - delta) * K(-1) + delta * I(-1);
G = gamma + Y;
T = tau + Y;
i = rstar / (1 + rstar) * r + pistar / (1 + pistar) * pi;
gamma = phi * gamma(-1) + epsgamma;
tau = phi * tau(-1) + epstau;
psi = pistar / (1 + pistar) * pi + r + tau - i - (2 * taustar * tau + pistar * taustar * (pi + tau) - gammastar * gamma) / ((2 + pistar) * taustar - gammastar);
pistar / (1 + pistar) * pi + tau - i - psi = (((1 + pistar) * taustar / istar) * (pistar / (1 + pistar) * pi(-1) + tau(-1) - i(-1)) + taustar / psistar * (tau(-1) - psi(-1)) - 2 * taustar * tau(-1) + gammastar * gamma(-1)) / (((1 + pistar) / istar + 1 / psistar - 2) * taustar - gammastar);
w = C + eta * L;
C(+1) - C = rstar / (rstar - delta + 1) * r(+1);
end;

shocks;
var epsgamma = 0.01;
var epstau = 0.01;
end;

stoch_simul Y K L r w C I pi;
