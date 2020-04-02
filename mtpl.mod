var A Y K L r w C I G T gamma tau mu pi MS_P pitilde zeta wtilde;
varexo epsA epsgamma epstau epsmu;

parameters alpha delta rho phi eta psi gammastar taustar;
alpha = 0.5;
delta = 0.01;
rho = 0.01;
phi = 0.5;
eta = 1;
psi = 0.5;
gammastar = 0.5;
taustar = 0.5;

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

parameters Istar Gstar Tstar mustar pistar MS_Pstar zetastar;
Istar = delta * Kstar;
Gstar = gammastar * Ystar;
Tstar = taustar * Ystar;
mustar = 0.01;
pistar = 0.01;
MS_Pstar = Tstar / rstar;
zetastar = (1 + rstar) * MS_Pstar + Tstar;
psi = 0.5;

model(linear);
A = phi * A(-1) + epsA;
Y = A + alpha * K + (1 - alpha) * L;
r = A + (alpha - 1) * K + (1-alpha) * L;
w = A + alpha * K - alpha * L;
Y = (Cstar * C + Istar * I + Gstar * G) / Ystar;
K = (1 - delta) * K(-1) + delta * I(-1);
G = gamma + Y;
T = tau + Y;
gamma = phi * gamma(-1) + epsgamma;
tau = phi * tau(-1) + epstau;
wtilde = C + eta * L;
C(+1) - C = rstar / (rstar - delta + 1) * r(+1);
mu = phi * mu(-1) + epsmu;
MS_P = mustar / (1 + mustar) * mu(-1) - pistar / (1 + pistar) * pi(-1) + MS_P(-1);
pistar / (1 + pistar) * pitilde = mustar / (1 + mustar) * mu + MS_P - zeta;
zeta = (1 + rstar) * MS_Pstar / zetastar * (rstar / (1 + rstar) * r + MS_P) + Tstar / zetastar * T;
pistar / (1 + pistar) * pi = psi * pistar / (1 + pistar) * pi(-1) + (1 - psi) * pistar / (1 + pistar) * pitilde;
w = pistar / (1 + pistar) * pi(-1) - pistar / (1 + pistar) * pitilde(-1) + wtilde;
end;

shocks;
var epsA = 0.01;
var epsgamma = 0.01;
var epstau = 0.01;
var epsmu = 0.01;
end;

stoch_simul Y K L r w C I pi;
