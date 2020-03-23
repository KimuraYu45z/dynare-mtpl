var A Y K L r w C I G T gamma tau MS_P;
varexo epsgamma epstau;

parameters alpha delta rho phi eta gammastar taustar;
alpha = 0.5;
delta = 0.01;
rho = 0.01;
phi = 0.9;
eta = 1;
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

parameters Istar Gstar Tstar MS_Pstar;
Istar = delta * Kstar;
Gstar = gammastar * Ystar;
Tstar = taustar * Ystar;
MS_Pstar = Tstar / rstar;

model(linear);
A = 0;
Y = A + alpha * K + (1 - alpha) * L;
r = A + (alpha - 1) * K + (1-alpha) * L;
w = A + alpha * K - alpha * L;
Y = (Cstar * C + Istar * I + Gstar * G) / Ystar;
K = (1 - delta) * K(-1) + delta * I(-1);
G = gamma + Y;
T = tau + Y;
gamma = phi * gamma(-1) + epsgamma;
tau = phi * tau(-1) + epstau;
w = C + eta * L;
C(+1) - C = rstar / (rstar - delta + 1) * r(+1);
MS_P(+1) = ((1 + rstar) * MS_Pstar * (rstar / (1 + rstar) * r + MS_P) - Tstar * T) / ((1 + rstar) * MS_Pstar - Tstar);
end;

shocks;
var epsgamma = 0.01;
var epstau = 0.01;
end;

stoch_simul Y K L r w C I MS_P;
