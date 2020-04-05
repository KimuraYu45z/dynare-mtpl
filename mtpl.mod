var A Y K L r w C I G T tau pi MS_P pitilde wtilde zeta1 zeta2;
varexo epsA epsG epstau;

parameters alpha delta rho phi eta psi gammastar taustar;
alpha = 0.5;
delta = 0.01;
rho = 0.01;
phi = 0.5;
eta = 1;
psi = 0.5;
gammastar = 0.3;
taustar = 0.25;

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

parameters Istar Gstar Tstar MS_Pstar zeta1star zeta2star pistar;
Istar = delta * Kstar;
Gstar = gammastar * Ystar;
Tstar = taustar * Ystar;
mustar = 0.01;
pistar = 0.01;
MS_Pstar = Tstar / rstar;
zeta1star = Gstar - Tstar + MS_Pstar;
zeta2star = (1 + rstar) * MS_Pstar - Tstar;
pistar = zeta1star / MS_Pstar - 1;

model(linear);
A = phi * A(-1) + epsA;
Y = A + alpha * K + (1 - alpha) * L;
r = A + (alpha - 1) * K + (1-alpha) * L;
w = A + alpha * K - alpha * L;
Y = (Cstar * C + Istar * I + Gstar * G) / Ystar;
K = (1 - delta) * K(-1) + delta * I(-1);
G = phi * G(-1) + epsG;
T = tau + Y;
tau = phi * tau(-1) + epstau;
wtilde = C + eta * L;
C(+1) - C = rstar * r(+1) / (rstar - delta + 1);
MS_P(+1) = zeta1 - pistar * pi / (1 + pistar);
pistar * pitilde / (1 + pistar) = zeta1 - zeta2;
zeta1 = (Gstar * G - Tstar * T + MS_Pstar * MS_P) / zeta1star;
zeta2 = ((1 + rstar) * MS_Pstar * (rstar * r / (1 + rstar) + MS_P) - Tstar * T) / zeta2star;
pi = psi * pi(-1) + (1 - psi) * pitilde;
w = pistar * pi(-1) / (1 + pistar) - pistar * pitilde(-1) / (1 + pistar) + wtilde;
end;

shocks;
var epsA = 0.01;
var epsG = 0.01;
var epstau = 0.01;
end;

stoch_simul Y K L r w C I pi;
