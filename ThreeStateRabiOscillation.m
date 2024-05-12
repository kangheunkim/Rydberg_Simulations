dm780 = 1.689538049319012e-29;
dm480 = -2.40325975371874e-32;
% calculated from arc in python

epsilon_0 = 8.8541878128e-12;
hbar = 1.0545718176461565e-34;
c = 299792458;

pow480 = 500e-3;
pow780 = 150e-6;
w0480 = 30e-6;
w0780 = 100e-6;

Delta_I = 2*pi*660e6;
Gamma_I = 2*pi*6e6;

Om780 = sqrt(4*pow780/epsilon_0/c/pi/w0780^2)/hbar*dm780;
Om480 = sqrt(4*pow480/epsilon_0/c/pi/w0480^2)/hbar*dm480;

acstark = (Om780^2/Delta_I-Om480^2/Delta_I);

detdd = acstark/2;

H_three = 0.5*[-detdd Om480 0; Om480 Delta_I Om780; 0 Om780 detdd];
L_three = [0 0 0; 0 0 0; 0 sqrt(Gamma_I) 0];

stepnum = 5e3;
tend = 2e-6; dt = tend/stepnum;
rhoi = zeros(3,3); rhoi(end,end) = 1;
rho = zeros(3,3,stepnum);
prob = zeros(3,stepnum);
lindcol = L_three;
lindcolctp = ctranspose(lindcol);

H1 = H_three; H2 = H_three; H3 = H_three;
rho(:,:,1) = rhoi; prob(:,1) = abs(diag(rho(:,:,1))).^2;

for tstep = 2:stepnum
    disp(['tstep : ' num2str(tstep) '/' num2str(stepnum)]);
    
    rhob = rho(:,:,tstep-1); rhobinit = rhob;
    
    k1 = -1i*H1*rhob + 1i*rhob*H1...
        + lindcol*rhob*lindcolctp - 1/2*lindcolctp*lindcol*rhob - 1/2*rhob*lindcolctp*lindcol;
    
    rhob = rhobinit + dt/2*k1;
    
    k2 = -1i*H2*rhob + 1i*rhob*H2...
        + lindcol*rhob*lindcolctp - 1/2*lindcolctp*lindcol*rhob - 1/2*rhob*lindcolctp*lindcol;
    
    rhob = rhobinit + dt/2*k2;
    
    k3 = -1i*H2*rhob + 1i*rhob*H2...
        + lindcol*rhob*lindcolctp - 1/2*lindcolctp*lindcol*rhob - 1/2*rhob*lindcolctp*lindcol;
    
    rhob = rhobinit + dt*k3;
    
    k4 = -1i*H3*rhob + 1i*rhob*H3...
        + lindcol*rhob*lindcolctp - 1/2*lindcolctp*lindcol*rhob - 1/2*rhob*lindcolctp*lindcol;
    
    rho(:,:,tstep) = rhobinit + 1/6*dt*(k1+2*k2+2*k3+k4);
    prob(:,tstep) = abs(diag(rho(:,:,tstep)));
    
end

disp(Om780/2/pi/1e6)
disp(Om480/2/pi/1e6)

disp(Om780*Om480/Delta_I/2/pi/1e6)

figure; plot((1:stepnum)*dt,prob);


