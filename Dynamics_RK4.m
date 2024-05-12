tic
numatoms = 2;
rhoi = zeros(2^numatoms,2^numatoms); rhoi(2^numatoms,2^numatoms) = 1;

px = [0 1; 1 0]; pz = [1 0; 0 -1]; py = [0 1i; -1i 0]; nz = [1 0; 0 0];

posX = 1:numatoms; posY = zeros(1,numatoms); 
% posX = [1 2 100]; posY = zeros(1,numatoms);

rab = 1e6; det = 0; lat = 8e-6; v = 1023000/(lat*1e6)^6*1e6;
stepnum = 1e3; tend = 5e-6; dt = tend/stepnum;
gammaind = 0.025e6*ones(1,numatoms);
gammacol = 0.100e6;

rho = zeros(2^numatoms,2^numatoms,stepnum);
prob = zeros(2^numatoms,stepnum);

distarr = zeros(numatoms,numatoms);
for ii = 1:numatoms; for jj = 1:numatoms; distarr(ii,jj) = norm([posX(ii),posY(ii)]-[posX(jj),posY(jj)]); end; end

hrab = zeros(2^numatoms,2^numatoms);
for idx = 1:numatoms
    seed = 1; 
    seed = kron(seed,eye(2^(idx-1)));
    seed = kron(seed,px); 
    seed = kron(seed,eye(2^(numatoms-idx)));
    hrab = hrab + seed;
end

hdet = zeros(2^numatoms,2^numatoms);
lindind = zeros(2^numatoms,2^numatoms,numatoms);

for idx = 1:numatoms
    seed = 1; 
    seed = kron(seed,eye(2^(idx-1)));
    seed = kron(seed,pz); 
    seed = kron(seed,eye(2^(numatoms-idx)));
    hdet = hdet + seed;
    lindind(:,:,idx) = sqrt(2*pi*gammaind(idx)/2)*seed;
end

lindcol = sqrt(2*pi*gammacol/2)*hdet;

hv = zeros(2^numatoms,2^numatoms);
for idx = 1:numatoms
    for idxidx = idx+1:numatoms
        seed = 1;
        seed = kron(seed,eye(2^(idx-1)));
        seed = kron(seed,nz);
        seed = kron(seed,eye(2^(idxidx-idx-1)));
        seed = kron(seed,nz);
        seed = kron(seed,eye(2^(numatoms-idxidx)));
        hv = hv + 1/distarr(idx,idxidx)^6*seed;
    end
end

rabarr = rab*ones(1,2*stepnum);
detarr = det*ones(1,2*stepnum);
varr = v*ones(1,2*stepnum);

rho(:,:,1) = rhoi; prob(:,1) = abs(diag(rho(:,:,1))).^2;

for tstep = 2:stepnum
    
    disp(['tstep : ' num2str(tstep) '/' num2str(stepnum)]);
    
    rhob = rho(:,:,tstep-1); rhobinit = rhob;
    
    H1 = 2*pi*(0.5*rabarr(2*(tstep-1)-1)*hrab - 0.5*detarr(2*(tstep-1)-1)*hdet + varr(2*(tstep-1)-1)*hv);
    H2 = 2*pi*(0.5*rabarr(2*(tstep-1))*hrab - 0.5*detarr(2*(tstep-1))*hdet + varr(2*(tstep-1))*hv);
    H3 = 2*pi*(0.5*rabarr(2*(tstep)-1)*hrab - 0.5*detarr(2*(tstep)-1)*hdet + varr(2*(tstep)-1)*hv);

    k1 = -1i*H1*rhob + 1i*rhob*H1...
        + lindcol*rhob*lindcol - 1/2*lindcol*lindcol*rhob - 1/2*rhob*lindcol*lindcol...
        + sum(prod3arr(prod3arr(lindind,repmat(rhob,1,1,numatoms)),lindind),3)...
        - 0.5*sum(prod3arr(prod3arr(repmat(rhob,1,1,numatoms),lindind),lindind),3)...
        - 0.5*sum(prod3arr(prod3arr(lindind,lindind),repmat(rhob,1,1,numatoms)),3);
    
    rhob = rhobinit + dt/2*k1;
    
    k2 = -1i*H2*rhob + 1i*rhob*H2...
        + lindcol*rhob*lindcol - 1/2*lindcol*lindcol*rhob - 1/2*rhob*lindcol*lindcol...
        + sum(prod3arr(prod3arr(lindind,repmat(rhob,1,1,numatoms)),lindind),3)...
        - 0.5*sum(prod3arr(prod3arr(repmat(rhob,1,1,numatoms),lindind),lindind),3)...
        - 0.5*sum(prod3arr(prod3arr(lindind,lindind),repmat(rhob,1,1,numatoms)),3);
    
    rhob = rhobinit + dt/2*k2;
    
    k3 = -1i*H2*rhob + 1i*rhob*H2...
        + lindcol*rhob*lindcol - 1/2*lindcol*lindcol*rhob - 1/2*rhob*lindcol*lindcol...
        + sum(prod3arr(prod3arr(lindind,repmat(rhob,1,1,numatoms)),lindind),3)...
        - 0.5*sum(prod3arr(prod3arr(repmat(rhob,1,1,numatoms),lindind),lindind),3)...
        - 0.5*sum(prod3arr(prod3arr(lindind,lindind),repmat(rhob,1,1,numatoms)),3);
    
    rhob = rhobinit + dt*k3;
    
    k4 = -1i*H3*rhob + 1i*rhob*H3...
        + lindcol*rhob*lindcol - 1/2*lindcol*lindcol*rhob - 1/2*rhob*lindcol*lindcol...
        + sum(prod3arr(prod3arr(lindind,repmat(rhob,1,1,numatoms)),lindind),3)...
        - 0.5*sum(prod3arr(prod3arr(repmat(rhob,1,1,numatoms),lindind),lindind),3)...
        - 0.5*sum(prod3arr(prod3arr(lindind,lindind),repmat(rhob,1,1,numatoms)),3);
    
    rho(:,:,tstep) = rhobinit + 1/6*dt*(k1+2*k2+2*k3+k4);
    prob(:,tstep) = abs(diag(rho(:,:,tstep)));
end

toc
figure; plot((1:stepnum)*dt,prob); legend(flip(dec2bin(0:2^numatoms-1))); legend('off')

function C = prod3arr(A,B)
arrsize = size(A);
if numel(arrsize) == 2
    C = A*B;
else
    C = zeros(arrsize);
    for i = 1:arrsize(3)
        C(:,:,i) = A(:,:,i)*B(:,:,i);
    end
end
end