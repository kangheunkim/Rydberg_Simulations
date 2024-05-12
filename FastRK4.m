tic
numatoms = 2;
rhob = zeros(2^numatoms,2^numatoms); rhob(2^numatoms,2^numatoms) = 1;

px = [0 1; 1 0]; pz = [1 0; 0 -1]; py = [0 1i; -1i 0]; nz = [1 0; 0 0];

posX = 1:numatoms; posY = zeros(1,numatoms); 
% posX = [1 2 100]; posY = zeros(1,numatoms);

rab = 1e6; det = 0; lat = 8e-6; v = 1023000/(lat*1e6)^6*1e6;
stepnum = 1e3; tend = 5e-6; dt = tend/stepnum;
gammaind = 0.025e6*2*pi/2;
gammacol = 0.100e6*2*pi/2;

% rho = zeros(2^numatoms,2^numatoms,stepnum);
prob = zeros(2^numatoms,stepnum);

distarr = zeros(numatoms,numatoms);
for ii = 1:numatoms; for jj = 1:numatoms; distarr(ii,jj) = norm([posX(ii),posY(ii)]-[posX(jj),posY(jj)]); end; end

hrab = sparse(0);
for idx = 1:numatoms
    seed = sparse(1); 
    seed = kron(seed,eye(2^(idx-1)));
    seed = kron(seed,px); 
    seed = kron(seed,eye(2^(numatoms-idx)));
    hrab = hrab + seed;
end

pzlist = zeros(numatoms,2^numatoms);

for i = 1:numatoms
    pzlist(i,:) = (-ones(1,2^numatoms)).^(2^(numatoms-i)-1<rem(((1:2^numatoms)-1),2^(numatoms-i+1)));
end

hdetdiag = sum(pzlist,1);

% nnlist = zeros(numatoms*(numatoms-1)/2,2^numatoms);
hvdiag = sparse(0);
% cnt = 1;
for idx = 1:numatoms
    for idxidx = idx+1:numatoms
        seed = sparse(1);
        seed = kron(seed,eye(2^(idx-1)));
        seed = kron(seed,nz);
        seed = kron(seed,eye(2^(idxidx-idx-1)));
        seed = kron(seed,nz);
        seed = kron(seed,eye(2^(numatoms-idxidx)));
%         hv = hv + 1/distarr(idx,idxidx)^6*seed;
        hvdiag = hvdiag + 1/distarr(idx,idxidx)^6*diag(seed)';
%         nnlist(cnt,:) = diag(seed);
%         cnt = cnt + 1;
    end
end

% hv = diag(hv)

% hv

rabarr = rab*ones(1,2*stepnum);
detarr = det*ones(1,2*stepnum);
varr = v*ones(1,2*stepnum);

% rho(:,:,1) = rhoi; 
prob(:,1) = abs(diag(rhob)).^2;

for tstep = 2:stepnum
    
    disp(['tstep : ' num2str(tstep) '/' num2str(stepnum)]);
    
%     rhob = rho(:,:,tstep-1); 
    rhobinit = rhob;
    
    H1rab = 2*pi*(0.5*rabarr(2*(tstep-1)-1))*hrab;
    H1diags =  - 2*pi*0.5*detarr(2*(tstep-1)-1)*hdetdiag + 2*pi*varr(2*(tstep-1)-1)*hvdiag;
    
    H2rab = 2*pi*0.5*rabarr(2*(tstep-1))*hrab;
    H2diags = - 2*pi*0.5*detarr(2*(tstep-1))*hdetdiag + 2*pi*varr(2*(tstep-1))*hvdiag;
    
    H3rab = 2*pi*0.5*rabarr(2*(tstep)-1)*hrab;
    H3diags = - 2*pi*0.5*detarr(2*(tstep)-1)*hdetdiag + 2*pi*varr(2*(tstep)-1)*hvdiag;

    k1 = cal_k_fast(rhob,H1rab,H1diags,sqrt(gammacol)*hdetdiag,pzlist,gammaind,numatoms);
    
    rhob = rhobinit + dt/2*k1;
    
    k2 = cal_k_fast(rhob,H2rab,H2diags,sqrt(gammacol)*hdetdiag,pzlist,gammaind,numatoms);
    
    rhob = rhobinit + dt/2*k2;
    
    k3 = cal_k_fast(rhob,H2rab,H2diags,sqrt(gammacol)*hdetdiag,pzlist,gammaind,numatoms);
    
    rhob = rhobinit + dt*k3;
    
    k4 = cal_k_fast(rhob,H3rab,H3diags,sqrt(gammacol)*hdetdiag,pzlist,gammaind,numatoms);
        
%     rho(:,:,tstep) = rhobinit + 1/6*dt*(k1+2*k2+2*k3+k4);
%     prob(:,tstep) = real(diag(rho(:,:,tstep)));

        rhob = rhobinit + 1/6*dt*(k1+2*k2+2*k3+k4);
    prob(:,tstep) = real(diag(rhob));

end

toc
figure; plot((1:stepnum)*dt,prob); legend(flip(dec2bin(0:2^numatoms-1))); legend('off')

function k = cal_k_fast(rho,Hrab,Hdiags,lindcoldiag,pzlist,gammaind,Nqub) % Hdiags is 1 times 2^numatoms matrix
    k = -1i*Hrab*rho + 1i*rho*Hrab + -1i*Hdiags'.*rho + 1i*rho.*Hdiags ...
        + lindcoldiag'.*rho.*lindcoldiag - 0.5*(lindcoldiag.^2)'.*rho - 0.5*rho.*(lindcoldiag.^2) ...
        - gammaind*Nqub*rho; % 0.5*gammaind*Nqub*rhob times 2
    for i = 1:Nqub
        k = k + gammaind*pzlist(i,:)'.*rho.*pzlist(i,:);
    end
end