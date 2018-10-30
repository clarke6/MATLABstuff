function [M] = MHOout(Ph,Vpol,Zr,k0)

S1 = zeros(size(Vpol));
ires = Ph(:,5)+Ph(:,6)+Ph(:,7);
S1(:,1) = (Ph(:,5)+k0*ires)*Zr - Ph(:,2);
S1(:,2) = (Ph(:,6)+k0*ires)*Zr - Ph(:,3);
S1(:,3) = (Ph(:,7)+k0*ires)*Zr - Ph(:,4);
S1(:,4) = (Ph(:,5)-Ph(:,6))*Zr - (Ph(:,2)-Ph(:,3));
S1(:,5) = (Ph(:,6)-Ph(:,7))*Zr - (Ph(:,3)-Ph(:,4));
S1(:,6) = (Ph(:,7)-Ph(:,5))*Zr - (Ph(:,4)-Ph(:,2));

M = zeros(size(Vpol));
for n = 1:length(Vpol(1,:))
    for k = 1:length(Vpol(:,1))
        if real(S1(k,n)*conj(Vpol(k,n))) > 0
            M(k,n) = 1;
        end
    end
end
end