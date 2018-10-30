function [MZp] = Zp(Ph,Vpol)

Vp = zeros(size(Vpol));
Vp(:,1) = Ph(:,2)-Vpol(:,1);
Vp(:,2) = Ph(:,3)-Vpol(:,2);
Vp(:,3) = Ph(:,4)-Vpol(:,3);
Vp(:,4) = (Ph(:,2)-Ph(:,3))-Vpol(:,4);
Vp(:,5) = (Ph(:,3)-Ph(:,4))-Vpol(:,5);
Vp(:,6) = (Ph(:,4)-Ph(:,2))-Vpol(:,6);

MZp = zeros(size(Vp));
ires = Ph(:,5)+Ph(:,6)+Ph(:,7);
for k = 1:3
    MZp(:,k) = Vp(:,k)./(Ph(:,k+4)+k0*ires);
end
MZp(:,4) = Vp(:,4)./(Ph(:,5)-Ph(:,6));
MZp(:,5) = Vp(:,5)./(Ph(:,6)-Ph(:,7));
MZp(:,6) = Vp(:,6)./(Ph(:,7)-Ph(:,5));