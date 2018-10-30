function [MZ] = Zapp(M,k0)
Ires = M(:,5)+M(:,6)+M(:,7);
Zag = M(:,2)./(M(:,5)+k0*Ires);
Zbg = M(:,3)./(M(:,6)+k0*Ires);
Zcg = M(:,4)./(M(:,7)+k0*Ires);
Zab = (M(:,2)-M(:,3))./(M(:,5)-M(:,6));
Zbc = (M(:,3)-M(:,4))./(M(:,6)-M(:,7));
Zca = (M(:,4)-M(:,2))./(M(:,7)-M(:,5));
MZ = [Zag,Zbg,Zcg,Zab,Zbc,Zca];
end