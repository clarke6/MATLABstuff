function [MV] = Vpol(M)

j=sqrt(-1);
MV = zeros(length(M),6);
MV(:,1) = (M(:,3)-M(:,4))*j;
MV(:,2) = (M(:,4)-M(:,2))*j;
MV(:,3) = (M(:,2)-M(:,3))*j;
MV(:,4) = M(:,4)*-j;
MV(:,5) = M(:,2)*-j;
MV(:,6) = M(:,3)*-j;
end