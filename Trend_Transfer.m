function  [a1,reg_a1,v1,d1,theta1] = Trend_Transfer(PC,EOF,V);

[ntimes,np] = size(PC);
[nx,ny ,np] = size(EOF);

d = nan(np-1,1);

a1     = PC;
reg_a1 = EOF;
v1     = V;

k = 0;
for i = 2:np
    k = k+1;
    
    pcr(:,1) = a1(:,1) - detrend(a1(:,1));
    pcr(:,2) = a1(:,i) - detrend(a1(:,i));
    pcr2 = diff(pcr,[],1);
    pcr2 = pcr2(1,:);
    d(k) = pcr2(2)/(pcr2(1)+eps);
        
    ap1 = 1/sqrt((1+d(k)^2)); ap2 =d(k)/sqrt((1+d(k)^2));
    ap12 = ap1^2;             ap22 = ap2^2;
    
    pc1 = ap1*a1(:,1) + ap2*a1(:,i);
    pc2 = ap1*a1(:,i) - ap2*a1(:,1);
    a1(:,1) = pc1;
    a1(:,i) = pc2;
    
    aa1 = ap1*reg_a1(:,:,1) + ap2*reg_a1(:,:,i);
    aa2 = ap1*reg_a1(:,:,i) - ap2*reg_a1(:,:,1);   
    reg_a1(:,:,1) = aa1;
    reg_a1(:,:,i) = aa2;
    
    V1 = ap12*v1(1) + ap22*v1(i);
    V2 = ap12*v1(i) + ap22*v1(1);   
    v1(1) = V1;
    v1(i) = V2;
end
d1 = d;
theta1 = atan(-d1)*180/pi;