function [E, C]=taskNetwork(Y,corthreshold)
          
[nV]=size(Y,2);
          
C=corrcoef(Y);
corthreshold = mean(C(:));
C(abs(C)<corthreshold)=0;
          
UC=triu(C,1); 

nzUC=find(UC~=0);
[E1,E2]=ind2sub([nV,nV],nzUC);
E=[E1,E2];


          
