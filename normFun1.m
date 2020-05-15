function [M]=normFun1(M) 
for i=1:length(M)         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   z=find(M(:,i)~=0);
   M(z,i)=M(z,i)./sum(M(z,i));
end  


