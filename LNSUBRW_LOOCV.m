clear;
clc;
warning('off');
%0.8874  

currentFolder = pwd;              
addpath(genpath(currentFolder));  

Wdd = load('diseasesimilarity.txt');
A = load('known_lncRNA_disease_interaction.txt');
Wrr_v=lncRNAfunsim(Wdd,A);
A_ori=A;

y_train=WKNKN( A, Wrr_v, Wdd, 5, 1 ); 

similairty_matrix2=LNS(y_train',0,50,'regulation2');    

        
similairty_matrix1=LNS(y_train,0,90,'regulation2');    %lncSimilairty

F_1_ori=BR(y_train,similairty_matrix1,similairty_matrix2,5,4,0.2);

index=find(A_ori==1);
for u=1:length(index)

    u
     A=A_ori ;
    A(index(u))=0;


     

     Wrr_v=lncRNAfunsim(Wdd,A);

y_train=WKNKN( A, Wrr_v, Wdd, 5, 1 ); 

similairty_matrix2=LNS(y_train',0,50,'regulation2');    

        
similairty_matrix1=LNS(y_train,0,90,'regulation2');    



F_1=BR(y_train,similairty_matrix1,similairty_matrix2,5,4,0.2);


    F_1_ori(index(u))=F_1(index(u));
    
    %A=A_ori ;
end
    pre_label_score = F_1_ori(:);
  
    label_y = A_ori(:);
    auc=roc_1(pre_label_score,label_y,'red');
  