 clc;clear;
warning('off');
%0.8632
currentFolder = pwd;              
addpath(genpath(currentFolder));  

Wdd = load('diseasesimilarity.txt');
A = load('known_lncRNA_disease_interaction.txt');
Wrr_v=lncRNAfunsim(Wdd,A);
A_ori=A;


y_train=WKNKN( A, Wrr_v, Wdd, 5, 1 ); 
similairty_matrix2=LNS(y_train',0,50,'regulation2');    

        
similairty_matrix1=LNS(y_train,0,90,'regulation2');    

F_1_ori=BR(y_train,similairty_matrix1,similairty_matrix2,5,4,0.2);
 
 F_1_ori_ori=F_1_ori;
 index=find(A_ori==1);


auc = zeros(1,100);
for i = 1:100
    i
    indices = crossvalind('Kfold', length(index), 5);
    A = A_ori;
     F_1_ori=F_1_ori_ori;
   
    for cv = 1:5
       cv;
       index_2 = find(cv == indices);
      
       A(index(index_2)) = 0;

  Wrr_v=lncRNAfunsim(Wdd,A);
y_train=WKNKN( A, Wrr_v, Wdd, 5, 1 ); 

similairty_matrix2=LNS(y_train',0,50,'regulation2');    

        
similairty_matrix1=LNS(y_train,0,90,'regulation2');    


F_1=BR(y_train,similairty_matrix1,similairty_matrix2,5,4,0.2);
      
      F_1_ori(index(index_2)) = F_1(index(index_2));
       A = A_ori;
    end

  
    pre_label_score = F_1_ori(:);
    label_y = A_ori(:);
    auc(i) = roc_1(pre_label_score,label_y,'red');
end
auc_ave = mean(auc);
auc_std = std(auc);
 

