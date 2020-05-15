clear;
clc;
warning('off');


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

[~,disease]=xlsread(['disease_178.xlsx']);
[~,lncRNA]=xlsread(['lncRNA_115.xlsx']);
%% ×îÖÕ½á¹û
allresult(lncRNA,disease,A,F_1_ori);  
