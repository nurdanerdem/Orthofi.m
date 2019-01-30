clear 
close all

%write the full path of IT list file for species 1
Spe1= readtable('/home/nerdem/IT Project/IT_list_Salmonella_vs_Bacillus_conservative.xlsx','Sheet','Salmonella');

%write the full path of IT list file for species 2
Spe2= readtable('/home/nerdem/IT Project/IT_list_Salmonella_vs_Bacillus_conservative.xlsx','Sheet','Bacillus');

%write the names of the species
species1='Salmonella'; 
species2='Bacillus';

%extract upstream gene name, IT sequence and downstream gene name from
%positive strand for species 1
pos_UpGeneName_Spe1=Spe1(find(strcmp('+',table2array(Spe1(:,4)))),8);
pos_IT_Spe1=Spe1(find(strcmp('+',table2array(Spe1(:,4)))),7);
pos_ITname_Spe1=Spe1(find(strcmp('+',table2array(Spe1(:,4)))),1);
pos_DwGeneName_Spe1=Spe1(find(strcmp('+',table2array(Spe1(:,4)))),12);

%extract upstream gene name, IT sequence and downstream gene name from
%negative strand for species 1
neg_UpGeneName_Spe1=Spe1(find(strcmp('-',table2array(Spe1(:,4)))),8);
neg_IT_Spe1=Spe1(find(strcmp('-',table2array(Spe1(:,4)))),7);
neg_ITname_Spe1=Spe1(find(strcmp('-',table2array(Spe1(:,4)))),1);
neg_DwGeneName_Spe1=Spe1(find(strcmp('-',table2array(Spe1(:,4)))),12);


%do the same for species 2
pos_UpGeneName_Spe2=Spe2(find(strcmp('+',table2array(Spe2(:,4)))),8);
pos_IT_Spe2=Spe2(find(strcmp('+',table2array(Spe2(:,4)))),7);
pos_ITname_Spe2=Spe2(find(strcmp('+',table2array(Spe2(:,4)))),1);
pos_DwGeneName_Spe2=Spe2(find(strcmp('+',table2array(Spe2(:,4)))),12);

neg_UpGeneName_Spe2=Spe2(find(strcmp('-',table2array(Spe2(:,4)))),8);
neg_IT_Spe2=Spe2(find(strcmp('-',table2array(Spe2(:,4)))),7);
neg_ITname_Spe2=Spe2(find(strcmp('-',table2array(Spe2(:,4)))),1);
neg_DwGeneName_Spe2=Spe2(find(strcmp('-',table2array(Spe2(:,4)))),12);


 

[p1p2 up_p1p2]= get_orthologs(pos_UpGeneName_Spe1,pos_UpGeneName_Spe2, pos_DwGeneName_Spe1, pos_DwGeneName_Spe2, pos_IT_Spe1, pos_IT_Spe2, pos_ITname_Spe1,pos_ITname_Spe2,species1,species2);
[n1n2 up_n1n2]= get_orthologs(neg_UpGeneName_Spe1,neg_UpGeneName_Spe2, neg_DwGeneName_Spe1, neg_DwGeneName_Spe2, neg_IT_Spe1, neg_IT_Spe2, neg_ITname_Spe1,neg_ITname_Spe2,species1,species2);
[p1n2 up_p1n2]= get_orthologs(pos_UpGeneName_Spe1,neg_UpGeneName_Spe2, pos_DwGeneName_Spe1, neg_DwGeneName_Spe2, pos_IT_Spe1, neg_IT_Spe2, pos_ITname_Spe1,neg_ITname_Spe2,species1,species2);
[n1p2 up_n1p2]= get_orthologs(neg_UpGeneName_Spe1,pos_UpGeneName_Spe2, neg_DwGeneName_Spe1, pos_DwGeneName_Spe2, neg_IT_Spe1, pos_IT_Spe2, neg_ITname_Spe1,pos_ITname_Spe2,species1,species2);
    
    
    AllCases=vertcat(p1p2,n1n2,p1n2,n1p2);
    Allups=vertcat(array2table(up_p1p2),array2table(up_n1n2),array2table(up_p1n2),array2table(up_n1p2));
 


 output_file='/home/nerdem/IT Project/IT_list_Salmonella_vs_Bacillus_conservative.xlsx';
 if(istable(AllCases))
 writetable(AllCases,output_file,'Sheet','All Co-Orthologs');
 else
     cout='No co-orthologs at all!'
 end
 if(istable(p1p2))
 writetable(p1p2,output_file,'Sheet','pos1_vs_pos2');
 else
     cout='No co-orthologs in the pos vs pos'
 end
 if(istable(n1n2))
 writetable(n1n2,output_file,'Sheet','neg1_vs_neg2');
 else
     cout='No co-orthologs in the neg vs neg'
 end
 
 if(istable(p1n2)) 
 writetable(p1n2,output_file,'Sheet','pos1_vs_neg2');
 else
     cout='No co-orthologs in the pos vs neg'
 end
 
 if(istable(n1p2))
 writetable(n1p2,output_file,'Sheet','neg1_vs_pos2');
 else
     cout='No co-orthologs in the neg vs pos'
 end
  if(istable(Allups))
 writetable(Allups,output_file,'Sheet','All_UpGeneOrthologs');
 else
     cout='No co-orthologs at all!'
 end

function [DU up]= get_orthologs(up_spe1,up_spe2, dw_spe1, dw_spe2, IT_spe1, IT_spe2, ITname_spe1,ITname_spe2, species1,species2)
[up, i1, i2] = intersect(up_spe1,up_spe2);
D=table('Size',size(up),'VariableTypes',"string" );

table_is_not_empty=0;
for i=1:size(up)
      
    if( isequal(subsref(dw_spe1(i1,:), struct('type', '()', 'subs', {{i, 1}})),subsref(dw_spe2(i2,:), struct('type', '()', 'subs', {{i, 1}}))))
      
        D(i,:)=subsref(dw_spe1(i1,:), struct('type', '()', 'subs', {{i, 1}})); %down
        U(i,:)=up(i,1); %up
        IT_1(i,:)= IT_spe1( i1(i,1),1); %IT of 1st species
        IT_2(i,:)= IT_spe2( i2(i,1),1); %IT of 1st species
        IT_names1(i,:)=ITname_spe1(i1(i,1),1);
        IT_names2(i,:)=ITname_spe2(i2(i,1),1);
        table_is_not_empty=1;

    end

end

if(table_is_not_empty)
    D=delete_empty_rows(D);
    U=delete_empty_rows(U);
    IT_1=delete_empty_rows(IT_1);
    IT_2=delete_empty_rows(IT_2);
    IT_names1=delete_empty_rows(IT_names1);
    IT_names2=delete_empty_rows(IT_names2);
    DU=[];
    rownames={species1; species2};
    rowlabels={};
    for n=1:size(U)
    
         DU_tmp=[IT_names1(n,1) U(n,1) IT_1(n,1) D(n,1); IT_names2(n,1) U(n,1) IT_2(n,1) D(n,1)];
         DU=[DU; DU_tmp];
         rowlabels=[rowlabels; rownames];
    end

        DU=[rowlabels, DU]; 
     
else
  DU=zeros(size(up)); 
    
end 


end 

 function y = delete_empty_rows(x)
 id_D=all(cellfun(@isempty,x{:,:}),2);
x(id_D,:)=[];
y=x;
end
 