function [A,B]=GRNMF(inputdata,params)
% Xiao, Q.,Luo J.W. et al, (2017),A Graph Regularized Non-negative Matrix Factorization  
% Method for Identifying MicroRNA-disease Associations, Bioinformatics.
% College of Computer Science and Electronics Engineering,Hunan University
% hnuyldf@hnu.edu.cn
% 2017/09/07

K=params.K; 
r=params.r;
p=params.p;
MD_mat=inputdata.MD_mat;
MM_mat=inputdata.MM_mat;
DD_mat=inputdata.DD_mat;
miRNAs_list=inputdata.miRNAs_list;
diseases_list=inputdata.diseases_list;

Y=WKNKN(MD_mat, MM_mat, DD_mat,K,r); 
m_graph_mat = Graph( MM_mat ,miRNAs_list, p ); 
d_graph_mat = Graph( DD_mat ,diseases_list, p );
m_mat_new = m_graph_mat.* MM_mat ;  
d_mat_new = d_graph_mat.* DD_mat ;
clear K r p;
 
k=params.k;
iterate=params.iterate; 
lamda=params.lamda;   
lamda_m=params.lamda_m;
lamda_d=params.lamda_d; 
fprintf('k=%d  maxiter=%d  lamda=%d  lamda_m=%d lamda_d=%d\n', k, iterate,lamda, lamda_m,lamda_d); 

[rows,cols] = size(Y);
A=abs(rand(rows,k));        
B=abs(rand(cols,k));

diag_m = diag(sum(m_mat_new,2));
diag_d = diag(sum(d_mat_new,2));
L_m =diag_m -m_mat_new; 
L_d =diag_d -d_mat_new;


fid = fopen( 'RunResult.txt','wt+');
for step=1:iterate
        YB = Y*B;
        BB =  B'*B;
        ABB = A*BB;        
        if lamda > 0 && lamda_m >0
            SA = m_mat_new*A;
            DA = diag_m*A;            
            YB = YB + lamda_m*SA;
            ABB = ABB + lamda*A + lamda_m*DA;
        end
        A = A.*(YB./ABB);
        YA = Y'*A;
        AA = A'*A;
        BAA = B*AA;
        if lamda > 0 && lamda_d >0
            SB = d_mat_new*B;
            DB = diag_d*B;
            YA = YA + lamda_d*SB;
            BAA = BAA + lamda*B + lamda_d*DB;
        end
        B = B.*(YA./BAA);
        
        dY = Y-A*B';
        obj_NMF = sum(sum(dY.^2));
        ABF = sum(sum(A.^2))+sum(sum(B.^2));
        ALA = sum(sum((A'*L_m).*A'));
        BLB = sum(sum((B'*L_d).*B'));
        obj = obj_NMF+lamda*ABF+lamda_m*ALA+lamda_d*BLB;
        error = mean(mean(abs(Y-A*B')))/mean(mean(Y));      
        fprintf(fid,'%s\n',[sprintf('step = \t'),int2str(step),...
            sprintf('\t obj = \t'),num2str(obj),...
		    sprintf('\t error = \t'), num2str(error)]);
        fprintf('step=%d  obj=%d  error=%d\n',step, obj, error);   
        if error< 10^(-5)
            fprintf('step=%d\n',step);
            break;
        end            
       
end
fclose(fid);

end

 
