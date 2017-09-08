function [ graph ] = Graph( network , name_list, p )
% demo: [ d_graph_mat ] = Graph( DD_mat , diseases_list,p );
%       [ m_graph_mat ] = Graph( MM_mat , miRNAs_list,p );
tic    
[rows, cols] = size( network );    
          
    %----------------------clustering--------------------%
     Clusters = Cluster( network,name_list(:,2));
     num_Clu=size(Clusters,1);
     Clu_mat = zeros(rows, cols);
     for i = 1 : num_Clu
             linestr = regexp(Clusters{i,1},' ','split');
             num_id = length(linestr);
             for j = 1 : num_id-1 
                 [~,idx01]=intersect( name_list(:,2),linestr(1,j));
                 for k = j+1 : num_id
                    [~,idx02]=intersect( name_list(:,2),linestr(1,k));
                    Clu_mat(idx01,idx02) = 1;
                    Clu_mat(idx02,idx01) = 1;
                 end                 
             end         
     end    
        
    %------------------p nearest neighbors--------------%
    network= network-diag(diag(network)); 
    PNN = zeros(rows, cols);
    graph = zeros(rows, cols);
    [sort_network,idx]=sort(network,2,'descend');
    for i = 1 : rows
        PNN(i,idx(i,1:p))=sort_network(i,1:p);
    end    
    for i = 1 : rows
        idx_i=find(PNN(i,:));
        for j = 1 : rows            
            idx_j=find(PNN(j,:));
            if ismember(j,idx_i) && ismember(i,idx_j) && isequal(Clu_mat(i,j),1)               
                graph(i,j)=1;
            elseif ~ismember(j,idx_i) && ~ismember(i,idx_j) && ~isequal(Clu_mat(i,j),1)  
                graph(i,j)=0;
            else
                graph(i,j)=0.5;               
            end       
       end
    end
    
 toc

end



