function [ Cluster_D ] = Cluster( sim_mat,name_list)
% demo: [ ClustermiRNAs ] = Cluster( MM_mat,miRNAs_list(:,2)); 

dn = size(sim_mat,1);
tnSim = size(dn,dn);

Sim_mat2= sim_mat-diag(diag(sim_mat)); 
Sim_mat2(Sim_mat2<0.3)=0;
for i = 1:dn-1
    for j = 2:dn
        com = Sim_mat2(i,:)&Sim_mat2(j,:);
        tnSim(i,j) = sum(com);
        tnSim(j,i) = sum(com);
    end
end
for i = 1:dn
   tnSim(i,i) = nnz(Sim_mat2(i,:));
end

dfid = fopen('Net.txt','wt');
for i = 2:dn
    for j = 1:i-1
        numShared = tnSim(i,j);
        if(numShared>0)
            fname =char(name_list(i));
            fprintf(dfid,'%s',fname);
            fprintf(dfid,'\t');
            sname = char(name_list(j));
            fprintf(dfid,'%s',sname);
            fprintf(dfid,'\t');
            fprintf(dfid,'%d',numShared);
            fprintf(dfid,'\n'); 
        end  
    end
end
fclose(dfid);
clc;

cdfid = fopen('ClusterResult.txt','w');
fclose(cdfid);
diary('ClusterResult.txt');
diary on
!java -jar "cluster_one-1.0.jar"  "Net.txt" -F csv
diary off

% %extract cluster results from the ClusterResult.txt;
ClusterDiseases = cell(100,1);
ClusterDiseasesQuality = size(100,1);
cdfid = fopen('ClusterResult.txt','r');
flag = 0;
diseasecluster_n = 0;
while(true)
    tline = fgetl(cdfid);
    if(tline==-1)
        break;
    end
    if(flag==1)
        while(true)
            tline = fgetl(cdfid);
            if(tline==-1)
                break;
            end
            linestr = regexp(tline,',','split');
            quality = linestr(6);
            pvalue = linestr(7);
            quality = str2double(quality);
            pvalue = str2double(pvalue);
            
            if(pvalue<0.001)
                m_clusterstr =char(linestr(8));
                m_clusterstr = strtrim(m_clusterstr);
                m_clusterstr(1) = [];
                cstrlen = size(m_clusterstr,2);
                m_clusterstr(cstrlen) = [];
                
                diseasecluster_n = diseasecluster_n+1;
                ClusterDiseases{diseasecluster_n} = m_clusterstr;
                ClusterDiseasesQuality(diseasecluster_n) = quality;
            end  
                Cluster_D=cell(diseasecluster_n,1);
                Cluster_D(1:diseasecluster_n,1) = ClusterDiseases(1:diseasecluster_n,1);                
        end
    end
    if(flag==1)
        break;
    end
    if(size(strfind(tline,'Detected'),1)~=0)
        flag = 1;
    end
end
fclose(cdfid);

end

