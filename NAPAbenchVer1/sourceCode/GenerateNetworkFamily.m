function GenerateNetworkFamily(M,tree_file,Na,growth_model,out_path, grow_parameters)
% Network Family Geneating algorithm
%
%  This function generates a set of evolutionary related networks based on
%  the algorithm presented in the following paper:
%  - S.M.E. Sahraeian and B.J. Yoon, “A Network Synthesis Model for
%  Generating Protein Interaction Network Families”, PLoS ONE, 2012, to
%  appear.
%
%   Usage: 
%   GenerateNetworkFamily(M,tree_file,Na,growth_model,grow_parameters,out_path)
%
%   inputs  M               Number of networks
%           tree_file       The underlying phylogentic tree of the networks
%                           in Newick format
%           Na              Number of nodes in the ancestral network
%           growth_model    The duplication algorithm to be used to grow
%                           the networks. The options are 'DMC', 'DMR', and
%                           'CG' schemes.
%           out_path        Output directory where the generated networks
%                           files will be placed.
%           grow_parameters The parameters of the growth model: 
%                           for DMC: grow_parameters=[q_con, q_mod] 
%                           for DMR: grow_parameters=[q_new, q_del] 
%                           for CG: grow_parameters=delta
%   output files:
%           A.net,..        Network files. These files define the
%                           structure of each of generated networks.
%           A.fo,..       Functional annotation Network files. These
%                           files include the functional ortholgy group of
%                           each network node.
%           A-B.sim,..   The similarity scores of nodes across the
%                           networks.
%
%   Example:
%        GenerateNetworkFamily(5,'test\my_tree.txt',200,'DMC','test\out')
%
%   Input tree file format (Newick format):
%       For genrating the following tree:
%
%               /\
%           100/  \50
%             A   /\
%             150/  \110
%               /    \
%              B     /\
%                 60/  \70
%                  C   /\
%                   30/  \80
%                    D    E
%
%       The tree file should be written as:
%       (A:100,(B:150,(C:60,(D:30,E:80):70):110):50)
%       Then, if the Number of nodes in the ancestral network (Na) is 200,
%       the networks will be of sizes:
%       |A|=300, |B|=400, |C|=420, |D|=460, |E|=510
%
%   Onput file format:
%       Network files: [e.g A.net]
%               a1	a2
%               a3	a1
%               a4	a2
%               a2	a3
%       Functional annotation files: [e.g. A.fo]
%               a1	FO:1
%               a2	FO:2
%               a1	FO:2
%               a4	FO:3
%       Similarity score files: [e.g. A-B.sim]
%               a1	b1	153
%               a1	b3	55
%               a1	b7	49
%               a2	b3	444
%               a3	b3	211
%               a3	b4	122
%               a4	b5	251
%               a4	b8	71
%
%
% For more information on the algorithms, please see:
%
% S.M.E. Sahraeian and B.J. Yoon, “A Network Synthesis Model for Generating
% Protein Interaction Network Families”, PLoS ONE, 2012, to appear.
%
%  By Sayed Mohammad Ebrahim Sahraeian and Byung-Jun Yoon
%  July 2012
%  Contact: msahraeian@tamu.edu, bjyoon@ece.tamu.edu
%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Default parameters
%--------------------------------------------------------------------------
% DMC parameters
q_con=0.1;
q_mod=0.6;

% DMR parameters
q_new=0.12 ;
q_del=0.635;

% CG parameters
delta=4;

% Gamma Distributions parameters
ko=0.72;
to=226;
kn=0.85;
tn=73;

% functional inheritance and sequence similarity control parameters
lambda_max=.1;
p_fo=.9;
beta=1.6;

%Similarity score threshold
Ts=50;

if nargin==7
    if (strcmp(growth_model,'DMC'))
        q_con=grow_parameters(1);
        q_mod=grow_parameters(2);
    elseif (strcmp(growth_model,'DMR'))
        q_new=grow_parameters(1);
        q_del=grow_parameters(2);
    elseif (strcmp(growth_model,'CG'))
        delta=grow_parameters(1);
    end;
end;


tic
%Read tree file
[tree, branch_len,tree_dictionary]= Read_Tree(M,tree_file);


%Generate the initial network seed
[G0,F0]=Generate_seed(growth_model);

%Generate the ancestral network by growing the seed
[Ga,Fa]=Grow_seed(G0,F0,Na-size(G0,1),q_con,q_mod,q_new,q_del,delta,growth_model);

%Network Family sets
G_set=cell(1,2*M-1);
FO_set=cell(1,2*M-1);
S_set=cell(2*M-1,2*M-1);

%Traverse the phylogenetic tree to generate all the networks
internal_node=2*M-1;
processed_int_nodes=[];
G_set{2*M-1}=Ga;
FO_set{2*M-1}=Fa;
processed_nodes=2*M-1;
disp('Traverse the phylogenetic tree to generate all the networks')
for i=1:M-1    
    left_branch=min(tree(M-i,:));
    right_branch=max(tree(M-i,:));
    [G_set,S_set,FO_set,internal_node,processed_int_nodes,processed_nodes]=Traverse(G_set,S_set,FO_set,branch_len(left_branch),branch_len(right_branch),left_branch,right_branch,internal_node,processed_int_nodes,processed_nodes,q_con,q_mod,q_new,q_del,delta,growth_model,lambda_max,p_fo,to,ko,tn,kn,M,Ts,beta);
    disp(sprintf('%d/%d internal node(s) traversed',i,M-1))
end;
G_set=G_set(1:M);
S_set=S_set(1:M,1:M);
FO_set=FO_set(1:M);


disp('Writing output files');
if exist(out_path)~=7
    mkdir(out_path);
end;
%Write Output files
for i=1:M
    Write_Net(G_set{i},char(96+i),strcat(out_path,'\',char(64+i),'.net'));
    Write_FO(FO_set{i},char(96+i),strcat(out_path,'\',char(64+i),'.fo'));
end;    
for i=1:M
    for j=i+1:M
        Write_Sim(S_set{i,j},char(96+i),char(96+j),strcat(out_path,'\',char(64+i),'-',char(64+j),'.sim'));
    end;
end;
Write_log(tree_dictionary,G_set,M,out_path,growth_model,q_con,q_mod,q_new,q_del,delta)
t2=toc;
disp(sprintf('Total elapsed time: %f seconds',t2));
Disp_log(tree_dictionary,G_set,M,growth_model,q_con,q_mod,q_new,q_del,delta);
end
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%Read tree file
%--------------------------------------------------------------------------
function [tree,branch_len,tree_dictionary]=Read_Tree(M,tree_file)
disp('Reading the tree file')
tree=zeros(M-1,2);
branch_len=zeros(1,2*M-2);
f_tree=fopen(tree_file);
line=fgets(f_tree);

cnt_leaf=M;
cnt_internal=2*M-2;
cnt_tree=M-1;
queue={line};
tree_dictionary={};
while ~isempty(queue)
    curr_br=queue{1};
    queue=queue(2:end);
    tmp_br=curr_br(2:end);
    st=1;
    if (tmp_br(1)=='(')
        pr=1;
        st=2;
        while (pr~=0)
            f=strfind(tmp_br(st:end),'(');
            b=strfind(tmp_br(st:end),')');
            if ~isempty(f)
                f=f(1);
            else
                f=0;
            end;
            if ~isempty(b)
                b=b(1);
            else
                b=0;
            end;
            if (((f<b) && ~isempty(f)) || isempty(b))
                st=st+f;
                pr =pr+1;
            else
                st=st+b;
                pr =pr-1;
            end;
        end;
    end;
    c1=strfind(curr_br,'(');
    c1=c1(1);
    c2=strfind(curr_br(st:end),',');
    c2=c2(1);
    c3=strfind(curr_br,')');
    c3=c3(end);
    left=curr_br(c1+1:st+c2-2);
    right=curr_br(st+c2:c3-1);
    
    c4=strfind(left,':');
    c4=c4(end);
    c5=strfind(right,':');
    c5=c5(end);
    
    len_l=str2num(left(c4+1:end));
    len_r=str2num(right(c5+1:end));
    if left(1)=='('
        branch_len(cnt_internal)=len_l;
        tree(cnt_tree,1)=cnt_internal;
        cnt_internal=cnt_internal-1;
        queue=[queue, left(1:c4-1)];
    else
        branch_len(cnt_leaf)=len_l;
        tree(cnt_tree,1)=cnt_leaf;
        tree_dictionary{cnt_leaf}=left(1:c4-1);
        cnt_leaf=cnt_leaf-1;
    end;
    if right(1)=='('
        branch_len(cnt_internal)=len_r;
        tree(cnt_tree,2)=cnt_internal;
        cnt_internal=cnt_internal-1;
        queue=[queue, right(1:c5-1)];
    else
        branch_len(cnt_leaf)=len_r;
        tree(cnt_tree,2)=cnt_leaf;
        tree_dictionary{cnt_leaf}=right(1:c5-1);
        cnt_leaf=cnt_leaf-1;
    end;
    cnt_tree=cnt_tree-1;
end;
[cc,ii]=sort(tree_dictionary);
trans_dict=ii;
trans_dict=[trans_dict,M+1:2*M-2];
tree2=zeros(M-1,2);
branch_len2=zeros(1,2*M-2);
tree_dictionary2=cell(1,length(tree_dictionary));
for i=1:M
    tree_dictionary2{i}=tree_dictionary{trans_dict(i)};
end;

for i=1:M-1
    for j=1:2
        tree2(i,j)=find(trans_dict==tree(i,j));
    end;
end;


for i=1:(2*M-2)
    branch_len2(i)=branch_len(trans_dict(i));
end;

tree=tree2;
branch_len=branch_len2;
tree_dictionary=tree_dictionary2;
end



%--------------------------------------------------------------------------
%Generate the initial network seed
%--------------------------------------------------------------------------
function [G0,F0]=Generate_seed(growth_model)
disp('Generating the network seed')
if ~(strcmp(growth_model,'CG'))
    G0=zeros(17);
    G0(1:10,1:10)=1;
    G0(11:17,11:17)=1;
    for i=1:17
        G0(i,i)=0;
    end;
    for i=1:10
        for j=11:17
            G0(i,j)=rand(1)<0.1;
            G0(j,i)=G0(i,j);
        end;
    end;
    G_tmp=zeros(50);
    G_tmp(1:17,1:17)=G0;
    G0=G_tmp;
    for i=18:50
        if (rand(1)<.5)
            j=randperm(10);
            j=j(1);
            G0(j,i)=1;
            G0(i,j)=1;
        else
            j=randperm(7);
            j=j(1)+10;
            G0(j,i)=1;
            G0(i,j)=1;
        end;
    end;
    F0=[ones(1,10),ones(1,7)*2,3:35];
else
    G0=rand(4)<.1;
    G0=G0+G0';
    G0(G0>0)=1;
    F0=[1:size(G0,1)];
end;
G0=sparse(G0);
end
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
%Generate the ancestral network by growing the seed
%--------------------------------------------------------------------------
function [G_out,F_out]=Grow_seed(G_in,F_in,n,q_con,q_mod,q_new,q_del,delta,growth_model)
disp('Generate the ancestral network')
I=1;
G_out=G_in;
F_out=F_in;
if (strcmp(growth_model,'DMC'))
    for i=1:n
        r=randperm(length(G_out));
        r=r(1);
        E=G_out(r,:);
        f1=find(E==1);
        E0=E;
        rnd=(rand(1,length(f1))<.5);
        while (1==1)
            E=E0;
            E2=E;
            E(f1)=(rand(1,length(f1))<(1-q_mod)).*rnd+(1-rnd);
            E2(f1)=(rand(1,length(f1))<(1-q_mod)).*(1-rnd)+rnd;
            if ((sum(E)>0)&&(sum(E2)>0))
                break
            end
        end;
        G_out(r,:)=E;
        G_out(:,r)=E';
        G_out=[G_out,E2';E2,0];
        G_out(r,end)=rand(1)<q_con;
        G_out(end,r)=G_out(r,end);
        mm=max(F_out)+1;
        F_out=[F_out mm];
    end;
elseif (strcmp(growth_model,'DMR'))
    for i=1:n
        r=randperm(length(G_out));
        r=r(1);
        E=G_out(r,:);
        f0=find(E==0);
        f1=find(E==1);
        E0=E;
        while (1==1)
            E=E0;
            E(f0)=rand(length(f0),1)<(q_new/length(E));
            E(f1)=rand(length(f1),1)<(1-q_del);
            if (sum(E)>0)
                break
            end
        end;
        G_out=[G_out,E';E,0];
        mm=max(F_out)+1;
        F_out=[F_out mm];
    end;
elseif (strcmp(growth_model,'CG'))
    for i=1:n
        if sum(find(I==i))~=0
            C=Find_modules(G_out);
            I=ceil(size(G_out,1)*0.1)+i;
        end;
        pnew=1/length(C);
        nn=size(G_out,1);
        if (rand(1)<pnew)
            pd=nn+1-sum(G_out);
            pd=pd/sum(pd);
            Pd=cumsum(pd);
            t=sparse(1,size(G_out,1));
            while (sum(t)<delta)
                t(find(rand(1)<Pd,1))=1;
            end;
            C=[C {nn+1}];
        else
            t=sparse(1,size(G_out,1));
            f0=-1;
            while (sum(t)<delta)
                rnd=randperm(length(C));
                rnd=rnd(1);
                if (f0<0)
                    f0=rnd;
                end;
                G2=G_out(C{rnd},C{rnd});
                pd=size(G2,1)+1-sum(G2);
                pd=pd/sum(pd);
                Pd=cumsum(pd);
                anch=find(rand(1)<Pd,1);
                t(C{rnd}(anch))=1;
                nnn=find(G2(anch,:));
                if length(nnn)<=(delta-sum(t))
                    t(C{rnd}(nnn))=1;
                else
                    f2=randperm(length(nnn));
                    f2=f2(1:(delta-sum(t)));
                    t(C{rnd}(nnn(f2)))=1;
                end;
            end;
            C{f0}=[C{f0} nn+1];
        end;
        G_out=[G_out,t';t,0];
        mm=max(F_out)+1;
        F_out=[F_out mm];
    end;
end;
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%Traverse the phylogenetic tree to generate all the networks
%--------------------------------------------------------------------------
function [G_set,S_set,FO_set,internal_node,processed_int_nodes,processed_nodes]=Traverse(G_set,S_set,FO_set,d1,d2,left_branch,right_branch,internal_node,processed_int_nodes,processed_nodes,q_con,q_mod,q_new,q_del,delta,growth_model,lambda_max,p_fo,to,ko,tn,kn,M,Ts,beta);
pp=internal_node(1);
internal_node=internal_node(2:end);
processed_int_nodes=[processed_int_nodes,pp];
Gi=G_set{pp};
Fi=FO_set{pp};
G1=Gi;
G2=Gi;
F1=Fi;
F2=Fi;
for i=1:length(processed_nodes)
    if pp<processed_nodes(i)
        S=S_set{pp,processed_nodes(i)};
    else
        S=S_set{processed_nodes(i),pp}' ;
    end;
    if left_branch<processed_nodes(i)
        S_set{left_branch,processed_nodes(i)}=S;
    else
        S_set{processed_nodes(i),left_branch}=S';
    end;
    if right_branch<processed_nodes(i)
        S_set{right_branch,processed_nodes(i)}=S;
    else
        S_set{processed_nodes(i),right_branch}=S';
    end;
end;
S_set{left_branch,right_branch}=Generate_Similarities(G1,G2,to,ko,tn,kn,F1,F2,Ts,beta);
processed_nodes=[processed_nodes,left_branch,right_branch];
[G1,S_set,F1]=Functional_Grow(G1,S_set,setdiff(processed_nodes,processed_int_nodes),F1,d1,q_con,q_mod,q_new,q_del,delta,growth_model,left_branch,lambda_max,p_fo,kn,tn,Ts,beta);
[G2,S_set,F2]=Functional_Grow(G2,S_set,setdiff(processed_nodes,processed_int_nodes),F2,d2,q_con,q_mod,q_new,q_del,delta,growth_model,right_branch,lambda_max,p_fo,kn,tn,Ts,beta);
FO_set{left_branch}=F1;
FO_set{right_branch}=F2;
G_set{left_branch}=G1;
G_set{right_branch}=G2;
if (left_branch>M)
    internal_node=[internal_node,left_branch];
end;
if (right_branch>M)
    internal_node=[internal_node,right_branch];
end;
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Generate the similarity scores across the networks
%--------------------------------------------------------------------------
function S=Generate_Similarities(G,T,to,ko,tn,kn,F1,F2,Ts,beta);
N1=size(G,1);
N2=size(T,1);
rnd1=round(rand(N1,1).^(1/(1-beta)));
rnd2=round(rand(1,N2).^(1/(1-beta)));
S=sparse(N1,N2);
perm1=randperm(N1);
for i=1:N1
    perm2=randperm(N2);
    if (F1(perm1(i))>0)
        orth=find(F1(perm1(i))==F2(perm2));
    else
        orth=[];
    end;
    S(perm1(i),perm2(orth))=gamrnd(ko,to,length(orth),1)+Ts;
    perm2(orth)=[];
    n_nonorth=min(max(rnd1(perm1(i))-length(orth),0),length(perm2));
    perm2=perm2(sum(S(:,perm2)>0)<rnd2(perm2));
    perm2=perm2(randperm(length(perm2)));
    if (n_nonorth>length(perm2))
        n_nonorth=length(perm2);
    end;
    S(perm1(i),perm2(1:n_nonorth))=gamrnd(kn,tn,n_nonorth,1)+Ts;
end;
S=sparse(S);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Grow the netwoks considering the functional annotations
%--------------------------------------------------------------------------
function [G,S,F]=Functional_Grow(G,S,processed_nodes,F,d,q_con,q_mod,q_new,q_del,delta,growth_model,left_branch,lambda_max,p_fo,kn,tn,Ts,beta)
I=1;
processed_nodes=sort(processed_nodes);
S_set=S(processed_nodes,processed_nodes);
M=length(S_set);
ii=find(processed_nodes==left_branch);
if (strcmp(growth_model,'DMC'))
    for i=1:d
        r=randperm(length(G));
        r=r(1);
        E=G(r,:);
        f1=find(E==1);
        E0=E;
        rnd=(rand(1,length(f1))<.5);
        while (1==1)
            E=E0;
            E2=E;
            E(f1)=(rand(1,length(f1))<(1-q_mod)).*rnd+(1-rnd);
            E2(f1)=(rand(1,length(f1))<(1-q_mod)).*(1-rnd)+rnd;
            if ((sum(E)>0)&&(sum(E2)>0))
                break
            end
        end;
        G(r,:)=E;
        G(:,r)=E';
        G=[G,E2';E2,0];
        G(r,end)=rand(1)<q_con;
        G(end,r)=G(r,end);
        for j=1:ii-1
            snew=S_set{j,ii}(:,r).*(1-rand(size(S_set{j,ii},1),1)*lambda_max);
            snew=snew.*(snew>Ts*0.9);
            S_set{j,ii}=[S_set{j,ii},snew];
        end;
        for j=ii+1:M
            snew=S_set{ii,j}(r,:).*(1-rand(1,size(S_set{ii,j},2))*lambda_max);
            snew=snew.*(snew>Ts*0.9);
            S_set{ii,j}=[S_set{ii,j};snew];
        end;
        F=[F (rand(1)<p_fo)*F(r)];
        
    end;
elseif (strcmp(growth_model,'DMR'))
    for i=1:d
        rnd=randperm(length(G));
        rnd=rnd(1);
        E=G(rnd,:);
        f0=find(E==0);
        f1=find(E==1);
        E0=E;
        while (1==1)
            E=E0;
            E(f0)=rand(length(f0),1)<(q_new/length(E));
            E(f1)=rand(length(f1),1)<(1-q_del);
            if (sum(E)>0)
                break
            end
        end;
        G=[G,E';E,0];
        for j=1:ii-1
            snew=S_set{j,ii}(:,rnd).*(1-rand(size(S_set{j,ii},1),1)*lambda_max);
            snew=snew.*(snew>Ts*0.9);
            S_set{j,ii}=[S_set{j,ii},snew];
        end;
        for j=ii+1:M
            snew=S_set{ii,j}(rnd,:).*(1-rand(1,size(S_set{ii,j},2))*lambda_max);
            snew=snew.*(snew>Ts*0.9);
            S_set{ii,j}=[S_set{ii,j};snew];
        end;
        F=[F (rand(1)<p_fo)*F(rnd)];
    end;
elseif (strcmp(growth_model,'CG'))
    for i=1:d
        if sum(find(I==i))~=0
            C=Find_modules(G);
            I=ceil(size(G,1)*0.1)+i;
        end;
        n=size(G,1);
        pnew=1/length(C);
        if (rand(1)<pnew)
            pd=n+1-sum(G);
            pd=pd/sum(pd);
            Pd=cumsum(pd);
            t=sparse(1,size(G,1));
            while (sum(t)<delta)
                t(find(rand(1)<Pd,1))=1;
            end;
            Fnew=0;
            for j=1:ii-1
                nn=size(S_set{j,ii},1);
                ff=min(round(rand(1).^(1/(1-beta))),nn);
                gg=gamrnd(kn,tn,ff,1)+Ts;
                ff2=randperm(nn);
                ff2=ff2(1:ff);
                snew=sparse(nn,1);
                snew(ff2)=gg;
                S_set{j,ii}=[S_set{j,ii},snew];
            end;
            for j=ii+1:M
                nn=size(S_set{ii,j},2);
                ff=min(round(rand(1).^(1/(1-beta))),nn);
                gg=gamrnd(kn,tn,1,ff)+Ts;
                ff2=randperm(nn);
                ff2=ff2(1:ff);
                snew=sparse(1,nn);
                snew(ff2)=gg;
                S_set{ii,j}=[S_set{ii,j};snew];
            end;
            C=[C {n+1}] ;
        else
            t=sparse(1,size(G,1));
            anch_flag=-1;
            f0=-1;
            while (sum(t)<delta)
                rnd=randperm(length(C));
                rnd=rnd(1);
                G2=G(C{rnd},C{rnd});
                pd=size(G2,1)+1-sum(G2);
                pd=pd/sum(pd);
                Pd=cumsum(pd);
                anch=find(rand(1)<Pd,1);
                
                if (anch_flag<0)
                    anch_flag=C{rnd}(anch);
                    Fnew=(rand(1)<p_fo)*F(anch_flag);
                    f0=rnd;
                end;
                t(C{rnd}(anch))=1;
                nn=find(G2(anch,:));
                if length(nn)<=(delta-sum(t))
                    t(C{rnd}(nn))=1;
                else
                    f2=randperm(length(nn));
                    f2=f2(1:(delta-sum(t)));
                    t(C{rnd}(nn(f2)))=1;
                end;
            end;
            for j=1:ii-1
                snew=S_set{j,ii}(:,anch_flag).*(1-rand(size(S_set{j,ii},1),1)*lambda_max);
                snew=snew.*(snew>Ts*0.9);
                S_set{j,ii}=[S_set{j,ii},snew];
            end;
            for j=ii+1:M
                snew=S_set{ii,j}(anch_flag,:).*(1-rand(1,size(S_set{ii,j},2))*lambda_max);
                snew=snew.*(snew>Ts*0.9);
                S_set{ii,j}=[S_set{ii,j};snew];
            end;
            C{f0}=[C{f0} n+1]  ;
        end;
        G=[G,t';t,0];
        F=[F Fnew];
    end;
end;
S(processed_nodes,processed_nodes)=S_set;
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Find the network modules
%--------------------------------------------------------------------------
function C=Find_modules(G)
N=size(G,1);
C=cell(1,N);
for i=1:N
    C{i}=i;
end;
M=nnz(G);
E=(G/2+diag(diag(G)/2))/M;
A=sum(E);
Q=sum(diag(E)'-A.^2);
while (1==1)
    [a,b,j]=find(triu(E,1));
    if isempty(a)
        break
    end;
    DQ=[2*j-A(a)'.*A(b)'];
    [mi,ii]=max(DQ);
    if mi<0
        break
    end;
    Q=Q+DQ(ii);
    C{a(ii)}=([C{a(ii)},C{b(ii)}]);
    C(b(ii):end-1)=C(b(ii)+1:end);
    C=C(1:end-1);
    E(a(ii),:)=E(a(ii),:)+E(b(ii),:);
    E(b(ii),:)=[];
    E(:,a(ii))=E(:,a(ii))+E(:,b(ii));
    E(:,b(ii))=[];
    A(a(ii))=A(a(ii))+A(b(ii));
    A(b(ii))=[];
    
end;
end
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
% Write the Network output file
%--------------------------------------------------------------------------
function Write_Net(G,a,my_file)
fid = fopen(my_file,'w');
[i,j,c]=find(G);
for k=1:length(i)
    if (i(k)<=j(k))
        fprintf(fid, '%s%d\t%s%d\n', a,i(k),a,j(k));
    end;
end;
fclose(fid);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Write the Network Functional annotation
%--------------------------------------------------------------------------
function Write_FO(F,a,my_file)
fid = fopen(my_file,'w');

for i=1:length(F)
    if F(i)>0
        fprintf(fid, '%s%d\tFO%d\n', a,i,F(i));
    end;
end;
fclose(fid);
end
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Write the Networks similarity scores
%--------------------------------------------------------------------------
function Write_Sim(S,a,b,my_file)
fid = fopen(my_file,'w');

[i,j,c]=find(S);
for k=1:length(i)
    fprintf(fid, '%s%d\t%s%d\t%f\n', a,i(k),b,j(k),c(k));
end;
fclose(fid);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Write the log file
%--------------------------------------------------------------------------
function Write_log(tree_dictionary,G_set,M,out_path,growth_model,q_con,q_mod,q_new,q_del,delta)
	fid=fopen(strcat(out_path,'\log_file.txt'),'w');
	fprintf(fid,'Generated Network Family Information\n');
	fprintf(fid,'-------------------------------------\n');
	fprintf(fid,'Number of Generated Networks:\t%d\n',M);
	fprintf(fid,'Growth Model:\t%s\n',growth_model);
	fprintf(fid,'Growth Parameters:\t');
	if strcmp(growth_model,'DMC')
		fprintf(fid,'q_con=%f\t q_mod=%f\t\n',q_con,q_mod);
    elseif strcmp(growth_model,'DMR')
		fprintf(fid,'q_new=%f\t q_del=%f\t\n',q_new,q_del);
    elseif strcmp(growth_model,'CG')
		fprintf(fid,'delta=%d\t\n',delta);
    end;
	fprintf(fid,'\nList of Generated Networks:\n');
	fprintf(fid,'Name\tLetter\tSize\n');
	for i=1:M
		fprintf(fid,'%s\t\t%s\t\t%d\n',tree_dictionary{i},char(64+i),size(G_set{i},1));
    end;
	fclose(fid);
end
   
%--------------------------------------------------------------------------
%Display output log	
%--------------------------------------------------------------------------
function Disp_log(tree_dictionary,G_set,M,growth_model,q_con,q_mod,q_new,q_del,delta)
	disp('-------------------------------------');
	disp('Generated Network Family Information');
	disp('-------------------------------------');
	disp(sprintf('Number of Generated Networks:\t%d',M));
	disp(sprintf('Growth Model:\t%s',growth_model));
	s='Growth Parameters:';
	if strcmp(growth_model,'DMC')
		s=sprintf('%s\tq_con=%f\t q_mod=%f\t',s,q_con,q_mod);
    elseif strcmp(growth_model,'DMR')
		s=sprintf('%s,q_new=%f\t q_del=%f\t',s,q_new,q_del);
    elseif strcmp(growth_model,'CG')
		s=sprintf('%sdelta=%d\t',s,delta);
    end;
    disp(s);
	disp('List of Generated Networks:');
	disp(sprintf('Name\tLetter\tSize'));
	for i=1:M
		disp(sprintf('%s\t\t%s\t\t%d',tree_dictionary{i},char(64+i),size(G_set{i},1)));
    end;
end
