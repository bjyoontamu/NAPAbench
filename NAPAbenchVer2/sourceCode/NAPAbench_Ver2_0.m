function varargout = NAPAbench_Ver2_0(varargin)
% NAPABENCH_VER2_0 MATLAB code for NAPAbench_Ver2_0.fig
%      NAPABENCH_VER2_0, by itself, creates a new NAPABENCH_VER2_0 or raises the existing
%      singleton*.
%
%      H = NAPABENCH_VER2_0 returns the handle to a new NAPABENCH_VER2_0 or the handle to
%      the existing singleton*.
%
%      NAPABENCH_VER2_0('CALLBACK',hObject,eventData,handles,...) calls the
%      local
%      function named CALLBACK in NAPABENCH_VER2_0.M with the given input arguments.
%
%      NAPABENCH_VER2_0('Property','Value',...) creates a new NAPABENCH_VER2_0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NAPAbench_Ver2_0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NAPAbench_Ver2_0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NAPAbench_Ver2_0

% Last Modified by GUIDE v2.5 08-Jan-2019 09:54:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NAPAbench_Ver2_0_OpeningFcn, ...
                   'gui_OutputFcn',  @NAPAbench_Ver2_0_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NAPAbench_Ver2_0 is made visible.
function NAPAbench_Ver2_0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NAPAbench_Ver2_0 (see VARARGIN)

% Choose default command line output for NAPAbench_Ver2_0
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Define default parameter
params.treeFileName = {'2way.txt', '5way.txt', '8way.txt'};
params.SS = [2000, 1000, 700];
params.IP = [0.15, 0.4; 0.2, 0.45; 4, 0; 0.75, 50];
params.IPLabel = {{'q_con', 'q_del'}, {'q_new', 'q_del'}, {'delta'}, {'s_del', 's_f'}};
setappdata(hObject, 'parameters', params);

defaultAll(hObject, handles)


% UIWAIT makes NAPAbench_Ver2_0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NAPAbench_Ver2_0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in generate.
function generate_Callback(hObject, eventdata, handles)
% hObject    handle to generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Enable', 'Off');
drawnow;
numberOfNetworks = str2num(get(handles.eTheNumberOfNetworks, 'String'));
numberOfFamilies = str2num(get(handles.eTheNumberOfDataset, 'String'));
contents = cellstr(get(handles.pNetworkGrowthModel,'String')); 
networkGrowthModel = contents{get(handles.pNetworkGrowthModel,'Value')};
sizeOfSeedNetwork = str2num(get(handles.eSeedNetworkSize, 'String'));
intraParams(1) = str2num(get(handles.eIP1, 'String'));
intraParams(2) = str2num(get(handles.eIP2, 'String'));
treeFileName = get(handles.bPyloGenTreeFilePath, 'String');

if mod(numberOfNetworks, 1) ~= 0 || numberOfNetworks < 2
    message = sprintf('Error: The number of networks should be a positive integer larger than 1');
    uiwait(warndlg(message));
    set(hObject, 'Enable', 'On');
    drawnow;
    return
end

if mod(numberOfFamilies, 1) ~= 0 || numberOfFamilies < 1
    message = sprintf('Error: The number of dataset should be a positive integer larger than 0');
    uiwait(warndlg(message));
    set(hObject, 'Enable', 'On');
    drawnow;
    return
end

if mod(sizeOfSeedNetwork, 1) ~= 0 || sizeOfSeedNetwork < 1
    message = sprintf('Error: The number of dataset should be a positive integer larger than 0');
    uiwait(warndlg(message));
    set(hObject, 'Enable', 'On');
    drawnow;
    return
end

for f = 1: numberOfFamilies
    pathForFamily = strcat('family', num2str(f), '/');
    if (numberOfFamilies == 1)
        pathForFamily = '';
    end
    outputPath = strcat('output/', networkGrowthModel, '/', pathForFamily);
    GenerateNetworkFamily(f, numberOfNetworks, treeFileName, sizeOfSeedNetwork, networkGrowthModel, outputPath, intraParams);
end
set(hObject, 'Enable', 'On');
drawnow;


function GenerateNetworkFamily(iteration, M, tree_file, Na, growth_model, out_path, grow_parameters)
%--------------------------------------------------------------------------
%Default parameters
%--------------------------------------------------------------------------
% DMC parameters
q_con = 0.15;
q_mod = 0.4;

% DMR parameters
q_new= 0.2 ;
q_del = 0.45;

% CG parameters
delta=4;

% STICKY parameters
s_del = 0.75;
s_f = 50;

% Gamma Distributions parameters
ko=0.94;
to=169.49;
kn=0.86;
tn=42.00;

% functional inheritance and sequence similarity control parameters
lambda_max=.1;
p_fo=.9;
beta = 1.7;

%Similarity score threshold
Ts = 45;


if (strcmp(growth_model,'DMC'))
    q_con=grow_parameters(1);
    q_mod=grow_parameters(2);
elseif (strcmp(growth_model,'DMR'))
    q_new=grow_parameters(1);
    q_del=grow_parameters(2);
elseif (strcmp(growth_model,'CG'))
    delta=grow_parameters(1);
elseif (strcmp(growth_model,'STICKY'))
    s_del = grow_parameters(1);
    s_f = grow_parameters(2);
end


tic
%Read tree file
[tree, nodeSizes] = Read_Tree(tree_file, M, Na);
traverseOrder = bfsearch(tree, 'S');
[~, nodeSizeIndices] = ismember(traverseOrder, nodeSizes{1});
nodeSizeInTraverseOrder = [Na; nodeSizes{2}(nodeSizeIndices(2:end))];
assert(sum(sort(nodeSizeInTraverseOrder) ~= nodeSizeInTraverseOrder)==0);
treeDepth = 0;
treeLevels = [0];
for i = 2: length(traverseOrder)
   treeDepth = max(treeDepth, length(shortestpath(tree, 'S', traverseOrder{i}))-1);
   treeLevels = [treeLevels length(shortestpath(tree, 'S', traverseOrder{i}))-1];
end
outNodeIndices = [];
for i = 1: length(traverseOrder)
    if isempty(successors(tree, traverseOrder{i}))
        outNodeIndices = [outNodeIndices i];
    end
end
G_set = cell(1, length(traverseOrder));
FO_set = cell(1, length(traverseOrder));
S_set = cell(length(traverseOrder), length(traverseOrder));
if iteration == 1
    figure; 
    plot(tree);
end
%Generate the initial network seed
[G0, F0, seeds] = Generate_seed(growth_model);

%Generate the ancestral network by growing the seed
[Ga, Fa, seeds] = Grow_seed(G0, F0, Na-size(G0,1), q_con, q_mod, q_new, q_del, delta, s_del, s_f, growth_model, seeds);
G_set{1}=Ga; 
FO_set{1}=Fa;

disp('Traverse the phylogenetic tree to generate all the networks')
[tree, G_set, S_set, FO_set] = Traverse(tree, traverseOrder, nodeSizeInTraverseOrder, treeDepth, treeLevels, ...
                                            G_set, S_set, FO_set, ...
                                            q_con, q_mod, ...   % Parmeter for DMC
                                            q_new, q_del, ...   % Parmeter for DMR
                                            delta, ...          % Parmeter for CG
                                            s_del, s_f, ...     % Parmeter for STICKY
                                            growth_model, ...
                                            lambda_max, p_fo, to, ko, tn, kn, M, Ts, beta);
G_set = G_set(outNodeIndices);
S_set = S_set(outNodeIndices, outNodeIndices);
FO_set = FO_set(outNodeIndices);

disp('Writing output files');
if exist(out_path)~=7
    mkdir(out_path);
end

% Write Output files
for i=1:M
    Write_Net(G_set{i}, lower(traverseOrder{outNodeIndices(i)}),strcat(out_path, upper(traverseOrder{outNodeIndices(i)}),'.net'));
    Write_FO(FO_set{i}, lower(traverseOrder{outNodeIndices(i)}),strcat(out_path, upper(traverseOrder{outNodeIndices(i)}),'.fo'));
end
    
for i=1:M
    for j=i+1:M
        Write_Sim(S_set{i,j}, ...
            lower(traverseOrder{outNodeIndices(i)}), lower(traverseOrder{outNodeIndices(j)}), ...
            strcat(out_path, upper(traverseOrder{outNodeIndices(i)}), '-', upper(traverseOrder{outNodeIndices(j)}), '.sim'));
    end
end
Write_log(traverseOrder(outNodeIndices),G_set,M,out_path,growth_model,q_con,q_mod,q_new,q_del,delta,s_del,s_f, Na)
t2=toc;                                          
disp(sprintf('Total elapsed time: %f seconds',t2));
Disp_log(traverseOrder(outNodeIndices),G_set,M,growth_model,q_con,q_mod,q_new,q_del,delta,s_del,s_f, Na);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%Read tree file
%--------------------------------------------------------------------------
function [tree, nodeSizes]=Read_Tree(tree_file, M, Na)
disp('Reading the tree file')

fTree = fopen(tree_file);
treeInMarkov = textscan(fTree, '%c%*2s%c:%d', 1000000, 'commentStyle', '//');
fclose(fTree);

tree = digraph(cellstr(treeInMarkov{1})', cellstr(treeInMarkov{2})');
nodeSizes = {cellstr(treeInMarkov{2}),treeInMarkov{3}};
assert(isdag(tree));
assert(sum(Na > nodeSizes{2})==0);


%--------------------------------------------------------------------------
%Generate the initial network seed
%--------------------------------------------------------------------------
function [G0,F0,seeds]=Generate_seed(growth_model)
disp('Generating the network seed')
if strcmp(growth_model,'DMC')
    growth_type = 1;    % duplication type
elseif strcmp(growth_model,'DMR')
    growth_type = 1;
elseif strcmp(growth_model,'CG')
    growth_type = 0;
elseif strcmp(growth_model,'STICKY')
    growth_type = 0;
end

if growth_type == 1
    seeds = [];
    
    % if ~(strcmp(growth_model,'CG'))
    G0=zeros(17);
    G0(1:10,1:10)=1;
    G0(11:17,11:17)=1;
    for i=1:17
        G0(i,i)=0;
    end
    % randomly connect two cliques
    for i=1:10
        for j=11:17
            G0(i,j)=rand(1)<0.1;
            G0(j,i)=G0(i,j);
        end
    end
    G_tmp=zeros(50);
    G_tmp(1:17,1:17)=G0;
    G0=G_tmp;
    % randomly connect new nodes to one of clique
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
        end
    end
    % functional annotation
    F0=[ones(1,10),ones(1,7)*2,3:35];
elseif growth_type == 0  % seed network for CG model
    seeds = [];
    
    % need to check the network is empty or not
    G0 = zeros(4,4);
    while(length(find(sum(G0))) ~= 4)
        G0=rand(4)<.1;
        G0=G0+G0';
        G0(G0>0)=1;
        F0=[1:size(G0,1)];
    end
end
G0=sparse(G0);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%Generate the ancestral network by growing the seed
%--------------------------------------------------------------------------
function [G_out,F_out,seeds] = Grow_seed(G_in, F_in, n,q_con, q_mod, q_new, q_del, delta, s_del, s_f, growth_model,seeds)
disp('Generate the ancestral network')
I=1;
G_out=G_in;
F_out=F_in;
if (strcmp(growth_model,'DMC'))
    seeds = [];
    for i=1:n
        % randomly select one of seed node
        r=randperm(length(G_out));
        r=r(1);
        % neighbor of the selected seed node
        E=G_out(r,:);
        f1=find(E==1);
        E0=E;
        % randomly select neighbor nodes
        rnd=(rand(1,length(f1))<.5);
        while (1==1)
            E=E0;
            E2=E;
            % random perturbation between u-v
            E(f1)=(rand(1,length(f1))<(1-q_mod)).*rnd+(1-rnd);
            % random link between u-v'
            E2(f1)=(rand(1,length(f1))<(1-q_mod)).*(1-rnd)+rnd;
            % modified model (to remove singleton)
            if ((sum(E)>0)&&(sum(E2)>0))
                break
            end
        end
        G_out(r,:)=E;
        G_out(:,r)=E';
        G_out=[G_out,E2';E2,0];
        G_out(r,end)=rand(1)<q_con;
        G_out(end,r)=G_out(r,end);
        mm=max(F_out)+1;
        % assign new functional annotation for the new node
        F_out=[F_out mm];
    end
elseif (strcmp(growth_model,'DMR'))
    seeds = [];
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
        end
        G_out=[G_out,E';E,0];
        mm=max(F_out)+1;
        F_out=[F_out mm];
    end
elseif (strcmp(growth_model,'CG'))
    seeds = [];
    for i=1:n
        if sum(find(I==i))~=0
            C=Find_modules(G_out);
            I=ceil(size(G_out,1)*0.1)+i;
        end
        pnew=1/length(C);
        nn=size(G_out,1);
        if (rand(1)<pnew)
            pd=nn+1-sum(G_out);
            pd=pd/sum(pd);
            Pd=cumsum(pd);
            t=sparse(1,size(G_out,1));
            while (sum(t)<delta)
                t(find(rand(1)<Pd,1))=1;
            end
            C=[C {nn+1}];
        else
            t=sparse(1,size(G_out,1));
            f0=-1;
            while (sum(t)<delta)
                rnd=randperm(length(C));
                rnd=rnd(1);
                if (f0<0)
                    f0=rnd;
                end
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
                end
            end
            C{f0}=[C{f0} nn+1];
        end
        G_out=[G_out,t';t,0];
        mm=max(F_out)+1;
        F_out=[F_out mm];
    end
elseif (strcmp(growth_model,'STICKY'))
    seeds = [];
    beta = 1.6;
    seed_size = size(G_out,1);
    
    % assign sticky index
    for i=1:n
        x = 1:(seed_size + i - 1);
        y = x.^(-beta);
        y = y./sum(y);
        y_cdf = cumsum(y);
        deg_new = find(histcounts(rand, full([0, y_cdf])));
        
        deg_dist = [sum(G_out,2); deg_new];
        sticky_ind = deg_dist./(sqrt(sum(deg_dist)));
        sticky_val = sticky_ind(end);
        
        while 1
            % random link between u-v'
            E = rand(length(sticky_ind)-1,1) < (sticky_ind(1:end-1)*sticky_val);
            nbr_new = find(E);
            if length(nbr_new) > 1
                E(nbr_new) = rand(length(nbr_new), 1) > s_del;
                nbr_new = find(E);
            end
            if ~isempty(nbr_new)
                break;
            end
        end
        G_out=[G_out, E; E',0];
        
        mm = max(F_out) + 1;
        % assign new functional annotation for the new node
        F_out=[F_out mm];
    end
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% the phylogenetic tree to generate all the networks
%--------------------------------------------------------------------------
function [tree, G_set, S_set, FO_set] = Traverse(tree, traverseOrder, nodeSizeInTraverseOrder, treeDepth, treeLevels, ...
                                                G_set, S_set, FO_set, ...
                                                q_con, q_mod, ...   % Parmeter for DMC
                                                q_new, q_del, ...   % Parmeter for DMR
                                                delta, ...          % Parmeter for CG
                                                s_del, s_f, ...     % Parmeter for STICKY
                                                growth_model, ...
                                                lambda_max, p_fo, to, ko, tn, kn, M, Ts, beta)
for l = 1: treeDepth
    % Find node whose depth is equal to l
    nodeIndices = find(treeLevels == l);
    assert(~mod(length(nodeIndices), 2));
    
    while ~isempty(nodeIndices)
        lIndex = nodeIndices(1);
        rIndex = nodeIndices(2);
        nodeIndices(1:2) = [];
        fprintf('Traversing nodes %s and %s\n', traverseOrder{lIndex}, traverseOrder{rIndex})
        
        % Error check l node and r node should have the same parent
        [~, parentIndexL] = ismember(predecessors(tree, traverseOrder{lIndex}), traverseOrder);
        [~, parentIndexR] = ismember(predecessors(tree, traverseOrder{rIndex}), traverseOrder);
        assert(parentIndexL == parentIndexR)
        
        Gp = G_set{parentIndexL};
        Fp = FO_set{parentIndexL};
        
        Gl = Gp; Fl = Fp;
        Gr = Gp; Fr = Fp;

        % Copy sim. score
        for j = 1: parentIndexL - 1
            S_set{j, lIndex} = S_set{j, parentIndexL};
            S_set{j, rIndex} = S_set{j, parentIndexL};
        end
        
        % Copy sim. score
        for j = parentIndexL + 1: lIndex - 1
            S_set{j, lIndex} = S_set{parentIndexL, j}';
            S_set{j, rIndex} = S_set{parentIndexL, j}';
        end
        
        S_set{lIndex, rIndex} = Generate_Similarities(Gl, Gr, to, ko, tn, kn, Fl, Fr, Ts ,beta);
        
        [Gl, S_set, Fl] = functionalGrowth(lIndex, Gl, S_set, Fl, ...
            double(nodeSizeInTraverseOrder(lIndex) - nodeSizeInTraverseOrder(parentIndexL)), ...
            q_con, q_mod, q_new, q_del, delta, s_del, s_f, ...
            growth_model, lambda_max, p_fo, kn, tn, Ts, beta);
        [Gr, S_set, Fr] = functionalGrowth(rIndex, Gr, S_set, Fr, ...
            double(nodeSizeInTraverseOrder(rIndex) - nodeSizeInTraverseOrder(parentIndexL)), ...
            q_con, q_mod, q_new, q_del, delta, s_del, s_f, ...
            growth_model, lambda_max, p_fo, kn, tn, Ts, beta);
        
        FO_set{lIndex} = Fl;
        G_set{lIndex} = Gl;
        FO_set{rIndex} = Fr;
        G_set{rIndex} = Gr;
    end
end


function S = Generate_Similarities(G, T, to, ko, tn, kn, F1, F2, Ts, beta)
N1=size(G,1);
N2=size(T,1);

x = 1:N2;
y = x.^(-beta);
y = y./sum(y);
y_cdf = cumsum(y);
[~, ~, rnd1] = histcounts(rand(N1,1), full([0, y_cdf]));

x = 1:N1;
y = x.^(-beta);
y = y./sum(y);
y_cdf = cumsum(y);
[~, ~, rnd2] = histcounts(rand(1,N2), full([0, y_cdf]));

S=sparse(N1,N2);
perm1=randperm(N1);
for i=1: N1
    perm2=randperm(N2);
    if (F1(perm1(i))>0)
        orth=find(F1(perm1(i))==F2(perm2)); % find orthologous proteins
    else
        orth=[];
    end
    S(perm1(i),perm2(orth))=gamrnd(ko,to,length(orth),1)+Ts;
    perm2(orth)=[];
    n_nonorth=min(  max(  rnd1(perm1(i))-length(orth),  0  ), length(perm2)  );
    perm2=perm2(sum(S(:,perm2)>0)<rnd2(perm2));
    perm2=perm2(randperm(length(perm2)));
    if (n_nonorth>length(perm2))
        n_nonorth=length(perm2);
    end
    S(perm1(i),perm2(1:round(n_nonorth/4))) = gamrnd(kn,tn,round(n_nonorth/4),1)+Ts;
end
S=sparse(S);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Grow the netwoks considering the functional annotations
%--------------------------------------------------------------------------
function [G, S, F] = functionalGrowth(nodeIndex, G, S, F, d, q_con, q_mod, q_new, q_del, delta, s_del, s_f, growth_model, lambda_max, p_fo, kn, tn, Ts, beta)

I=1;
M=length(S);
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
        end
        G(r,:)=E;
        G(:,r)=E';
        G=[G,E2';E2,0];
        G(r,end)=rand(1)<q_con;
        G(end,r)=G(r,end);
        
        for j = 1: nodeIndex - 1
            if ~isempty(S{j, nodeIndex})
                snew = S{j, nodeIndex}(:,r).*(1-rand(size(S{j, nodeIndex},1),1)*lambda_max);  % random scaling factor between [0 lambda_max]
                snew=snew.*(snew>Ts);
                S{j, nodeIndex}=[S{j, nodeIndex},snew];
            end
        end        
        for j = nodeIndex + 1: M
            if ~isempty(S{nodeIndex, j})
                snew = S{nodeIndex, j}(r,:).*(1-rand(1, size(S{nodeIndex, j},2))*lambda_max);  % random scaling factor between [0 lambda_max]
                snew=snew.*(snew>Ts);
                S{nodeIndex, j}=[S{nodeIndex, j};snew];
            end
        end
        F=[F (rand(1)<p_fo)*F(r)];       
    end
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
        end
        G=[G,E';E,0];
        for j = 1: nodeIndex - 1
            if ~isempty(S{j, nodeIndex})
                snew = S{j, nodeIndex}(:,rnd).*(1-rand(size(S{j, nodeIndex},1),1)*lambda_max);
                snew=snew.*(snew>Ts);
                S{j, nodeIndex}=[S{j, nodeIndex},snew];
            end
        end
        for j = nodeIndex + 1: M
            if ~isempty(S{nodeIndex, j})
                snew = S{nodeIndex, j}(rnd,:).*(1-rand(1, size(S{nodeIndex, j},2))*lambda_max);  % random scaling factor between [0 lambda_max]
                snew=snew.*(snew>Ts);
                S{nodeIndex, j}=[S{nodeIndex, j};snew];
            end
        end
        F=[F (rand(1)<p_fo)*F(rnd)];
    end
elseif (strcmp(growth_model,'CG'))
    for i=1:d
        if sum(find(I==i))~=0
            C=Find_modules(G);
            I=ceil(size(G,1)*0.1)+i;
        end
        n=size(G,1);
        pnew=1/length(C);
        if (rand(1)<pnew)
            pd=n+1-sum(G);
            pd=pd/sum(pd);
            Pd=cumsum(pd);
            t=sparse(1,size(G,1));
            while (sum(t)<delta)
                t(find(rand(1)<Pd,1))=1;
            end
            Fnew=0;
            for j = 1: nodeIndex - 1
                if ~isempty(S{j, nodeIndex})
                    nn=size(S{j,nodeIndex},1);
                    ff=min(round(rand(1).^(1/(1-beta))),nn);
                    gg=gamrnd(kn,tn,ff,1)+Ts;
                    ff2=randperm(nn);
                    ff2=ff2(1:ff);
                    snew=sparse(nn,1);
                    snew(ff2)=gg;
                    S{j,nodeIndex}=[S{j,nodeIndex},snew];
                end
            end
            for j = nodeIndex + 1: M
                if ~isempty(S{nodeIndex, j})
                    nn=size(S{nodeIndex, j},2);
                    ff=min(round(rand(1).^(1/(1-beta))),nn);
                    gg=gamrnd(kn,tn,1,ff)+Ts;
                    ff2=randperm(nn);
                    ff2=ff2(1:ff);
                    snew=sparse(1,nn);
                    snew(ff2)=gg;
                    S{nodeIndex, j}=[S{nodeIndex, j}; snew];
                end
            end
            C=[C {n+1}];
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
                end
                t(C{rnd}(anch))=1;
                nn=find(G2(anch,:));
                if length(nn)<=(delta-sum(t))
                    t(C{rnd}(nn))=1;
                else
                    f2=randperm(length(nn));
                    f2=f2(1:(delta-sum(t)));
                    t(C{rnd}(nn(f2)))=1;
                end
            end
            for j=1: nodeIndex - 1
                if ~isempty(S{j, nodeIndex})
                    snew=S{j, nodeIndex}(:,anch_flag).*(1-rand(size(S{j,nodeIndex},1),1)*lambda_max);
                    snew=snew.*(snew>Ts);
                    S{j,nodeIndex}=[S{j,nodeIndex},snew];
                end
            end
            for j = nodeIndex + 1: M
                if ~isempty(S{nodeIndex, j})
                    snew = S{nodeIndex, j}(anch_flag,:).*(1-rand(1, size(S{nodeIndex, j},2))*lambda_max);  % random scaling factor between [0 lambda_max]
                    snew=snew.*(snew>Ts);
                    S{nodeIndex, j}=[S{nodeIndex, j};snew];
                end
            end
            C{f0}=[C{f0} n+1];
        end
        G=[G,t';t,0];
        F=[F Fnew];
    end
elseif (strcmp(growth_model, 'STICKY'))
    beta = 1.6;
    init_n = size(G,1);
    
    for ind=1:d
        x = 1:(init_n + ind - 1);
        y = x.^(-beta);
        y = y./sum(y);
        y_cdf = cumsum(y);
        deg_new = find(histcounts(rand, full([0, y_cdf])));
        
        deg_dist = [sum(G,2); deg_new];
        sticky_ind = deg_dist./(sqrt(sum(deg_dist)));
        sticky_val = sticky_ind(end);
        
        while 1
            % random link between u-v'
            E = rand(length(sticky_ind)-1,1) < (sticky_ind(1:end-1)*sticky_val);
            nbr_new = find(E);
            if length(nbr_new) > 1
                E(nbr_new) = rand(length(nbr_new), 1) > s_del;
                nbr_new = find(E);
            end
            
            if ~isempty(nbr_new)
                break;
            end
        end
        G=[G, E; E',0];
        
        % find anchor node
        [v_anchor, order_anchor] = sort(sticky_ind(nbr_new), 'descend');
        
        % pick up to k candidates
        if length(order_anchor) < s_f
            numberOfCandidates = length(order_anchor);
        else
            numberOfCandidates = s_f;
        end
%         numberOfCandidates = length(order_anchor);
        
        v_anchor = v_anchor(1:numberOfCandidates);
        order_anchor = order_anchor(1:numberOfCandidates);
        
        sticky_candidates = v_anchor./sum(v_anchor);
        cdf_sticky_candidates = cumsum(sticky_candidates);
        index_f = find(histcounts(rand, full([0; cdf_sticky_candidates])));
        index_anchor = nbr_new(order_anchor(index_f));
        
        % assign similarity scores
        for j=1: nodeIndex - 1
            if ~isempty(S{j, nodeIndex})
                snew = S{j, nodeIndex}(:,index_anchor).*(1-rand(size(S{j, nodeIndex},1),1)*lambda_max);  % random scaling factor between [0 lambda_max]
                snew=snew.*(snew>Ts);
                S{j, nodeIndex}=[S{j, nodeIndex},snew];
            end
        end
        for j = nodeIndex + 1: M
            if ~isempty(S{nodeIndex, j})
                snew = S{nodeIndex, j}(index_anchor,:).*(1-rand(1, size(S{nodeIndex, j},2))*lambda_max);  % random scaling factor between [0 lambda_max]
                snew=snew.*(snew>Ts);
                S{nodeIndex, j}=[S{nodeIndex, j};snew];
            end
        end
        F = [F (rand(1)<p_fo)*F(index_anchor)];
    end
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
end
M=nnz(G);
E=(G/2+diag(diag(G)/2))/M;
A=sum(E);
Q=sum(diag(E)'-A.^2);
while (1==1)
    [a,b,j]=find(triu(E,1));
    if isempty(a)
        break
    end
    DQ=[2*j-A(a)'.*A(b)'];
    [mi,ii]=max(DQ);
    if mi<0
        break
    end
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
    end
end
fclose(fid);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Write the Network Functional annotation
%--------------------------------------------------------------------------
function Write_FO(F,a,my_file)
fid = fopen(my_file,'w');

for i=1:length(F)
    if F(i)>0
        fprintf(fid, '%s%d\tFO%d\n', a,i,F(i));
    end
end
fclose(fid);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Write the Networks similarity scores
%--------------------------------------------------------------------------
function Write_Sim(S,a,b,my_file)
fid = fopen(my_file,'w');

[i,j,c]=find(S);
for k=1:length(i)
    fprintf(fid, '%s%d\t%s%d\t%f\n', a,i(k),b,j(k),c(k));
end
fclose(fid);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Write the log file
%--------------------------------------------------------------------------
function Write_log(tree_dictionary,G_set,M,out_path,growth_model,q_con,q_mod,q_new,q_del,delta,s_del,s_f,Na)
fid=fopen(strcat(out_path,'log_file.txt'),'w');
fprintf(fid,'Generated Network Family Information\n');
fprintf(fid,'-------------------------------------\n');
fprintf(fid,'Number of Generated Networks:\t%d\n',M);
fprintf(fid,'Ancestral network size:\t%d\n', Na);
fprintf(fid,'Growth Model:\t%s\n',growth_model);
fprintf(fid,'Growth Parameters:\t');
if strcmp(growth_model,'DMC')
    fprintf(fid,'q_con=%f\t q_mod=%f\t\n',q_con,q_mod);
elseif strcmp(growth_model,'DMR')
    fprintf(fid,'q_new=%f\t q_del=%f\t\n',q_new,q_del);
elseif strcmp(growth_model,'CG')
    fprintf(fid,'delta=%d\t\n',delta);
elseif strcmp(growth_model,'STICKY')
    fprintf(fid,'s_del=%f\t s_f=%f\t\n',s_del,s_f);
end
fprintf(fid,'\nList of Generated Networks:\n');
fprintf(fid,'Name\tLetter\tSize\tInteractions\n');
for i=1:M
    fprintf(fid,'%s\t\t%s\t\t%d\t%d\n',tree_dictionary{i},char(64+i),size(G_set{i},1),sum(sum(triu(full(G_set{i}), 1)~=0)));
end
fclose(fid);


%--------------------------------------------------------------------------
%Display output log
%--------------------------------------------------------------------------
function Disp_log(tree_dictionary,G_set,M,growth_model,q_con,q_mod,q_new,q_del,delta,s_del,s_f,Na)
disp('-------------------------------------');
disp('Generated Network Family Information');
disp('-------------------------------------');
disp(sprintf('Number of Generated Networks:\t%d',M));
disp(sprintf('Ancestral network size:\t%d',Na));
disp(sprintf('Growth Model:\t%s',growth_model));
s='Growth Parameters:';
if strcmp(growth_model,'DMC')
    s=sprintf('%s\tq_con=%f\t q_mod=%f\t',s,q_con,q_mod);
elseif strcmp(growth_model,'DMR')
    s=sprintf('%s,q_new=%f\t q_del=%f\t',s,q_new,q_del);
elseif strcmp(growth_model,'CG')
    s=sprintf('%sdelta=%d\t',s,delta);
elseif strcmp(growth_model,'STICKY')
    s=sprintf('%s,s_del=%f\t s_f=%f\t',s,s_del,s_f);
end
disp(s);
disp('List of Generated Networks:');
disp(sprintf('Name\tLetter\tSize\tInteractions'));
for i=1:M
    disp(sprintf('%s\t%s\t%d\t%d',tree_dictionary{i},tree_dictionary{i},size(G_set{i},1),sum(sum(triu(full(G_set{i}), 1)~=0))));
end


function eTheNumberOfDataset_Callback(hObject, eventdata, handles)
% hObject    handle to eTheNumberOfDataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

numberOfFamilies = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function eTheNumberOfDataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eTheNumberOfDataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pNetworkGrowthModel.
function pNetworkGrowthModel_Callback(hObject, eventdata, handles)
% hObject    handle to pNetworkGrowthModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selected = get(hObject, 'Value');

defaultParameters = getappdata(handles.figure1, 'parameters');

set(handles.sIP1, 'String', defaultParameters.IPLabel{selected}{1});
if length(defaultParameters.IPLabel{selected}) > 1
    set(handles.sIP2, 'Visible', 'on');
    set(handles.eIP2, 'Visible', 'on');
    set(handles.sIP2, 'String', defaultParameters.IPLabel{selected}{2});
else
    set(handles.sIP2, 'Visible', 'off');
    set(handles.eIP2, 'Visible', 'off');
end
set(handles.eIP1, 'String', num2str(defaultParameters.IP(selected, 1)));
set(handles.eIP2, 'String', num2str(defaultParameters.IP(selected, 2)));


% --- Executes during object creation, after setting all properties.
function pNetworkGrowthModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pNetworkGrowthModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function eTheNumberOfNetworks_Callback(hObject, eventdata, handles)
% hObject    handle to eTheNumberOfNetworks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

defaultParameters = getappdata(handles.figure1, 'parameters');
selected = str2double(get(hObject,'String'));

if selected == 2
    predefinedValue = defaultParameters.SS(1);
    set(handles.eSeedNetworkSize, 'String', num2str(predefinedValue));
    set(handles.bPyloGenTreeFilePath, 'String', defaultParameters.treeFileName{1});
elseif selected == 5
    predefinedValue = defaultParameters.SS(2);
    set(handles.eSeedNetworkSize, 'String', num2str(predefinedValue));
    set(handles.bPyloGenTreeFilePath, 'String', defaultParameters.treeFileName{2});
elseif selected == 8
    predefinedValue = defaultParameters.SS(3);
    set(handles.eSeedNetworkSize, 'String', num2str(predefinedValue));
    set(handles.bPyloGenTreeFilePath, 'String', defaultParameters.treeFileName{3});
else
    set(handles.bPyloGenTreeFilePath, 'String', '');
end


% --- Executes during object creation, after setting all properties.
function eTheNumberOfNetworks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eTheNumberOfNetworks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function eInternetworkParameters_Callback(hObject, eventdata, handles)
% hObject    handle to eInternetworkParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eInternetworkParameters as text
%        str2double(get(hObject,'String')) returns contents of eInternetworkParameters as a double


% --- Executes during object creation, after setting all properties.
function eInternetworkParameters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eInternetworkParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function eIP1_Callback(hObject, eventdata, handles)
% hObject    handle to eIP1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eIP1 as text
%        str2double(get(hObject,'String')) returns contents of eIP1 as a double


% --- Executes during object creation, after setting all properties.
function eIP1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eIP1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eSeedNetworkSize_Callback(hObject, eventdata, handles)
% hObject    handle to eSeedNetworkSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eSeedNetworkSize as text
%        str2double(get(hObject,'String')) returns contents of eSeedNetworkSize as a double


% --- Executes during object creation, after setting all properties.
function eSeedNetworkSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eSeedNetworkSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eIP2_Callback(hObject, eventdata, handles)
% hObject    handle to eIP2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eIP2 as text
%        str2double(get(hObject,'String')) returns contents of eIP2 as a double


% --- Executes during object creation, after setting all properties.
function eIP2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eIP2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function defaultAll(hObject, handles)

defaultParameters = getappdata(hObject, 'parameters');

set(handles.eTheNumberOfNetworks, 'string', '2');
set(handles.eTheNumberOfDataset, 'string', '1');
set(handles.eIP1, 'string', defaultParameters.IP(1, 1));
set(handles.eIP2, 'string', defaultParameters.IP(1, 2));
set(handles.sIP1, 'String', defaultParameters.IPLabel{1}{1});
set(handles.sIP2, 'String', defaultParameters.IPLabel{1}{2});
set(handles.eSeedNetworkSize, 'string', num2str(defaultParameters.SS(1)));
set(handles.bPyloGenTreeFilePath, 'String', defaultParameters.treeFileName{1});

guidata(hObject, handles);


% --- Executes on button press in bPyloGenTreeFilePath.
function bPyloGenTreeFilePath_Callback(hObject, eventdata, handles)
% hObject    handle to bPyloGenTreeFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, filePath] = uigetfile({'*.txt'}, 'Search Phylogenetic Tree Files');
set(hObject, 'String', fileName);
