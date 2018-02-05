close all; clear all; clc
%% LOAD DATA
load('sc_neighb_cell_pnts_cott.mat');
fullData = sc_neighb_cell_pnts;
n = size(fullData,2);
X =[];
for i=1:n
   X = [X ; delteOutlier(unique(fullData{:,i},'rows'))];
end

%% SCATTER PLOT 2D
ClustNo=[]; 
[n,m,p]=size(X);
a=reshape(X,n,[],1);
b=reshape(a(:),n*m,[])';
c=unique(b,'rows')';
X=reshape(c,n,m,[]);
%% CLUSTERING
[CLUSTERSs, CentroidPoints] = func_GDD(X);

CLUSTERS = CLUSTERSs;
maxClustNo = size(CLUSTERS,2);
N=size(X,1);
D=size(X,2);
for i = 1: maxClustNo
    if(size(CLUSTERS{i},2) < 100)
        CLUSTERS{:,i} =[];
    end
end
CLUSTERS = CLUSTERS(~cellfun(@isempty, CLUSTERS));

maxClustNo = size(CLUSTERS,2);
N=size(X,1);
D=size(X,2);
clusArray = zeros(N,1);
for i=1: maxClustNo
    clusArray(CLUSTERS{i},1) = i;
end

if(size(CLUSTERS,2)==0 && size(CentroidPoints,2)==0)
    disp('No cluster found exiting');
    return;
end

%% FIGURES
figure;
p1 = gscatter(X(:,1),X(:,2),clusArray);
hold on;

%% FIND THE CENTRE
Par =[];
for i=1:size(CLUSTERS,2)
    dataCluster = X(clusArray==i,:);
    Par = [Par ;CircleFitByPratt(dataCluster)];
end
M = mode(Par(:,1:2));

%% PLOT AGAIN
p2 = plot(M(1),M(2),'*r','MarkerSize',12);
legend(p2,'Center of concentric circles');

%% CLEAR
clearvars -except X clusArray Mœ