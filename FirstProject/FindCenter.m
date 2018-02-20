function M_final = FindCenter(data, plotflag)
    switch nargin
        case 2
            
        case 1
            plotflag = 0;
        otherwise
            error('WRONG INPUT!');
    end
    %% LOAD DATA
    % load('sc_neighb_cell_pnts_cott.mat');
    n = size(data,2);
    center=[];
    for i=1:n
       X = delteOutlier(unique(data{:,i},'rows'));
       if(~isempty(X))
           [clusArray, M,~] = FindCenterCell(X);

            %% FIGURES
            if(plotflag)
                p1 = gscatter(X(:,1),X(:,2),clusArray);
                hold on;
            end
            if(size(M,2) == 2)
                center = [center ; M] ;
            end

       end
    end

%     center = delteOutlier(center);
    plot(center(:,1),center(:,2),'+r','MarkerSize',10)
    M_final(:,1) = median(center(:,1));
    M_final(:,2) = median(center(:,2));

   if(plotflag)
        p2 = plot(M_final(:,1),M_final(:,2),'*r','MarkerSize',10);
        legend(p2,'Center of concentric circles');
   end
    %% CLEAR
    % clearvars -except X clusArray M
end

function [clusArray, M,Par] = FindCenterCell(X)
    %% CLUSTERING
    [n,m,~]=size(X);
    a=reshape(X,n,[],1);
    b=reshape(a(:),n*m,[])';
    c=unique(b,'rows')';
    X=reshape(c,n,m,[]);
    [CLUSTERSs, CentroidPoints] = func_GDD(X);
    
    CLUSTERS = CLUSTERSs;
    maxClustNo = size(CLUSTERS,2);
    N=size(X,1);
    D=size(X,2);
    for i = 1: maxClustNo
        if(size(CLUSTERS{i},2) < 10)
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

    %% FIND THE CENTRE
    Par =[];
    for i=1:size(CLUSTERS,2)
        dataCluster = X(clusArray==i,:);
        if (~isempty(CircleFitByPratt(dataCluster)))
            Par = [Par ;CircleFitByPratt(dataCluster)];
        end
    end
    M = mode(Par(:,1:2));
    
end

function XX = delteOutlier(X)
    D = pdist(X);
    Z = squareform(D);
    for i=1:size(Z,2)
        Z(i,i)=9999;
    end
    minvector = (min(Z));
    idx = isoutlier(minvector);
    XX = X; XX(idx==1,:) = [];

end

function Par = CircleFitByPratt(XY)

%--------------------------------------------------------------------------
%  
%     Circle fit by Pratt
%      V. Pratt, "Direct least-squares fitting of algebraic surfaces",
%      Computer Graphics, Vol. 21, pages 145-152 (1987)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this fit does not use built-in matrix functions (except "mean"),
%           so it can be easily programmed in any programming language
%
%--------------------------------------------------------------------------

n = size(XY,1);      % number of data points

centroid = mean(XY);   % the centroid of the data set

%     computing moments (note: all moments will be normed, i.e. divided by n)

Mxx=0; Myy=0; Mxy=0; Mxz=0; Myz=0; Mzz=0;

for i=1:n
    Xi = XY(i,1) - centroid(1);  %  centering data
    Yi = XY(i,2) - centroid(2);  %  centering data
    Zi = Xi*Xi + Yi*Yi;
    Mxy = Mxy + Xi*Yi;
    Mxx = Mxx + Xi*Xi;
    Myy = Myy + Yi*Yi;
    Mxz = Mxz + Xi*Zi;
    Myz = Myz + Yi*Zi;
    Mzz = Mzz + Zi*Zi;
end
   
Mxx = Mxx/n;
Myy = Myy/n;
Mxy = Mxy/n;
Mxz = Mxz/n;
Myz = Myz/n;
Mzz = Mzz/n;

%    computing the coefficients of the characteristic polynomial

Mz = Mxx + Myy;
Cov_xy = Mxx*Myy - Mxy*Mxy;
Mxz2 = Mxz*Mxz;
Myz2 = Myz*Myz;

A2 = 4*Cov_xy - 3*Mz*Mz - Mzz;
A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz2 - Myz2 - Mz*Mz*Mz;
A0 = Mxz2*Myy + Myz2*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy;
A22 = A2 + A2;

epsilon=1e-12; 
ynew=1e+20;
IterMax=20;
xnew = 0;

%    Newton's method starting at x=0

for iter=1:IterMax
    yold = ynew;
    ynew = A0 + xnew*(A1 + xnew*(A2 + 4.*xnew*xnew));
    if (abs(ynew)>abs(yold))
        disp('Newton-Pratt goes wrong direction: |ynew| > |yold|');
        xnew = 0;
        break;
    end
    Dy = A1 + xnew*(A22 + 16*xnew*xnew);
    xold = xnew;
    xnew = xold - ynew/Dy;
    if (abs((xnew-xold)/xnew) < epsilon), break, end
    if (iter >= IterMax)
        disp('Newton-Pratt will not converge');
        xnew = 0;
    end
    if (xnew<0.)
        fprintf(1,'Newton-Pratt negative root:  x=%f\n',xnew);
        xnew = 0;
    end
end

%    computing the circle parameters

DET = xnew*xnew - xnew*Mz + Cov_xy;
Center = [Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;

Par = [Center+centroid , sqrt(Center*Center'+Mz+2*xnew)];

end    %    CircleFitByPratt

%GDD Clustering Algorithm

%--------------------------------------------------------------------------
%Author:  Emre G�ng�r, Ahmet �zmen
% Please Cite as (if you use in you study or work)
%       G�ng�r, Emre, and Ahmet �zmen. "Distance and density based 
%       clustering algorithm using Gaussian kernel." Expert Systems with 
%       Applications 69(2017): 10-20.
%
% Please Note:
%       Standard error of the mean is used in variables DevD,GPD_CPU,
%       DPD_CPU, which is sqrt(1/N)*standard deviation. It is used for  
%       estimation of the standard deviation of a sampling distribution.
%       Reason for that to correctly estimate relation between samples and
%       region that samples are resides.
%       
%       Also SIP deviation calculations are computed using 
%       std=sqrt((1/N)*(sum(xi^2)-(sum(xi)^2))/N);
%       for performance and coding ease. Other deviations can be calculated
%       this way for performance in loops instead of std() function. 
%   
%       Code is not optimized - And can contain errors please help us
%       if you see any errors by sending e-mail to us. Thanks in advance :)
%       
% Contact Info:
%       Emre G�ng�r: emregungor84@hotmail.com
%       Ahmet �zmen: ozmen@sakarya.edu.tr
%--------------------------------------------------------------------------


%% GDD Algorithm
function [CLUSTERS,centroids]=func_GDD(DataSet)
%  tic
 % CLUSTERS=[];
% centroids=[];
%DataSet=[y x];

%% FILTERING-----------------------------------------------------------------
X_CPU=DataSet;
N=size(X_CPU,1); %Number of Data Samples

%DUPLICATE REMOVAL
%     %ALL POINTS IN X Needed to be unique (No Duplicate Points)
%     [n,m,p]=size(X_CPU);
%     a=reshape(X_CPU,n,[],1);
%     b=reshape(a(:),n*m,[])';
%     % c=unique(b,'rows','stable')';
%     c=unique(b,'rows')';
%     % c=unique(b)';
%     X_CPU=reshape(c,n,m,[]);
A = table(X_CPU);
[C, indexA, indexC]=unique(A);
X_CPU=table2array(C);

% Return if N=2 as 2 cluster
if(size(X_CPU,1)==2 || size(X_CPU,1)<= size(DataSet,1)/100)
    for i=1:size(X_CPU,1)
        CLUSTERS{i}=[i];
        centroids(i)=i;
    end
    %Re-order Cluster Indices According to Original DuplicateDataSet if Any
        for i=1:size(CLUSTERS,2)
            ClustMemCount=size(CLUSTERS{i},2);
            for j=1:ClustMemCount
                OrgIndices=find(indexC(:)==CLUSTERS{i}(j));
                CLUSTERS{i}(j)=[];
                for k=1:size(OrgIndices,1)
                    CLUSTERS{i}(end+1)=OrgIndices(k);
                end
            end
            centroidSelections =find(indexC(:)==centroids(i));
            centroidPointSize=size(centroidSelections ,1);
            centroidPointIndex  =1;
            if(centroidPointSize>1)
                centroidPointIndex  =(round(centroidPointSize/2));
            end
            centroids(i)=centroidSelections (centroidPointIndex  );
        end
    return;
end

%MINIMUM AS GROUND ZERO 
%mean of min max interval to eliminate shifting artifacts
N=size(X_CPU,1);
D=size(X_CPU,2);
for j=1:D
    minX=min(X_CPU(:,j));
    for i=1:N
        X_CPU(i,j)=X_CPU(i,j)-minX;
    end
end 

N=size(X_CPU,1); %Number of Data Samples
D=size(X_CPU,2); %Dimension of Data Space
%--------------------------------------------------------------------------


%% MEAN & DEVIATION
%Mean for all dimensions
meanD=zeros(D,1);

%Deviation for all dimensions
DevD=zeros(D,1);

%Mean
for i=1:N 
    for j=1:D
        meanD(j)=meanD(j)+X_CPU(i,j);
    end   
end
for j=1:D
    meanD(j)=(meanD(j)/N);
end

%Dimensional Deviance Coeff (DevD)= Deviance * sqrt(1/N)
for i=1:N
    for j=1:D
        DevD(j)=DevD(j)+(X_CPU(i,j)-meanD(j))^2;
    end
end
for j=1:D
    DevD(j)=(1/N)*sqrt(DevD(j)); 
    %     DevD(j)=sqrt(1/N)*sqrt((1/N)*DevD(j)); 
end
%___________________________________________________________
%% REMOVE UNNECESSARY DIMENSIONS
EraseDimension=[];

for j=1:D 
    %--------------------
    if(DevD(j)==0) 
       %Erase Dimension Properties 
       EraseDimension(end+1)=j;
    end
end
%---------------------------------------------------------------------
%Erase Unnecessary Dimensions for Clustering--------------EXTREME CASES
for i=size(EraseDimension,2):-1:1
    X_CPU(:,EraseDimension(i))=[];
    meanD(EraseDimension(i))=[];
    DevD(EraseDimension(i))=[];
end
D=size(X_CPU,2); %Properties
% meanD
% DevD

%---------------------------------------------------------------------

D=size(X_CPU,2); %Dimension of Data Space


if(D==0) %If there is no data to be culstered close algorithm
    disp('No clustering needed. Recommendation: Either take all data one cluster or all sample points are clusters');
    CLUSTERS=[];
    centroids=[];
    return;
end
 
 
%% CLUSTERING 
 
%% STEP 1- STATISTICAL PROPERTY COMPUTATIONS
 
    DM_CPU=zeros(N,N);
    GM_CPU=zeros(N,N);
    AG_CPU=zeros(N,1)';
    GPM_CPU=zeros(N,1)';
    GPD_CPU=zeros(N,1)';
    DPM_CPU=zeros(N,1)';
    DPD_CPU=zeros(N,1)';
    for i=1:N
        for j=1:N
            tempDist=0;
            tempExp=0;
            for k=1:D
                tempDist=tempDist+(X_CPU(i,k)-X_CPU(j,k))^2;
                tempExp=tempExp+ ((X_CPU(i,k)-X_CPU(j,k))^2)/(2*( sqrt((meanD(k)*DevD(k))/(2*pi) ) )^2);
            end
            tempDist=sqrt(tempDist);
            DM_CPU(i,j)=tempDist;
            GM_CPU(i,j)=exp(-1*tempExp);
            AG_CPU(i)=AG_CPU(i)+GM_CPU(i,j);
        end
        GPM_CPU(i)=(AG_CPU(i)-1)/(N-1);
        DPM_CPU(i)=sum(DM_CPU(i,:))/(N-1);
    end
    
    for i=1:N
        for j=1:N 

                GPD_CPU(i)= GPD_CPU(i)+((GM_CPU(i,j)-GPM_CPU(i))^2); 
                DPD_CPU(i)= DPD_CPU(i)+((DM_CPU(i,j)-DPM_CPU(i))^2); 

        end
        GPD_CPU(i) = (1/(N))*sqrt(GPD_CPU(i));
        DPD_CPU(i) = (1/(N))*sqrt(DPD_CPU(i));
    end
    

    %---------------------------
    SpaceGaussDeviance=mean(GPD_CPU);
    SpaceGaussMean=mean(GPM_CPU);
    SpaceDistanceDeviance=mean(DPD_CPU);
    SpaceAGMean=mean(AG_CPU);
    %---------------------------    
  

%% STEP 2 - CLUSTERING PROCESS    
CLUSTERS={};

DataSetIndices=1:N;
centroids=[];
    %Centroid: Cluster Centroid
    %take maximum peak energy of gaussian to pinpoint maximum density in space
    [MaxGaussianMean,MaxGaussIndex] = max(GPM_CPU);
    MaxGaussianVal=MaxGaussianMean*N; 
    centroidIndex =MaxGaussIndex;
    centroidAccumulativeGaussian=MaxGaussianVal;
   
%LOOP - AFTER FINDING CLUSTER FIND NEXT CLUSTER THAT IS MOST DENSE
while(size(DataSetIndices,2)~=0)
    
    %Cluster Centroids
    centroids(end+1)=centroidIndex ;

    
% % %    % Centroid (Cluster Centroid) DISTANCE CALCULATIONS
%--------------------------------------------------------------------------------------------------------------

  % CLUSTER FIXED COMPONENTS (EQ. 8-9)
  FGDT=(GPM_CPU(centroidIndex )/(pi*sqrt(GPD_CPU(centroidIndex ))))-(pi*GPD_CPU(centroidIndex ));
  FDT=abs((DPM_CPU(centroidIndex ))/(log(GPM_CPU(centroidIndex )*N+DPD_CPU(centroidIndex ))*(log(GPD_CPU(centroidIndex )))));

  %DEVIANCE FOR EACH NEIGHBOUR TO ADJUST GRADUAL CHANGES
%--------------------------------------------------------------------------------------------------------------    
    %NSL creation and Centroid (Cluster Centroid) Assignments
    ClusterIndices=[];
    NSL=[]; %Neighbor Search List
    NSL(end+1)=centroidIndex ;
    ClusterIndices(end+1)=centroidIndex ;

    %Remove Cluster Centroid from Unclustered Sample List
    unclusteredCentroidID =find(DataSetIndices(:)==centroidIndex );
    DataSetIndices(unclusteredCentroidID )=[]; %Centroid Assigned to Cluster Remove from unclustered List

    %SPL Mean and Deviation Variables
    %----------------------------------------------------------------------
    %Variance of Neighborhood variances
    MeanBase=0; %Ground Zero Mean
    Mg=0;
    Vg=0; %Ground Zero Variance
    Diffcounter=0;
    %Gaussian CounterPart
    GMeanBase=0;GMg=0;GVg=0;
    %----------------------------------------------------------------------

% PAUSE for Benchmarking and Debugging process    
%     pause
    
    %Main Clustering Loop (figure 3)
    while (1)
        NSL_new=[];

        if (size(NSL,2)~=0)
            %Inner Loop of Flow-Chart (figure 3)
            for i=size(NSL,2):-1:1
                NPs=[]; %Neighbor Sample Points
                
                %SPL (Samples In Proximity List) Calculations______________
                %Note: SPL members can be found and calculations can be 
                %      done afterwards ....
                    PDistMean=0;
                    PGaussMean=0;
                    cntr=0;
                    
                    for j=1:N 
                        if(DM_CPU(NSL(i),j)<=FDT)
                            PDistMean=PDistMean+DM_CPU(NSL(i),j);
                            PGaussMean=PGaussMean+GM_CPU(NSL(i),j);
                            cntr=cntr+1;
                        end
                    end
                    if(cntr~=0)
                        PDistMean=PDistMean/cntr;
                        PGaussMean=PGaussMean/cntr;
                    end


                    cntr=0;
                    PDistDev=0;
                    PGaussDev=0;
                    for j=1:N 
                        if(DM_CPU(NSL(i),j)<=FDT)
                            PDistDev=PDistDev+(PDistMean-DM_CPU(NSL(i),j))^2;
                            PGaussDev=PGaussDev+(PGaussMean-GM_CPU(NSL(i),j))^2;
                            cntr=cntr+1;
                        end
                    end
                    if(cntr~=0)
                        PDistDev=sqrt(PDistDev/cntr);
                        PGaussDev=sqrt(PGaussDev/cntr);
                        %Variance of Neighborhood Variances
                        Diffcounter=Diffcounter+1;
                        MeanBase=MeanBase+PDistDev;
                        Vg=Vg+PDistDev^2;
                        
                        GMeanBase=GMeanBase+PGaussDev;
                        GVg=GVg+PGaussDev^2;
                    end
                        
                       %Std Deviance Using N loop
                       %-----------------------------
%                      %std=sqrt((1/N)*(sum(xi^2)-(sum(xi)^2))/N);
                       %----------------------------
                       %Calculations below are standard deviations of
                       %variances
                       Mg= (MeanBase)^2/Diffcounter;
                       StdVariance=sqrt((1/(Diffcounter))*abs(Vg-Mg));
                       MeanVariance=MeanBase/Diffcounter;
                       
                       GMg=(GMeanBase)^2/Diffcounter;
                       GStdVariance=sqrt((1/(Diffcounter))*abs(GVg-GMg));
             %___________________________________________________________________________________________
              %GRADIENT CHANGE EQUATIONS (EQ. 10 -11 )  
              GDT=(MeanVariance*(AG_CPU(NSL(i))/AG_CPU(centroidIndex )))+2*pi*StdVariance*(SpaceGaussDeviance/GPD_CPU(NSL(i)));
              GGDT=GStdVariance;
              
              %THRESHOLDS
              A=(FGDT-GGDT)/2;
              B=(FDT+GDT)/2;
              %___________________________________________________________________________________________
               
                for k=size(DataSetIndices,2) :-1:1
                    if ((GM_CPU(NSL(i),DataSetIndices(k))>=A )&& DataSetIndices(k)~=NSL(i))
                       if((DM_CPU(NSL(i),DataSetIndices(k))<=(B) ) && DataSetIndices(k)~=NSL(i))
                            NPs(end+1)=DataSetIndices(k); %point k is in the neighbourhood and same cluster
                            DataSetIndices(k)=[]; %(UnclusteredList) Remove Cluster Members from unclustered list 
                       end
                    end
                end

                %ADD Neighbors to Cluster
                if (size(NPs,2)~=0)
                    ClusterIndices=[ClusterIndices,NPs];
                    NSL_new=[NSL_new,NPs];
                end
            end
            NSL=NSL_new; %New Found Neighbors are now NSL
        else
            %There is no neighbours that fit similarity measures
            break;
        end
    end

    %FOR NEXT CLUSTERING - SETTING
    CLUSTERS{end+1}=ClusterIndices;

    %Set New Maximum Gauss Point - NEXT MOST DENSE POINT
    MaxGaussianMean=0;
    MaxGaussIndex=0;
    for i=1:size(DataSetIndices,2)
        if(GPM_CPU(DataSetIndices(i))>MaxGaussianMean)
            MaxGaussianMean=GPM_CPU(DataSetIndices(i));
            MaxGaussIndex=DataSetIndices(i);
        end
    end
    centroidIndex =MaxGaussIndex;
    
    %Special Case - No Centroid found & 1-Sample Remains
    if(centroidIndex ==0 && size(DataSetIndices,2)==1)
        %OPTION: Can add 1 sample to new cluster
        CLUSTERS{end+1}=DataSetIndices(1);
        centroids(end+1)=DataSetIndices(1);
        break; 
    end
    %Special Case - No Centroid found & N-Sample Remains
    if(centroidIndex ==0 && size(DataSetIndices,2)>1)
%        size(DataSetIndices,2) 
%         DataSetIndices %Sequential Cluster
        for n=1:size(DataSetIndices,2)
%             if(DataSetIndices(n)+1==DataSetIndices(n+1))
            CLUSTERS{end+1}=DataSetIndices(1);  
            centroids(end+1)=DataSetIndices(1);
            DataSetIndices(1)=[];
        end
    end
    if(centroidIndex ==0)
        break; 
    end
end %END CLUSTERING

% t_CPU=toc,
% disp('GDD Function Completed');


%% POST-PROCESS STEP
%Re-order Cluster Indices According to Original DuplicateDataSet if Any
% disp([ 'Cluster Number: ' num2str(size(CLUSTERS,2))]);
CLUSTERSnew=cell(size(CLUSTERS,2),1)';
for i=1:size(CLUSTERS,2)
    ClustMemCount=size(CLUSTERS{i},2);
%     disp([ 'Cluster Number: ' num2str(i) ' ClusterMember: ' num2str(ClustMemCount)]);
%     disp([ 'Total Index Size '  num2str(size(indexC)) ]);
    for j=1:ClustMemCount
        OrgIndices=find(indexC(:)==CLUSTERS{i}(j));
%         CLUSTERS{i}(j)=[];
        for k=1:size(OrgIndices,1)
            CLUSTERSnew{i}(end+1)=OrgIndices(k);
        end
    end
    centroidSelections =find(indexC(:)==centroids(i));
    centroidPointSize=size(centroidSelections ,1);
    centroidPointIndex  =1;
    if(centroidPointSize>1)
        centroidPointIndex  =(round(centroidPointSize/2));
    end
    centroids(i)=centroidSelections (centroidPointIndex  );
    
end
CLUSTERS=CLUSTERSnew;

end
    
    
    
    