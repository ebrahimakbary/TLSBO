clc
clear

% E. Akbari, M. Ghasemi, M. Gil, A. Rahimnejad, and S. A. Gadsden. 
% "Optimal Power Flow via Teaching-Learning-Studying-Based Optimization
%  Algorithm." Electric Power Components and Systems (2021): 1-18.
% doi: 10.1080/15325008.2021.1971331
% WebLink: https://www.tandfonline.com/doi/abs/10.1080/15325008.2021.1971331

disp('TLSBO (Teaching-Learning-Studying-Based Optimization Algorithm)');

%% TLSBO Parameters

is_TLSBO=1;      % 1 to run TLSBO,  0 to run TLBO

MaxIt = 1000;    % Maximum Number of Iterations

nPop =50;        % Population Size

%% Problem Definition

CostFunction=@(x) CostFunction(x);

VarMin=-100;
VarMax=100;

nVar=30;  % Number of Unknown Variables

%% Initialization

VarSize = [1 nVar]; % Unknown Variables Matrix Size

% Empty Structure for Individuals
empty_individual.Position = [];
empty_individual.Cost = [];

% Initialize Population Array
pop = repmat(empty_individual, nPop, 1);

% Initialize Best Solution
BestSol.Cost = inf;

% Initialize Population Members
for i=1:nPop
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    pop(i).Cost = CostFunction(pop(i).Position);
    
    if pop(i).Cost < BestSol.Cost
        BestSol = pop(i);
    end
end
BestCosts = zeros(MaxIt,1);

if is_TLSBO
    PPosition=zeros(nPop,nVar);
end

%% TLBO Main Loop

for  iter= 1:MaxIt
    
    Mean = 0;
    for i=1:nPop
        Mean = Mean + pop(i).Position;
    end
    Mean = Mean/nPop;
    
    % Select Teacher
    Teacher = pop(1);
    for i=2:nPop
        if pop(i).Cost < Teacher.Cost
            Teacher = pop(i);
        end
    end
    
    %% Teacher Phase
    for i=1:nPop
        newsol = empty_individual;
        
        % Teaching Factor
        TF = randi([1 2]);
        
        if is_TLSBO
            A = 1:nPop;
            for k=1:nVar
                jj = A(randi(nPop-1));
                if pop(i).Cost < pop(jj).Cost
                    PPosition (i,k) = rand*( pop(i).Position(1,k)-pop(jj).Position(1,k));
                else
                    PPosition (i,k) = rand*( pop(jj).Position(1,k)-pop(i).Position(1,k));
                end
            end
            kk1=randn;
            % Teaching
            newsol.Position = pop(i).Position ...
                + rand(VarSize).*( Teacher.Position - TF*Mean)+kk1*PPosition(i,:);
        else
            newsol.Position = pop(i).Position ...
                + rand(VarSize).*( Teacher.Position - TF*Mean); %#ok<*UNRCH>
        end
        
        newsol.Position = max(newsol.Position, VarMin);
        newsol.Position = min(newsol.Position, VarMax);
        
        newsol.Cost = CostFunction(newsol.Position);
        
        if newsol.Cost<pop(i).Cost
            pop(i) = newsol;
            if pop(i).Cost < BestSol.Cost
                BestSol = pop(i);
            end
        end
        
    end
    
    %% Learner Phase
    for i=1:nPop
        
        A = 1:nPop;
        A(i)=[];
        j = A(randi(nPop-1));
        
        Step = pop(i).Position - pop(j).Position ;
        if pop(j).Cost < pop(i).Cost
            Step = -Step;
        end
        
        
        newsol = empty_individual;
        
        if is_TLSBO
            kk=randn;
            newsol.Position = pop(i).Position + rand(VarSize).*Step+kk*PPosition(i,:);
        else
            newsol.Position = pop(i).Position + rand(VarSize).*Step;
        end
        
        newsol.Position = max(newsol.Position, VarMin);
        newsol.Position = min(newsol.Position, VarMax);
        
        newsol.Cost = CostFunction(newsol.Position);
        
        %% Update The Best Solution
        
        if newsol.Cost<pop(i).Cost
            pop(i) = newsol;
            if pop(i).Cost < BestSol.Cost
                BestSol = pop(i);
            end
        end
        
    end
    
    % Store Record for Current Iteration
    BestCosts(iter) = BestSol.Cost+(0);
    
    % Show Iteration Information
    disp(['Iteration ' num2str(iter) ': Best Cost = ' num2str(BestCosts(iter))]);
    
    
end

%% Results

plot((BestCosts),'b--');
