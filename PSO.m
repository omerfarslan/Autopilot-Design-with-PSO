function out = PSO(problem,param,des_performance,weights)

%% Problem Definition

weights_cost = weights;
des_perfomance_cost = des_performance;
CostFunction = problem.CostFunction;        % Cost Function

nVar = problem.nVar;            % Number of Decision Variables

VarSize = [1 nVar];   % Size of Decision Variables Matrix

VarMinf = problem.VarMinf;         % Lower Bound of Variables
VarMaxf = problem.VarMaxf;         % Upper Bound of Variables


%% PSO Parameters

MaxIt = param.MaxIt;      % Maximum Number of Iterations

nPop = param.nPop;        % Population Size (Swarm Size)

% PSO Parameters
w = param.w;            % Inertia Weight
wdamp = param.wdamp;     % Inertia Weight Damping Ratio
c1 = param.c1;         % Personal Learning Coefficient
c2 = param.c2;         % Global Learning Coefficient


%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=inf;

for i=1:nPop
    
    x = unifrnd(-0.02,0); % from paper ka is our independent variable
    VarMin = VarMinf;
    VarMax = double(vpa(subs(VarMaxf,x),4));
    particle(i).Position = vpa(unifrnd(VarMin,VarMax),4);
    
    %Initialize Position
    particle(i).Position = vpa(unifrnd(VarMin,VarMax),4); % randomly distrubeted
    
    %Initialize Velocity
    particle(i).Velocity = zeros(VarSize);
    
    % Evaluation
    particle(i).Cost = CostFunction(double(particle(i).Position),des_perfomance_cost,weights_cost);
    
    % Update Personal Best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest = particle(i).Best;
        
    end
    
end

BestCost=zeros(MaxIt,1);

%% PSO Main Loop

bestcost_array =[99999];

for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Velocity Limits
        VelMax = 0.1*(VarMax-VarMin);
        VelMin = -VelMax;
        
        %         Velocity limits changes at every particle since we have position
        %         dependency from upper ki bound
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Check that position will not go away with high velocity update
        
        % Evaluation
        particle(i).Cost = CostFunction(double(particle(i).Position),des_perfomance_cost,weights_cost);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest = particle(i).Best;
                
            end
           
        end
        
    end
    
    BestCost(it)=GlobalBest.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    bestcost_array = [bestcost_array BestCost(it)];
    bestsoln(it) = GlobalBest;
    
    w=w*wdamp;
    
    % Divide Max_iteration to 5 and reset PSO in case of local minima
    % convergence
    
    if it == MaxIt/5 || it == 2*MaxIt/5 || it == 3*MaxIt/5 || it == 4*MaxIt/5
        
        w = param.w;
        GlobalBest.Cost=inf;
        
        for i=1:nPop
            
            x = unifrnd(-0.02,0); % from paper ka is our independent variable
            VarMin = VarMinf;
            VarMax = double(vpa(subs(VarMaxf,x),4));
            particle(i).Position = vpa(unifrnd(VarMin,VarMax),4);
            
            %Initialize Position
            particle(i).Position = vpa(unifrnd(VarMin,VarMax),4); % randomly distrubeted
            
            %Initialize Velocity
            particle(i).Velocity = zeros(VarSize);
            
            % Evaluation
            particle(i).Cost = CostFunction(double(particle(i).Position),des_perfomance_cost,weights_cost);
            
            % Update Personal Best
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest = particle(i).Best;
                
            end
            
        end
        
        BestCost=zeros(MaxIt-it,1);
        
    end
  
end

out.BestCost_Array = bestcost_array;
[out.BestCost ind] = min(bestcost_array);
out.BestSoln = bestsoln(ind-1);
end


