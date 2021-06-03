clear variables
close all

reversalparameters %load parameters associated with cells

plotflag=1; %if 1 plot trajectories, to increase speed of simulation set flag =0

figure

Ntimesteps = tN/dt; %Number of time steps between t0 and tN

X=zeros(Ncells,Ntimesteps);
Y=X;
theta=X;
entry=[];

%for each cell run simulation WITH BIAS
for cell=1:Ncells
    
    %intialise cell orientation
    Theta =rand*2*pi;
    
    %for each time step update velocity and position of cell
    for step=2:Ntimesteps
        
            pR = rand; % pick pR from U[0,1] to determine if cell reverses or not
            if pR < lambdaRb*exp(-vs*ChiR*sin(Theta))*dt % Cell reverses
                entry=[entry,Theta]; %store entry angle
                if plotflag==1
                   plot(X(cell,step-1),Y(cell,step-1),'b*')  %plot reversal as blue star
                end
                hold on
                Theta=Theta+pi;  %reverse direction
            end 
        
    Theta=Theta+randn*(2*Dr*dt)^(1/2); %add rotational noise 
    vx=vs*cos(Theta); % Update x component of velocity
    vy=vs*sin(Theta); % Update y component of velocity
    X(cell,step)=X(cell,step-1)+ vx*dt; % Store cell position at step using new velocity
    Y(cell,step)=Y(cell,step-1)+ vy*dt; % Store cell position at step using new velocity
    theta(cell,step)=Theta; %Store theta at step    
    end
    if plotflag==1  %plot whole track
        plot(X(cell,:), Y(cell,:),'LineWidth', 1.5)
        axis equal
        axis square
        xlabel('X')
        ylabel('Y')
    end
   
end


%%

figure
entry=mod(entry,2*pi);
% histogram(entry);
% axis([0,2*pi,0,800]);


histogram(entry,'Normalization','pdf','BinLimits',[0 2*pi]);
hold on
%plot line at f lambda??
xlabel('Entry Angle $\theta$','Interpreter','latex');
ylabel('Number of reversals away from $\theta$','Interpreter','latex');

%%
%Calculate ftheta for all angles at the end of the simulation
beta=vs*ChiR;
vec=linspace(0,2*pi,Ncells);
for i=1:length(vec)
thdistribution(i)=ftheta(vec(i),beta,lambdaRb);
end 

% Plot theoretical f compared to histogram of simulated angles at last time
% point
figure
Thetanew=mod(theta,2*pi);
histogram(Thetanew(:,[end-100:end]),10,'Normalization','pdf')
axis([0 2*pi -inf inf])
xlabel('$\theta$','Interpreter','latex');
ylabel('$f(\theta)$','Interpreter','latex');
hold on
plot(vec,thdistribution,'r','LineWidth',1.5)

%% Simulation with no bias

figure
%for each cell run simulation
for cell=1:Ncells
    
    %intialise cell orientation
    Theta =rand*2*pi;
    
    %for each time step update velocity and position of cell
    for step=2:Ntimesteps
        
            pR = rand; % pick pR from U[0,1] to determine if cell reverses or not
            if pR < lambdaRb*dt % Cell reverses
                if plotflag==1
                   plot(X(cell,step-1),Y(cell,step-1),'b*')  %plot reversal as blue star
                end
                hold on
                Theta=Theta+pi;  %reverse direction
            end 
        
    Theta=Theta+randn*(2*Dr*dt)^(1/2); %add rotational noise 
    vx=vs*cos(Theta); % Update x component of velocity
    vy=vs*sin(Theta); % Update y component of velocity
    X(cell,step)=X(cell,step-1)+ vx*dt; % Store cell position at step using new velocity
    Y(cell,step)=Y(cell,step-1)+ vy*dt; % Store cell position at step using new velocity
    theta(cell,step)=Theta; %Store theta at step    
    end
    if plotflag==1  %plot whole track
        plot(X(cell,:), Y(cell,:),'LineWidth', 1.5)
        axis equal
        axis square
        xlabel('X')
        ylabel('Y')
    end
   
end

%%
%Add sublabels
%htext=sublabel(h,0,-16);
 
%% %Save figure
%print('reversalsfig2', '-djpeg', '-r500')




%%

% T=[dt:dt:tN];
% meanY= mean(Y,1);
% coeffsU=polyfit(T,meanY,1); %line fit, drift is gradient
% Ydrift=polyval(coeffsU,T);
% drift=coeffsU(1)
% 
% %calculate mean squared displacement and estimate diffusion coefficient
% diff=(Y-meanY).^2;
% MSD=mean(diff,1);
% coeffsD=polyfit(T,MSD,1); %line fit to estimate diffusion coefficient 
% Ydiff=polyval(coeffsD,T);
% diffusion=coeffsD(1)
% 
% if plotflag==1 
%     figure
% end
% %plot mean Y position against time
% plot(T,meanY) %cells start at origin
% hold on
% plot(T,Ydrift) %plot line fit for drift estimate
% xlabel('t')
% ylabel('Mean position (y)')
% 
% if plotflag==1 
%     figure
% end
% 
% %plot MSD
% plot(T,MSD)
% hold on 
% plot(T,Ydiff)
% xlabel('t')
% ylabel('Mean Squared Displacement')
% 
% %Plot snapshot histograms of orientation distribution
% Thetanew=mod(theta,2*pi);
% %Thetanew=Thetamod-pi;
% figure
% subplot(2,2,1)
% histogram(Thetanew(:,10),5) %counts only
% xlabel('Theta')
% ylabel('f(theta) T=10')
% subplot(2,2,2)
% histogram(Thetanew(:,500),5) %counts only
% xlabel('Theta')
% ylabel('f(theta) T=500')
% subplot(2,2,3)
% histogram(Thetanew(:,1000),5) %counts only
% xlabel('Theta')
% ylabel('f(theta) T=1000')
% subplot(2,2,4)
% histogram(Thetanew(:,10000),5) %counts only
% xlabel('Theta')
% ylabel('f(theta) T=10000')
% 
% %%
% %Calculate ftheta for all angles at the end of the simulation
% beta=vs*ChiR;
% %for i=1:Ncells
% %distribution(i)=ftheta(Thetanew(i,end),beta,lambdaRb);
% %end 
% %vec=linspace(-pi,pi,Ncells);
% vec=linspace(0,2*pi,Ncells);
% for i=1:length(vec)
% thdistribution(i)=ftheta(vec(i),beta,lambdaRb);
% end 
% 
% % Plot theoretical f compared to histogram of simulated angles at last time
% % point
% figure
% histogram(Thetanew(:,end),10,'Normalization','pdf')
% axis([0 2*pi -inf inf])
% xlabel('Theta')
% ylabel('f(theta)')
% hold on
% plot(vec,thdistribution,'r','LineWidth',1.5)
% 
% 
% %Calculate py analytically
% pyanalytic=integral(@(theta)ftheta(theta,beta,lambdaRb).*sin(theta),0,2*pi);
% 
% %Calculate py from the data
% for i=1:Ncells
% sintheta(i)=sin(Thetanew(i,end));
% end
% pydata=sum(sintheta)/Ncells;
% 
% 
% 
