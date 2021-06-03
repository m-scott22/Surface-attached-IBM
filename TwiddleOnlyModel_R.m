clear variables
close all

twiddleparameters %load parameters associated with cells

figure

plotflag=1; %if 1 plot trajectories, to increase speed of simulation set flag =0

Ntimesteps = tN/dt; %Number of time steps between t0 and tN

X=zeros(Ncells,Ntimesteps);
Y=X;
theta=X;
entry=[];

for cell=1:Ncells
   
    %intialise cell orientation 
    Theta = rand*2*pi;
    theta(cell,1)=Theta; %Store theta at initial step 
    %for each time step update velocity and position of cell
    for step=2:Ntimesteps
                pT = rand; %pick pT from U[0,1] to determine if cell enters twiddle
                if pT < lambdaTb*exp(-vs*ChiT*sin(Theta))*dt %Cell twiddles
                    entry=[entry,Theta];
                    if plotflag==1
                        plot(X(cell,step-1),Y(cell,step-1),'r*') %plot twiddle as red star
                    end
                    hold on
                    %choose exit angle biased around theta=pi/2
                    theta_range=[0:0.001:2*pi];
                    weight=exp(kappa*cos(theta_range-pi/2)); %von-Mises for bias 
                    Theta = randsample(theta_range,1,true,weight);%%Choose new random angle out of twiddle
                end
       Theta=Theta+randn*(2*Dr*dt)^(1/2); %add rotational noise 
       vx=vs*cos(Theta); % Update x component of velocity
       vy=vs*sin(Theta); % Update y component of velocity
       X(cell,step)=X(cell,step-1)+ vx*dt; % Store cell position at step using new velocity
       Y(cell,step)=Y(cell,step-1)+ vy*dt; % Store cell position at step using new velocity
       theta(cell,step)=Theta; %Store theta at step
    end              
  if plotflag==1  %Plot whole track
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
histogram(entry,'Normalization','pdf','BinLimits',[0 2*pi]);
hold on
vec=linspace(0,2*pi,Ncells);
for i=1:length(vec)
    thdistribution(i)=von_mises(vec(i),kappa);
end 
plot(vec,thdistribution,'r','LineWidth',1.5)

%axis([0,2*pi,0,800]);
xlabel('Entry Angle $\theta$','Interpreter','latex');
ylabel('Number of twiddles away from $\theta$','Interpreter','latex');

%%
%Calculate ftheta for all angles at the end of the simulation
beta=vs*ChiT;
for i=1:length(vec)
    thdistribution(i)=ftheta_new(vec(i),beta,lambdaTb,kappa);
end 

% Plot theoretical f compared to histogram of simulated angles at last time
% point
figure
Thetanew=mod(theta,2*pi);
histogram(Thetanew(:,end),10,'Normalization','pdf','BinLimits',[0 2*pi])
%axis([0 2*pi -inf inf])
xlabel('$\theta$','Interpreter','latex');
ylabel('$f(\theta)$','Interpreter','latex');
hold on
plot(vec,thdistribution,'r','LineWidth',1.5)

% % Plot theoretical f compared to histogram over time course
% 
% plot(vec,thdistribution,'r','LineWidth',1.5)
% hold on
% axis([0 2*pi -inf inf])
% xlabel('$\theta$','Interpreter','latex');
% ylabel('$f(\theta)$','Interpreter','latex');
% for step=1:Ntimesteps
%     figure(4)
%     histogram(Thetanew(:,step),10,'Normalization','pdf','BinLimits',[0 2*pi])
%     pause(0.1)
% end






%% Simulation no bias
X=zeros(Ncells,Ntimesteps);
Y=X;
figure
for cell=1:Ncells
   
    %intialise cell orientation 
    Theta = rand*2*pi;
        
    %for each time step update velocity and position of cell
    for step=2:Ntimesteps
                pT = rand; %pick pT from U[0,1] to determine if cell enters twiddle
                if pT < lambdaTb*dt %Cell twiddles
                    if plotflag==1
                        plot(X(cell,step-1),Y(cell,step-1),'r*') %plot twiddle as red star
                    end
                    hold on
                    %choose exit angle at random from unit circle               
                Theta = rand*2*pi; %Choose new random angle out of twiddle
                end
       Theta=Theta+randn*(2*Dr*dt)^(1/2); %add rotational noise 
       vx=vs*cos(Theta); % Update x component of velocity
       vy=vs*sin(Theta); % Update y component of velocity
       X(cell,step)=X(cell,step-1)+ vx*dt; % Store cell position at step using new velocity
       Y(cell,step)=Y(cell,step-1)+ vy*dt; % Store cell position at step using new velocity
       theta(cell,step)=Theta; %Store theta at step
    end              
  if plotflag==1  %Plot whole track
        plot(X(cell,:), Y(cell,:),'LineWidth', 1.5)
        axis equal
        axis square
        xlabel('X')
        ylabel('Y')
  end
end 

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
% figure
% subplot(2,2,1)
% histogram(Thetanew(:,10),5)
% xlabel('Theta')
% ylabel('f(theta) T=10')
% subplot(2,2,2)
% histogram(Thetanew(:,500),5)
% xlabel('Theta')
% ylabel('f(theta) T=500')
% subplot(2,2,3)
% histogram(Thetanew(:,1000),5)
% xlabel('Theta')
% ylabel('f(theta) T=1000')
% subplot(2,2,4)
% histogram(Thetanew(:,10000),5)
% xlabel('Theta')
% ylabel('f(theta) T=10000')
