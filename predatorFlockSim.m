%{
Flocks and Predators
Self-propelled particle model of aggregation in two dimensions.
By Alessandro Piccolo
%}
clear all; close all;

%Set up movie
makemovie = 1;

J   = 100;       % Number of timestep t0 be used, 1000
t   = 1/J;       % Size of one time step
N1  = 20;        % Number of particles
N2  = 20;        % Other half of particle pop.
N   = N1+N2+1;   % Tot nr of particles
e   = 0.5;       % Eta is the noise parameter, whose maximum value is 2*pi
n1  = 4;         % Calculates closest n neighbors
n2  = 20;
R   = 200;       % Number of runs, 15
sprint_steps= 5; % For how many time steps can predator sprint
rest_steps  = 6; % Predator has to rest after sprint
rest_sprint_steps = sprint_steps + rest_steps + 2;
c_sprint = 1;     % Counter sprint
c_rest = 1;

%The radius of influence of a particle
L = 10;     % L is the size of the domain on which the particles can move
v = 0.5;    % v is the speed at which the particles move

n = ones(1, N-1)*n2;
index = randperm(N,N1);
n(index) = n1;
predator_index = N;
n(predator_index) = 1; % Predator is last
%alignment = runAlign2Dtask3(J, N, e, n, L, v, makemovie);

if (makemovie)
    fig = figure;
    movien = VideoWriter('FlocksPredator','Uncompressed AVI');
    open(movien);
end

% x & y coord. of the ith particle at time j (one column is a time % step)
x         = zeros(N,J+1);
x(:,1)    = L*rand(N,1);    % Define initial x coordiantes of all particles
y         = zeros(N,J+1);   % [1,N] row has x coord in column j, for time j
y(:,1)    = L*rand(N,1);    % Define initial y coordiantes of all particles
T         = zeros(N,J+1);
T(:,1)    = 2*pi*rand(N,1); % Define init direction of all particles in rad
r_kill    = 0.5;            % Distance to predator in which particles are
N1_pop    = zeros(1,J);
N2_pop    = zeros(1,J);
N1_pop(1) = N1;
N2_pop(1) = N2;


for j = 1:J % For all time steps
    N1_pop_temp = N1_pop(j);
    N2_pop_temp = N2_pop(j);
    for i = 1:N % For each particle
        % Finds n closest neigbors
        A(:,1)=((x(i,j)-x(:,j)).^2+(y(i,j)-y(:,j)).^2).^0.5;
        A(:,2)=((x(i,j)-x(:,j)-L).^2+(y(i,j)-y(:,j)).^2).^0.5;
        A(:,3)=((x(i,j)-x(:,j)).^2+(y(i,j)-y(:,j)-L).^2).^0.5;
        A(:,4)=((x(i,j)-x(:,j)+L).^2+(y(i,j)-y(:,j)).^2).^0.5;
        A(:,5)=((x(i,j)-x(:,j)).^2+(y(i,j)-y(:,j)+L).^2).^0.5;
        A(:,6)=((x(i,j)-x(:,j)+L).^2+(y(i,j)-y(:,j)+L).^2).^0.5;
        A(:,7)=((x(i,j)-x(:,j)+L).^2+(y(i,j)-y(:,j)-L).^2).^0.5;
        A(:,8)=((x(i,j)-x(:,j)-L).^2+(y(i,j)-y(:,j)+L).^2).^0.5;
        A(:,9)=((x(i,j)-x(:,j)-L).^2+(y(i,j)-y(:,j)-L).^2).^0.5;
        A_min = min(A'); % Save smallest range for each row
        if (i~=N) % Not a predator
            % Predator far away he does not matter
            if (A_min(N) > 3.5)
                % Find index of n closest particles
                B = zeros(N,1);
                A_min(N) = inf;         % Do not align with predator
                for k = 1:n(i)
                    [~,index] = min(A_min);
                    B(index) = 1;       % 1 equals to n closest neigbors
                    A_min(index) = inf; % Remove smallest value, find next
                end
                
                ss = nansum(sin(T(:,j)).*B)/nansum(B);
                sc = nansum(cos(T(:,j)).*B)/nansum(B);
                S  = atan2(ss,sc); % Avg new direction (angle) of particle
                
                T(i,j+1) = S+e*(rand-0.5);          % Adds noise to angle
            else % Predator close by uh oh! Run in opposite direction
                T(i,j+1) = -atan2(y(N,j)-y(i,j), x(N,j)-x(i,j));
                if (sqrt((y(N,j)-y(i,j))^2+(x(N,j)-x(i,j))^2) > L/2)
                    T(i,j+1) = T(i,j+1);
                end
            end
            x(i,j+1) = x(i,j)+v*cos(T(i,j+1));  % Update x-coordinate
            y(i,j+1) = y(i,j)+v*sin(T(i,j+1));  % Update y-coordinate
            x(i,j+1) = mod(x(i,j+1),L);         % Periodic BC left right
            y(i,j+1) = mod(y(i,j+1),L);         % Periodic BC top bottom
            
            % Kill particle if to close to predator
            if min(A(end,:)) < r_kill
                if (n(i) == n1)
                    N1_pop_temp = N1_pop_temp - 1;
                else
                    N2_pop_temp = N2_pop_temp - 1;
                end
                T(i,j+1) = NaN; % Remove particle
                x(i,j+1) = NaN;
                y(i,j+1) = NaN;
            end
        else % Predator
            % Look up closest target and sprint towards it
            [~,index]     = min(A_min);
            A_min(index)  = inf;
            [range,index] = min(A_min);
            % New direction (angle) of hunter
            T(i,j+1) = atan2(y(index,j)-y(N,j), x(index,j)-x(N,j));
            if (sqrt((y(index,j)-y(N,j))^2+(x(index,j)-x(N,j))^2) > L/2)
                T(i,j+1) = -T(i,j+1);
            end
            % Reset zero
            if (c_sprint + c_rest == rest_sprint_steps)
                c_sprint = 1;
                c_rest = 1;
            end
            % Sprint for sprint_steps
            if (c_sprint <= sprint_steps)
                vv = 1.6*v;
                c_sprint = c_sprint + 1;
            end
            % Rest for rest_steps after sprint
            if (c_rest <= rest_steps && c_sprint > sprint_steps)
                vv = v*0.4;
                c_rest = c_rest + 1;
            end
            x(i,j+1) = x(i,j)+vv*cos(T(i,j+1));   % Update x-coordinate
            y(i,j+1) = y(i,j)+vv*sin(T(i,j+1));   % Update y-coordinate
            x(i,j+1) = mod(x(i,j+1),L);           % Periodic BC left right
            y(i,j+1) = mod(y(i,j+1),L);           % Periodic BC top bottom
        end
        
        % Plot particles
        if (makemovie)
            if (n(i) == n1)
                if (abs(x(i,j)-x(i,j+1)) < v && abs(y(i,j)-y(i,j+1)) < v)
                    % Plots the first half of the particles in black
                    plot([x(i,j),x(i,j+1)], ...
                        [y(i,j),y(i,j+1)],'g-','markersize',4)
                    axis([0 L 0 L]);
                    hold on
                    plot(x(i,j+1) ,y(i,j+1),'g.','markersize',10)
                    xlabel('X position')
                    ylabel('Y position')
                end
            elseif (n(i) == n2)
                if (abs(x(i,j)-x(i,j+1)) < v && abs(y(i,j)-y(i,j+1)) < v)
                    % Plots the first half of the particles in black
                    plot([x(i,j),x(i,j+1)], ...
                        [y(i,j),y(i,j+1)],'b-','markersize',4)
                    axis([0 L 0 L]);
                    hold on
                    plot(x(i,j+1) ,y(i,j+1),'b.','markersize',10)
                    xlabel('X position')
                    ylabel('Y position')
                end
            else % Predator
                if (abs(x(i,j)-x(i,j+1)) < vv && abs(y(i,j)-y(i,j+1)) < vv)
                    plot(x(i,j+1) ,y(i,j+1),'r.','markersize',20)
                    xlabel('X position')
                    ylabel('Y position')
                    axis([0 L 0 L]);
                end
            end
        end
        
    end
    if (makemovie)
        hold off
        M(j) = getframe;         % Makes a movie fram from the plot
        writeVideo(movien,M(j)); % Adds this movie fram to the movie
    end
    N1_pop(j+1) = N1_pop_temp;
    N2_pop(j+1) = N2_pop_temp;
end
if (makemovie); close(movien); end;

figure(1)
plot(N1_pop, 'LineWidth', 2)
hlx = title(['Flocks and predator - n_1 = ' num2str(n1)]); 
set(hlx,'FontSize',14);
hlx = ylabel(['Surviving flock']); set(hlx,'FontSize',14);
hlx = xlabel('Time step'); set(hlx,'FontSize',14);
grid on

figure(2)
plot(N2_pop, 'LineWidth', 2)
hlx = title(['Flocks and predator - n_2 = '  num2str(n2)]); 
set(hlx,'FontSize',14);
hlx = ylabel(['Surviving flock']); set(hlx,'FontSize',14);
hlx = xlabel('Time step'); set(hlx,'FontSize',14);
grid on
