%%
%how tw looks like
% close all
% x = [1;0;0];
% M = [0.9, 0,0;0.1, 0.9,0;0,0.1,0.9];
% for i = 2:100
% x(:,i) = M*x(:,i-1);
% end
% figure
% plot(x')
% 
% figure
% scatter3(x(1,:),x(2,:),x(3,:))
% 
% 
% u,s,v = svd(x);


%%

% 
% close all
% x = [1;0.5;-2];
% %M = randn(3,3);
% M = [1,0,0;0,0.99,0;0,0,1.01];%[1,0.1,-0.5;0.1,1.2,0.4;-0.5,0.4,0.8]*0.63;
% for i = 2:100
% x(:,i) = M*x(:,i-1);
% end
% figure
% plot(x')
% 
% figure
% scatter3(x(1,:),x(2,:),x(3,:))


%u,s,v = svd(x);

%%
close all

N = 3;
 delta = 0.002;%0.125/2;
 f = 10.0;
 fs = 1000;
 v = 0.5;
 nT = 92;
 
 x = zeros(N, nT);
 cum_dist = 0.;
 
 for i = 1:N
     for t = 0:nT-1
         x(i, t+1) = cum_dist/v+t/fs;
     end
     cum_dist = cum_dist+delta;
 end
 
 s = sin(2*pi*(f)*x);
 
 [u,ss,v] = svd(s,'econ');
 %s = u'*s;
 
 figure
 plot(s(1,:)', 'LineWidth',2, 'Color','b')
 hold on
 plot(s(2,:)', 'LineWidth',2,  'Color','m')
 hold on
 plot(s(3,:)', 'LineWidth',2,  'Color','r')
 hold on
 legend('Sensor 3', 'Sensor 2', 'Sensor 1')
xlabel('time')
ylabel('magnitude')
c = linspace(1,10,nT);
figure
 scatter3(s(1,:),s(2,:),s(3,:),30, c,'filled')
 
 xlabel('Sensor 1')
 ylabel('Sensor 2')
 zlabel('Sensor 3')
 
 ds = s(:,2:end)-s(:,1:end-1);
 
 M = ds*pinv(s(:,1:end-1));
 
 err_full = norm(ds-M*s(:,1:end-1)-mean(ds,2),'fro');
 
 M_sym = travelling_matrix(s(:,1:end-1)',ds')';
 
 err_sym = norm(ds-M_sym*s(:,1:end-1),'fro');
 
 p2 = N^2;
 p1 = N*(N-1)/2+N;
 F_value = ((-err_full^2+err_sym^2)/(p2-p1))/((err_sym)^2/(nT-p2));

%%
% close all
% 
% N = 2;
% delta = 0.125/2;
% f = 10.0;
% fs = 1000;
% v = 0.5;
% nT = 100;
% 
% x = zeros(N, nT);
% cum_dist = 0.;
% 
% for i = 1:N
%     for t = 0:nT-1
%         x(i, t+1) = cum_dist/v+t/fs;
%     end
%     cum_dist = cum_dist+delta;
% end
% 
% s = sin(2*pi*(f)*x);
% 
% %[u,ss,v] = svd(s,'econ');
% %s = u'*s;
% 
% 
% 
% ds = s(:,2:end)-s(:,1:end-1);
% 
% M = ds*pinv(s(:,1:end-1));
% 
% 
% 
% t = linspace(-1,1,nT);
% figure
% %plot3(x(1,:),x(2,:),zeros(1,nT),'LineWidth',2)%,'Color',t,'interp')
% hold on
% surface([s(1,:);s(1,:)],[s(2,:);s(2,:)],[zeros(1,nT);zeros(1,nT)],[t ;t],...
%     'facecol','no','edgecol','interp','LineWidth', 2);
% hold on
% plot3([0,0],[-10,10],[0,0],'--','Color','k')
% hold on
% plot3([-10,10],[0,0],[0,0],'--','Color','k')
% 
% colormap default
% ylim([-1.2,1.2])
% xlim([-1.2,1.2])
% xlabel('x')
% ylabel('y')
% title('skew-symmetric M')









%%
% 
% close all
% 
% M = [-0.1,-0.15;-0.15,0.15]; % M = [-0.005,0.1;0,-0.005]%
% nT = 20;
% x = [];
% x(:,1) = [2;1];
% for i = 2:nT
%     x(:,i) = x(:,i-1) + M*x(:,i-1);
% end
% 
% t = linspace(-1,1,nT);
% figure
% %plot3(x(1,:),x(2,:),zeros(1,nT),'LineWidth',2)%,'Color',t,'interp')
% hold on
% surface([x(1,:);x(1,:)],[x(2,:);x(2,:)],[zeros(1,nT);zeros(1,nT)],[t ;t],...
%     'facecol','no','edgecol','interp','LineWidth', 2);
% hold on
% plot3([0,0],[-10,10],[0,0],'--','Color','k')
% hold on
% plot3([-10,10],[0,0],[0,0],'--','Color','k')
% 
% colormap default
% ylim([-0.5,2.5])
% xlim([-1.5,2.5])
% xlabel('x')
% ylabel('y')
% title('symmetric M')
% 


%%
% close all
% N = 500;
% X = zeros(N);
% 
% x = linspace(0.5,1,N);
% for i = 1:N
%     X(i,:) = x;
%     x = x-0.5/(N-1);
% end
% figure
% imshow(X)
% colormap 'jet'


%%
% close all
% N = 500;
% X = zeros(N);
% 
% x = linspace(0,1,N);
% for i = 1:N
%     X(i,:) = x;
%     x = x-1/(N-1);
% end
% figure
% imshow(abs(X))
% colormap 'jet'



