%修改人：陈颖频
%日期：2015.12.22
% 参数设置
N = 100; %粒子总数
Q = 0.005;      %过程噪声
R = 0.005;      %测量噪声
T = 300;     %测量帧数
%WorldSize = 120;    %世界大小
X = zeros(2, T);    %存储系统状态
%Z = zeros(2, T);    %存储系统的观测状态
L=load('D:\\test.txt');
Z=L(:,1:2)';
P = zeros(2, N);    %建立粒子群
PCenter = zeros(2, T);  %所有粒子的中心位置
w = zeros(N, 1);         %每个粒子的权重
err = zeros(1,T);     %误差
X(:,1)=Z(:, 1);
%X(:, 1) = [50; 20];     %初始系统状态
%Z(:, 1) = [50; 20] + wgn(2, 1, 10*log10(R));    %初始系统的观测状态
 %r = a + (b-a).*rand([m n]));
%初始化粒子群
for i = 1 : N
   P(:, i) = [39.9849+ (40.0082-39.9849).*rand; 116.3012+ (116.3249-116.3012).*rand];
   dist = norm(P(:, i)-Z(:, 1));     %与测量位置相差的距离
   w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
end
PCenter(:, 1) = sum(P, 2) / N;      %所有粒子的几何中心位置
 
%%
err(1) = norm(X(:, 1) - PCenter(:, 1));     %粒子几何中心与系统真实状态的误差
figure(1);
set(gca,'FontSize',12);
hold on
plot(X(1, 1), X(2, 1), 'r.', 'markersize',30)   %系统状态位置
%axis([0 100 0 100]);
axis([39.9800 40.0082 116.3000 116.3249]);
plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %各个粒子位置
plot(PCenter(1, 1), PCenter(2, 1), 'b.', 'markersize',25); %所有粒子的中心位置
legend('True State', 'Particles', 'The Center of Particles');
title('Initial State');
hold off


%%
%开始运动
for k = 2 : T
      
   X(:, k) = X(:, k-1) +wgn(2, 1,log10(Q));     %状态方程
   %Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程 
  
   %粒子滤波
   %预测
   for i = 1 : N
       %P(:, i) = P(:, i)  + wgn(2, 1, 10*log10(Q)); 
       P(:, i) = P(:, i)  + wgn(2, 1, log10(Q));             %只需要粒子的移动噪声参数与状态方程的噪声参数相同即可
       dist = norm(P(:, i)-Z(:, k));                                     %与测量位置相差的距离
       w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重，距离越大，权重越小，距离越小，权重越大，重要性越大
   end
%归一化权重
   wsum = sum(w);
   for i = 1 : N
       w(i) = w(i) / wsum;
   end
  
   %重采样（更新）
   for i = 1 : N
      wmax =  max(w) * 0.5;  %重要性阈值
       index = randi(N, 1);
       while(wmax > w(index))    
           wmax = wmax ;
           index = index + 1;      %忽略不重要的粒子，在本粒子滤波算法中，认为低于重要性阈值的粒子是不重要的，去掉
           if index > N
               index = 1;
           end          
       end
       P(:, i) = P(:, index);     %复制重要性粒子，这就是所谓的重采样，重新采样的意思就是重要的粒子复制，不重要的忽略
   end
  
   PCenter(:, k) = sum(P, 2) / N;      %所有粒子的中心位置，以重新采样的粒子中心作为跟踪的中心
  
   %计算误差
   err(k) = norm(X(:, k) - PCenter(:, k));     %粒子几何中心与系统真实状态的误差
  
   figure(2);
   set(gca,'FontSize',12);
   clf;
   hold on
   plot(X(1, k), X(2, k), 'r.', 'markersize',50);  %系统状态位置
   %axis([0 100 0 100]);
 axis([39.9800 40.0082 116.3000 116.3249]);
   plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %各个粒子位置
   plot(PCenter(1, k), PCenter(2, k), 'b.', 'markersize',25); %所有粒子的中心位置
   legend('True State', 'Particle', 'The Center of Particles');
   hold off
   pause(0.1);
end
 
%%
figure(3);
set(gca,'FontSize',12);
plot(X(1,:), X(2,:), 'r', Z(1,:), Z(2,:), 'g', PCenter(1,:), PCenter(2,:), 'b-');
%axis([0 100 0 100]);
axis([39.9800 40.0082 116.3000 116.3249]);
legend('True State', 'Measurement', 'Particle Filter');
xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);

%%
figure(4);
set(gca,'FontSize',12);
plot(err,'.-');
xlabel('t', 'FontSize', 20);
title('The err');                           %绘制T帧的误差曲线