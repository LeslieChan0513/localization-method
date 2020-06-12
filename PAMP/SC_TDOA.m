% SDP form 
tic;
clear all;
M = 8; % number of the sensor nodes 传感器节点的数量
% delta=0.0005
delta = 1e-6; % the penalty factor 惩罚因子
sigma_n = 10^0.5; % the measurement noise power 噪声功率的测量值
Q = sigma_n^2*(ones(M-1)+eye(M-1)); % the covariance 协方差
s_true=[400 400 -400 -400 800 800 -800 -800;
400 -400 400 -400 800 -800 800 -800];
t=[3000 10].';% the actual target position 实际的目标位置
Y_sdp = [];
Y = [];
N = 100; % number of MC runs MC的运行数量
for n = 1:N
n
s = s_true;
q = gauss_samples(zeros(M-1,1),Q,1);
dd = zeros(M-1,1); % the TDOA measurements TDOA的测量
for i = 2:M
dd(i-1) = norm(t-s_true(:,i))-norm(t-s_true(:,1))+q(i-1);
end

G = [-ones(M-1,1) eye(M-1)];
FA=G.'*inv(Q)*G;
Fb = -G.'*inv(Q)*dd;
Fc = dd'*inv(Q)*dd;
FF = [FA Fb;Fb' Fc];
% x: SDP solution SDP的结果
cvx_solver sedumi
cvx_begin sdp
cvx_quiet(1)
variable x(2);
variable z;
variable tt(M);
variable T(M,M) symmetric;
minimize (trace([T tt;tt.' 1]*FF)+ delta*sum(sum(T)))
subject to
for i=1:M
T(i,i) == [s(:,i);-1].'*[eye(2) x;x.' z]*[s(:,i);-1];
for j=1:M 
if j>i
T(i,j) >= abs([s(:,i);-1].'*[eye(2) x;x.' z]*[s(:,j);-1]); 
end
end
end
[T tt;tt.' 1] >= 0;
[eye(2) x;x.' z]>=0;
cvx_end
x
Y_sdp = [Y_sdp x]; % the set of SDP solutions SDP的解决方案
x00 = x;
xb = [-10^6 10^6;
-10^6 10^6];
x0 = [x00 xb];
[y jh] = solnp(x0,FF,s); % the local search routine 本地搜索程序
Y = [Y y]; % the set of final solutions 最终解决方案
end
for i=1:N
temp(i)=norm(Y(:,i)-t)^2;
end
RMSE=sqrt((1/N)*sum(temp)); % RMSE
toc;