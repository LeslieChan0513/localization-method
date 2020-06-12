% SDP form 
tic;
clear all;
M = 8; % number of the sensor nodes �������ڵ������
% delta=0.0005
delta = 1e-6; % the penalty factor �ͷ�����
sigma_n = 10^0.5; % the measurement noise power �������ʵĲ���ֵ
Q = sigma_n^2*(ones(M-1)+eye(M-1)); % the covariance Э����
s_true=[400 400 -400 -400 800 800 -800 -800;
400 -400 400 -400 800 -800 800 -800];
t=[3000 10].';% the actual target position ʵ�ʵ�Ŀ��λ��
Y_sdp = [];
Y = [];
N = 100; % number of MC runs MC����������
for n = 1:N
n
s = s_true;
q = gauss_samples(zeros(M-1,1),Q,1);
dd = zeros(M-1,1); % the TDOA measurements TDOA�Ĳ���
for i = 2:M
dd(i-1) = norm(t-s_true(:,i))-norm(t-s_true(:,1))+q(i-1);
end

G = [-ones(M-1,1) eye(M-1)];
FA=G.'*inv(Q)*G;
Fb = -G.'*inv(Q)*dd;
Fc = dd'*inv(Q)*dd;
FF = [FA Fb;Fb' Fc];
% x: SDP solution SDP�Ľ��
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
Y_sdp = [Y_sdp x]; % the set of SDP solutions SDP�Ľ������
x00 = x;
xb = [-10^6 10^6;
-10^6 10^6];
x0 = [x00 xb];
[y jh] = solnp(x0,FF,s); % the local search routine ������������
Y = [Y y]; % the set of final solutions ���ս������
end
for i=1:N
temp(i)=norm(Y(:,i)-t)^2;
end
RMSE=sqrt((1/N)*sum(temp)); % RMSE
toc;