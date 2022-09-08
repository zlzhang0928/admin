clear all;
%矩阵A1没有模型误差，名义模型系数
A1=[0.9802 0.0196;0 0.9802];%+[0.0198;0.0]*(-0.8508)*[0.0 5.0];  %改变0.0196的值，相对的改变了模型误差的影响
%矩阵A2没有模型误差，名义模型系数
A2=[-0.021 0.006;0 -0.021];%+[0.0198;0.0]*(-0.8508)*[0.0 5.0]; %改变0.006的值，相对的改变了模型误差的影响%矩阵B1没有模型误差，名义模型系数%矩阵B1没有模型误差，名义模型系数
A3=[0.9902 0.0196;0 0.9902];%+[0.0198;0.0]*(-0.8508)*[0.0 5.0];  %改变0.0196的值，相对的改变了模型误差的影响
%矩阵A2没有模型误差，名义模型系数
A4=[-0.024 0.006;0 -0.024];
B1=[0.0196 0;0 0.0196];%+[0.0198;0.0]*(-0.8508)*[0.0 5.0]; %改变0.00的值，相对的改变了模型误差的影响
BB1=[0.0196 0;0 0.0196];
B2=[1 0;0 1];
C=[1 -1];
n=2;
m=1;
d=2;
bar_A =A1 ; %%创造状态矩阵Ai
bar_A2=A3 ; %%创造状态矩阵Ai
bar_B1=B1;      %%创造状态矩阵B1i代表状态输入
bar_BB1=BB1;
bar_B2=B2;           %%创造状态矩阵B2i代表过程误差状态矩阵               
bar_C=C ;                %%创造状态矩阵Ci代表测量输出状态矩阵 
Q=[1.9608 0.0195;0.0195 1.9605];
bar_Q=Q ;
R=1;
PI0=eye(2);
epsilon12=cell(500,1);
for i=1:500
      while 1
     epsilon12{i}=normrnd(0,1);    %正态分布  %修改随机产生模型误差的方差
  if  abs(epsilon12{i})<=1       %不同实验随机产生一个值
     break;
  end
  end
       
end
z=cell(500,1);%先前变量Z(k)
for i=1:500
    z{i}=cell(1003,1);
    for t=1:1003
        z{i}{t}=0;
    end
end
u=cell(500,1);%隶属函数
for i=1:500
    u{i}=cell(1003,1);
    for t=1:1003
        u{i}{t}=cell(2,1);
       
            u{i}{t}{1}=0.8;
            u{i}{t}{2}=0.2;
    end
end
xx=cell(500,1);%X(k)真实状态
for i=1:500
    xx{i}=cell(1003,1);
    for t=1:1003
        xx{i}{t}=[0;0];
    end
end
% es=[1.0;0.1];
                           %产生500次实验数据
% ---------------------------------------------------------------------------%          
    X1           = zeros(2,1003,500);
    Y1           = zeros(1,1003,500);
    X2           = zeros(2,1003,500);
    Y2           = zeros(1,1003,500);
Real_A = cell(500,1);
Real_B = cell(500,1);
Real_A2 = cell(500,1);

for i=1:500 %500次实验
for t=1:1003
   esk=cell2mat(epsilon12(i,1));
if t == 1
X1(:,t,i)=[0;0];
Y1(:,t,i)=[ 0 ];%定义初始输出
X2(:,t,i)=[0;0];
Y2(:,t,i)=[ 0 ];%定义初始输出
else
Real_A1_Temp = A1+[0.0298;0.0]*esk*[0.0 5.0]; 
Real_A2_Temp = A2+[0.0298;0.0]*esk*[0.0 5.0]; 
Real_A3_Temp = A3+[0.0298;0.0]*esk*[0.0 5.0]; 
Real_A4_Temp = A4+[0.0298;0.0]*esk*[0.0 5.0];
Real_B1_Temp = B1+[0.0298;0.0]*esk*[0.0 5.0];
Real_BB1_Temp = BB1+[0.0298;0.0]*esk*[0.0 5.0];
Real_A_Temp=Real_A1_Temp ;
Real_B_Temp=Real_B1_Temp ;
Real_BB_Temp=Real_BB1_Temp ; 

Real_A2_Temp=Real_A3_Temp ;

Real_A{i}{t-1}=Real_A_Temp;
Real_A2{i}{t-1}=Real_A2_Temp;
Real_B{i}{t-1}=Real_B_Temp;

X1(:,t,i)=Real_A_Temp* X1(:,t-1,i)+u{i}{t}{1}*bar_B2*Q^(1/2)*[normrnd(0,1);normrnd(0,1)]+u{i}{t}{1}*Real_B_Temp*[1.0;0.1];
Y1(:,t,i)=bar_C* X1(:,t,i)+u{i}{t}{1}*normrnd(0,1);

X2(:,t,i)=Real_A2_Temp* X2(:,t-1,i)+u{i}{t}{2}*bar_B2*Q^(1/2)*[normrnd(0,1);normrnd(0,1)]+u{i}{t}{2}*Real_BB_Temp*[1.0;0.1];
Y2(:,t,i)=bar_C* X2(:,t,i)+u{i}{t}{1}*normrnd(0,1);
xx{i}{t}=X1(1:2,t,i)+X2(1:2,t,i);
z{i}{t}=(xx{i}{t}(1)*3.1416/180-0.0098*xx{i}{t}(2)*3.1416/180);
u{i}{t}{1}=1-1/(1+exp(3.1416/2-z{i}{t}));
u{i}{t}{2}=1-u{i}{t}{1};

end
%每个时刻的X为6个元素的列向量，每个时刻的Y为1个元素的向量
end
end
 X_hat = zeros(2,1001,500);
 X_hat2 = zeros(2,1001,500);
 fprintf('轨迹数据已生成\n');
                      

                             %对500次实验数据分别运用灵敏度惩罚算法
                                  %  gamma =0.0305
% ---------------------------------------------------------------------------%
S=[0,0;0,0.099];
T1=[0,0;0,0.099];      
T2=[0,0;0,0];
    gamma1 =(1-0.035)/0.035; %权重系数
    Pi1i= zeros(2);
    for i=1:500
    for t=1:1001
        if t == 1
            Pi1i        =(PI0^(-1) + bar_C'*R^(-1)*bar_C)^(-1);
            Xili_hat    = Pi1i*bar_C'*R^(-1)*Y1(:,t,i); 
        else
            Pi1i_hat	=(Pi1i^(-1) + gamma1*S'*S)^(-1);
            T2_hat      =T2-gamma1*S*Pi1i_hat*S'*T2;
            Q_hat       =(Q^(-1)+gamma1*T2'*(eye(2)+gamma1*S*Pi1i*S')^(-1)*T2)^(-1);
            B2_hat      = u{i}{t}{1}*bar_B2-gamma1*bar_A*Pi1i_hat*S'*T2;
            A_hat       =(bar_A-gamma1*B2_hat*Q_hat*T2'*S)*(eye(2)-gamma1*Pi1i_hat*S'*S);
            B1_hat      =u{i}{t}{1}*bar_B1-gamma1*(bar_A*Pi1i_hat*S'+B2_hat*Q_hat*T2_hat')*T1; 
            
            Pit1li      = bar_A*Pi1i_hat*bar_A' + B2_hat*Q_hat*B2_hat';
            Reit1       = R + bar_C*Pit1li*bar_C';
            Pit1lit1	= Pit1li - Pit1li*bar_C'*Reit1^(-1)*bar_C*Pit1li;
            xit1lit1_hat= B1_hat*[1.0;0.1]+A_hat*Xili_hat+ Pit1lit1*bar_C'*R^(-1)*(Y1(:,t,i)-bar_C*(B1_hat*[1.0;0.1]+ A_hat*Xili_hat));

            Xili_hat	= xit1lit1_hat;
            Pi1i        = Pit1lit1;
        end
            X_hat(:,t,i)= Xili_hat;
    end
    end
    
    for i=1:500
    for t=1:1001
        if t == 1
            Pi1i        =(PI0^(-1) + bar_C'*R^(-1)*bar_C)^(-1);
            Xili_hat    = Pi1i*bar_C'*R^(-1)*Y2(:,t,i); 
        else
            Pi1i_hat	=(Pi1i^(-1) + gamma1*S'*S)^(-1);
            T2_hat      =T2-gamma1*S*Pi1i_hat*S'*T2;
            Q_hat       =(Q^(-1)+gamma1*T2'*(eye(2)+gamma1*S*Pi1i*S')^(-1)*T2)^(-1);
            B2_hat      = u{i}{t}{2}*bar_B2-gamma1*bar_A2*Pi1i_hat*S'*T2;
            A_hat       =(bar_A2-gamma1*B2_hat*Q_hat*T2'*S)*(eye(2)-gamma1*Pi1i_hat*S'*S);
            B1_hat      =u{i}{t}{2}*bar_B1-gamma1*(bar_A2*Pi1i_hat*S'+B2_hat*Q_hat*T2_hat')*T1; 
            
            Pit1li      = bar_A2*Pi1i_hat*bar_A2' + B2_hat*Q_hat*B2_hat';
            Reit1       = R + bar_C*Pit1li*bar_C';
            Pit1lit1	= Pit1li - Pit1li*bar_C'*Reit1^(-1)*bar_C*Pit1li;
            xit1lit1_hat= B1_hat*[1.0;0.1]+A_hat*Xili_hat+ Pit1lit1*bar_C'*R^(-1)*(Y2(:,t,i)-bar_C*(B1_hat*[1.0;0.1]+ A_hat*Xili_hat));

            Xili_hat	= xit1lit1_hat;
            Pi1i        = Pit1lit1;
        end
            X_hat2(:,t,i)= Xili_hat;
    end
    end
Method1_Deviation_SUM       = zeros(1,1000);
Method1_Deviation_SUM2      = zeros(1,1000);
Method1_Standard_Deviation1	= zeros(1,1000);
Method1_Standard_Deviation2	= zeros(1,1000);
Method1_Standard_Deviation	= zeros(1,1000);
Method1_STD_dB              = zeros(1,1000);
for t=1:1000
     for i=1:500
       Method1_Deviation_SUM(t) = Method1_Deviation_SUM(t) + ...
                                     ( X_hat(1,t+1,i) - X1(1,t+1,i) )^2+( X_hat(2,t+1,i) - X1(2,t+1,i) )^2;
     
       Method1_Deviation_SUM2(t) = Method1_Deviation_SUM2(t) + ...
                                     ( X_hat2(1,t+1,i) - X2(1,t+1,i) )^2+( X_hat2(2,t+1,i) - X2(2,t+1,i) )^2;
     end
     Method1_Standard_Deviation1(t) = Method1_Deviation_SUM(t) / 500;
     Method1_Standard_Deviation2(t) = Method1_Deviation_SUM2(t) / 500;
     Method1_Standard_Deviation(t) = (Method1_Standard_Deviation1(t)^2+Method1_Standard_Deviation2(t)^2)^(1/2);
     Method1_STD_dB(t)             = 10 * log10( Method1_Standard_Deviation(t) ); 
end
n1=Method1_STD_dB(500);
fprintf('预测模型1+鲁棒状态估计方法已完成\n');

 %对500次实验数据分别运用灵敏度惩罚算法
                                  %  gamma =0.0305
% ---------------------------------------------------------------------------%
S=[0,0;0,0.099];
T1=[0,0;0,0.099];      
T2=[0,0;0,0];
    gamma2 =(1-0.289)/0.289; %权重系数
    Pi1i= zeros(2);
    for i=1:500
    for t=1:1001
        if t == 1
            Pi1i        =(PI0^(-1) + bar_C'*R^(-1)*bar_C)^(-1);
            Xili_hat    = Pi1i*bar_C'*R^(-1)*Y1(:,t,i); 
        else
            Pi1i_hat	=(Pi1i^(-1) + gamma2*S'*S)^(-1);
            T2_hat      =T2-gamma2*S*Pi1i_hat*S'*T2;
            Q_hat       =(Q^(-1)+gamma2*T2'*(eye(2)+gamma2*S*Pi1i*S')^(-1)*T2)^(-1);
            B2_hat      = u{i}{t}{1}*bar_B2-gamma2*bar_A*Pi1i_hat*S'*T2;
            A_hat       =(bar_A-gamma2*B2_hat*Q_hat*T2'*S)*(eye(2)-gamma2*Pi1i_hat*S'*S);
            B1_hat      =u{i}{t}{1}*bar_B1-gamma2*(bar_A*Pi1i_hat*S'+B2_hat*Q_hat*T2_hat')*T1; 
            
            Pit1li      = bar_A*Pi1i_hat*bar_A' + B2_hat*Q_hat*B2_hat';
            Reit1       = R + bar_C*Pit1li*bar_C';
            Pit1lit1	= Pit1li - Pit1li*bar_C'*Reit1^(-1)*bar_C*Pit1li;
            xit1lit1_hat= B1_hat*[1.0;0.1]+A_hat*Xili_hat+ Pit1lit1*bar_C'*R^(-1)*(Y1(:,t,i)-bar_C*(B1_hat*[1.0;0.1]+ A_hat*Xili_hat));

            Xili_hat	= xit1lit1_hat;
            Pi1i        = Pit1lit1;
        end
            X_hat(:,t,i)= Xili_hat;
    end
    end
    
    for i=1:500
    for t=1:1001
        if t == 1
            Pi1i        =(PI0^(-1) + bar_C'*R^(-1)*bar_C)^(-1);
            Xili_hat    = Pi1i*bar_C'*R^(-1)*Y2(:,t,i); 
        else
            Pi1i_hat	=(Pi1i^(-1) + gamma2*S'*S)^(-1);
            T2_hat      =T2-gamma2*S*Pi1i_hat*S'*T2;
            Q_hat       =(Q^(-1)+gamma2*T2'*(eye(2)+gamma2*S*Pi1i*S')^(-1)*T2)^(-1);
            B2_hat      = u{i}{t}{2}*bar_B2-gamma2*bar_A2*Pi1i_hat*S'*T2;
            A_hat       =(bar_A2-gamma2*B2_hat*Q_hat*T2'*S)*(eye(2)-gamma2*Pi1i_hat*S'*S);
            B1_hat      =u{i}{t}{2}*bar_B1-gamma2*(bar_A2*Pi1i_hat*S'+B2_hat*Q_hat*T2_hat')*T1; 
            
            Pit1li      = bar_A2*Pi1i_hat*bar_A2' + B2_hat*Q_hat*B2_hat';
            Reit1       = R + bar_C*Pit1li*bar_C';
            Pit1lit1	= Pit1li - Pit1li*bar_C'*Reit1^(-1)*bar_C*Pit1li;
            xit1lit1_hat= B1_hat*[1.0;0.1]+A_hat*Xili_hat+ Pit1lit1*bar_C'*R^(-1)*(Y2(:,t,i)-bar_C*(B1_hat*[1.0;0.1]+ A_hat*Xili_hat));

            Xili_hat	= xit1lit1_hat;
            Pi1i        = Pit1lit1;
        end
            X_hat2(:,t,i)= Xili_hat;
    end
    end
Method2_Deviation_SUM       = zeros(1,1000);
Method2_Deviation_SUM2       = zeros(1,1000);
Method2_Standard_Deviation1	= zeros(1,1000);
Method2_Standard_Deviation2	= zeros(1,1000);
Method2_Standard_Deviation	= zeros(1,1000);
Method2_STD_dB              = zeros(1,1000);
for t=1:1000
     for i=1:500
       Method2_Deviation_SUM(t) = Method2_Deviation_SUM(t) + ...
                                     ( X_hat(1,t+1,i) - X1(1,t+1,i) )^2+( X_hat(2,t+1,i) - X1(2,t+1,i) )^2;
     
       Method2_Deviation_SUM2(t) = Method2_Deviation_SUM2(t) + ...
                                     ( X_hat2(1,t+1,i) - X2(1,t+1,i) )^2+( X_hat2(2,t+1,i) - X2(2,t+1,i) )^2;
     end
     Method2_Standard_Deviation1(t) = Method2_Deviation_SUM(t) / 500;
     Method2_Standard_Deviation2(t) = Method2_Deviation_SUM2(t) / 500;
     Method2_Standard_Deviation(t) = (Method2_Standard_Deviation1(t)^2+Method2_Standard_Deviation2(t)^2)^(1/2);
     Method2_STD_dB(t)             = 10 * log10( Method2_Standard_Deviation(t) ); 
end
n2=Method2_STD_dB(500);
fprintf('预测模型2+鲁棒状态估计方法已完成\n');
%对500次实验数据分别运用灵敏度惩罚算法
                                  %  gamma =0.0305
% ---------------------------------------------------------------------------%
S=[0,0;0,0.099];
T1=[0,0;0,0.099];      
T2=[0,0;0,0];
gamma3 =(1-0.75)/0.75; %权重系数
    Pi1i= zeros(2);
    for i=1:500
    for t=1:1001
        if t == 1
            Pi1i        =(PI0^(-1) + bar_C'*R^(-1)*bar_C)^(-1);
            Xili_hat    = Pi1i*bar_C'*R^(-1)*Y1(:,t,i); 
        else
            Pi1i_hat	=(Pi1i^(-1) + gamma3*S'*S)^(-1);
            T2_hat      =T2-gamma3*S*Pi1i_hat*S'*T2;
            Q_hat       =(Q^(-1)+gamma3*T2'*(eye(2)+gamma3*S*Pi1i*S')^(-1)*T2)^(-1);
            B2_hat      = u{i}{t}{1}*bar_B2-gamma3*bar_A*Pi1i_hat*S'*T2;
            A_hat       =(bar_A-gamma3*B2_hat*Q_hat*T2'*S)*(eye(2)-gamma3*Pi1i_hat*S'*S);
            B1_hat      =u{i}{t}{1}*bar_B1-gamma3*(bar_A*Pi1i_hat*S'+B2_hat*Q_hat*T2_hat')*T1; 
            
            Pit1li      = bar_A*Pi1i_hat*bar_A' + B2_hat*Q_hat*B2_hat';
            Reit1       = R + bar_C*Pit1li*bar_C';
            Pit1lit1	= Pit1li - Pit1li*bar_C'*Reit1^(-1)*bar_C*Pit1li;
            xit1lit1_hat= B1_hat*[1.0;0.1]+A_hat*Xili_hat+ Pit1lit1*bar_C'*R^(-1)*(Y1(:,t,i)-bar_C*(B1_hat*[1.0;0.1]+ A_hat*Xili_hat));

            Xili_hat	= xit1lit1_hat;
            Pi1i        = Pit1lit1;
        end
            X_hat(:,t,i)= Xili_hat;
    end
    end
    
    for i=1:500
    for t=1:1001
        if t == 1
            Pi1i        =(PI0^(-1) + bar_C'*R^(-1)*bar_C)^(-1);
            Xili_hat    = Pi1i*bar_C'*R^(-1)*Y2(:,t,i); 
        else
            Pi1i_hat	=(Pi1i^(-1) + gamma3*S'*S)^(-1);
            T2_hat      =T2-gamma3*S*Pi1i_hat*S'*T2;
            Q_hat       =(Q^(-1)+gamma3*T2'*(eye(2)+gamma3*S*Pi1i*S')^(-1)*T2)^(-1);
            B2_hat      = u{i}{t}{2}*bar_B2-gamma3*bar_A2*Pi1i_hat*S'*T2;
            A_hat       =(bar_A2-gamma3*B2_hat*Q_hat*T2'*S)*(eye(2)-gamma3*Pi1i_hat*S'*S);
            B1_hat      =u{i}{t}{2}*bar_B1-gamma3*(bar_A2*Pi1i_hat*S'+B2_hat*Q_hat*T2_hat')*T1; 
            
            Pit1li      = bar_A2*Pi1i_hat*bar_A2' + B2_hat*Q_hat*B2_hat';
            Reit1       = R + bar_C*Pit1li*bar_C';
            Pit1lit1	= Pit1li - Pit1li*bar_C'*Reit1^(-1)*bar_C*Pit1li;
            xit1lit1_hat= B1_hat*[1.0;0.1]+A_hat*Xili_hat+ Pit1lit1*bar_C'*R^(-1)*(Y2(:,t,i)-bar_C*(B1_hat*[1.0;0.1]+ A_hat*Xili_hat));

            Xili_hat	= xit1lit1_hat;
            Pi1i        = Pit1lit1;
        end
            X_hat2(:,t,i)= Xili_hat;
    end
    end
Method3_Deviation_SUM       = zeros(1,1000);
Method3_Deviation_SUM2       = zeros(1,1000);
Method3_Standard_Deviation1	= zeros(1,1000);
Method3_Standard_Deviation2	= zeros(1,1000);
Method3_Standard_Deviation	= zeros(1,1000);
Method3_STD_dB              = zeros(1,1000);
for t=1:1000
     for i=1:500
       Method3_Deviation_SUM(t) = Method3_Deviation_SUM(t) + ...
                                     ( X_hat(1,t+1,i) - X1(1,t+1,i) )^2+( X_hat(2,t+1,i) - X1(2,t+1,i) )^2;
     
       Method3_Deviation_SUM2(t) = Method3_Deviation_SUM2(t) + ...
                                     ( X_hat2(1,t+1,i) - X2(1,t+1,i) )^2+( X_hat2(2,t+1,i) - X2(2,t+1,i) )^2;
     end
     Method3_Standard_Deviation1(t) = Method3_Deviation_SUM(t) / 500;
     Method3_Standard_Deviation2(t) = Method3_Deviation_SUM2(t) / 500;
     Method3_Standard_Deviation(t) = (Method3_Standard_Deviation1(t)^2+Method3_Standard_Deviation2(t)^2)^(1/2);
     Method3_STD_dB(t)             = 10 * log10( Method3_Standard_Deviation(t) ); 
end
n3=Method3_STD_dB(500);
fprintf('预测模型3+鲁棒状态估计方法已完成\n');

%对500次实验数据分别运用灵敏度惩罚算法
                                  %  gamma =0.0305
% ---------------------------------------------------------------------------%
S=[0,0;0,0.099];
T1=[0,0;0,0.099];      
T2=[0,0;0,0];
    gamma4 =(1-0.949)/0.949; %权重系数
    Pi1i= zeros(2);
    for i=1:500
    for t=1:1001
        if t == 1
            Pi1i        =(PI0^(-1) + bar_C'*R^(-1)*bar_C)^(-1);
            Xili_hat    = Pi1i*bar_C'*R^(-1)*Y1(:,t,i); 
        else
            Pi1i_hat	=(Pi1i^(-1) + gamma4*S'*S)^(-1);
            T2_hat      =T2-gamma4*S*Pi1i_hat*S'*T2;
            Q_hat       =(Q^(-1)+gamma4*T2'*(eye(2)+gamma4*S*Pi1i*S')^(-1)*T2)^(-1);
            B2_hat      = u{i}{t}{1}*bar_B2-gamma4*bar_A*Pi1i_hat*S'*T2;
            A_hat       =(bar_A-gamma4*B2_hat*Q_hat*T2'*S)*(eye(2)-gamma4*Pi1i_hat*S'*S);
            B1_hat      =u{i}{t}{1}*bar_B1-gamma4*(bar_A*Pi1i_hat*S'+B2_hat*Q_hat*T2_hat')*T1; 
            
            Pit1li      = bar_A*Pi1i_hat*bar_A' + B2_hat*Q_hat*B2_hat';
            Reit1       = R + bar_C*Pit1li*bar_C';
            Pit1lit1	= Pit1li - Pit1li*bar_C'*Reit1^(-1)*bar_C*Pit1li;
            xit1lit1_hat= B1_hat*[1.0;0.1]+A_hat*Xili_hat+ Pit1lit1*bar_C'*R^(-1)*(Y1(:,t,i)-bar_C*(B1_hat*[1.0;0.1]+ A_hat*Xili_hat));

            Xili_hat	= xit1lit1_hat;
            Pi1i        = Pit1lit1;
        end
            X_hat(:,t,i)= Xili_hat;
    end
    end
    
    for i=1:500
    for t=1:1001
        if t == 1
            Pi1i        =(PI0^(-1) + bar_C'*R^(-1)*bar_C)^(-1);
            Xili_hat    = Pi1i*bar_C'*R^(-1)*Y2(:,t,i); 
        else
            Pi1i_hat	=(Pi1i^(-1) + gamma4*S'*S)^(-1);
            T2_hat      =T2-gamma4*S*Pi1i_hat*S'*T2;
            Q_hat       =(Q^(-1)+gamma4*T2'*(eye(2)+gamma4*S*Pi1i*S')^(-1)*T2)^(-1);
            B2_hat      = u{i}{t}{2}*bar_B2-gamma4*bar_A2*Pi1i_hat*S'*T2;
            A_hat       =(bar_A2-gamma4*B2_hat*Q_hat*T2'*S)*(eye(2)-gamma4*Pi1i_hat*S'*S);
            B1_hat      =u{i}{t}{2}*bar_B1-gamma4*(bar_A2*Pi1i_hat*S'+B2_hat*Q_hat*T2_hat')*T1; 
            
            Pit1li      = bar_A2*Pi1i_hat*bar_A2' + B2_hat*Q_hat*B2_hat';
            Reit1       = R + bar_C*Pit1li*bar_C';
            Pit1lit1	= Pit1li - Pit1li*bar_C'*Reit1^(-1)*bar_C*Pit1li;
            xit1lit1_hat= B1_hat*[1.0;0.1]+A_hat*Xili_hat+ Pit1lit1*bar_C'*R^(-1)*(Y2(:,t,i)-bar_C*(B1_hat*[1.0;0.1]+ A_hat*Xili_hat));

            Xili_hat	= xit1lit1_hat;
            Pi1i        = Pit1lit1;
        end
            X_hat2(:,t,i)= Xili_hat;
    end
    end
Method4_Deviation_SUM       = zeros(1,1000);
Method4_Deviation_SUM2       = zeros(1,1000);
Method4_Standard_Deviation1	= zeros(1,1000);
Method4_Standard_Deviation2	= zeros(1,1000);
Method4_Standard_Deviation	= zeros(1,1000);
Method4_STD_dB              = zeros(1,1000);
for t=1:1000
     for i=1:500
       Method4_Deviation_SUM(t) = Method4_Deviation_SUM(t) + ...
                                     ( X_hat(1,t+1,i) - X1(1,t+1,i) )^2+( X_hat(2,t+1,i) - X1(2,t+1,i) )^2;
     
       Method4_Deviation_SUM2(t) = Method4_Deviation_SUM2(t) + ...
                                     ( X_hat2(1,t+1,i) - X2(1,t+1,i) )^2+( X_hat2(2,t+1,i) - X2(2,t+1,i) )^2;
     end
     Method4_Standard_Deviation1(t) = Method4_Deviation_SUM(t) / 500;
     Method4_Standard_Deviation2(t) = Method4_Deviation_SUM2(t) / 500;
     Method4_Standard_Deviation(t) = (Method4_Standard_Deviation1(t)^2+Method4_Standard_Deviation2(t)^2)^(1/2);
     Method4_STD_dB(t)             = 10 * log10( Method4_Standard_Deviation(t) ); 
end
n4=Method4_STD_dB(500);
fprintf('预测模型4+鲁棒状态估计方法已完成\n');

 
    
  
                          % 基于名义参数的Kalman 
% -------------------------------------------------------------------------- %
xn	= zeros(2, 1);                   % x-negetive（x-）先验估计
xp	= zeros(2, 1);                   % x-positive（x+）后验估计
pn	= zeros(2);                   % p-negetive（p-）先验估计方差
pp	= zeros(2);                   % p-positive（p+）后验估计方差
g	= zeros(2, 1);                   % 稳定增益矩阵

for i=1:500
    for t=1:1001
            if t == 1
                pp = ( PI0^(-1) + bar_C' * R^(-1) * bar_C )^(-1);
                xp = pp * bar_C' * R^(-1) * Y1(:,t,i);
                
            else
                xn = bar_A * xp+u{i}{t}{1}*bar_B1*[1;0.1];
                pn = bar_A * pp * bar_A' + bar_Q;
                g  = (pn * bar_C') * inv( bar_C * pn * bar_C' + R);
                xp = xn + g * (Y1(:,t,i) - bar_C * xn);
                pp = (eye(2) - g * bar_C) * pn;   
                
                
            end
            X_hat(:,t,i)	= xp';
         
    end
end
for i=1:500
    for t=1:1001
            if t == 1
                pp = ( PI0^(-1) + bar_C' * R^(-1) * bar_C )^(-1);
                xp = pp * bar_C' * R^(-1) * Y2(:,t,i);
                
            else
                xn = bar_A2 * xp+u{i}{t}{2}*bar_BB1*[1;0.1];
                pn = bar_A2 * pp * bar_A2' + bar_Q;
                g  = (pn * bar_C') * inv( bar_C * pn * bar_C' + R);
                xp = xn + g * (Y2(:,t,i) - bar_C * xn);
                pp = (eye(2) - g * bar_C) * pn;   
                
                
            end
            X_hat2(:,t,i)	= xp';
         
    end
end
Method5_Deviation_SUM       = zeros(1,1000);
Method5_Deviation_SUM2      = zeros(1,1000);
Method5_Standard_Deviation1	= zeros(1,1000);
Method5_Standard_Deviation2	= zeros(1,1000);
Method5_Standard_Deviation	= zeros(1,1000);
Method5_STD_dB              = zeros(1,1000);
for t=1:1000
     for i=1:500
       Method5_Deviation_SUM(t) = Method5_Deviation_SUM(t) + ...
                                     ( X_hat(1,t+1,i) - X1(1,t+1,i) )^2+( X_hat(2,t+1,i) - X1(2,t+1,i) )^2;
     
       Method5_Deviation_SUM2(t) = Method5_Deviation_SUM2(t) + ...
                                     ( X_hat2(1,t+1,i) - X2(1,t+1,i) )^2+( X_hat2(2,t+1,i) - X2(2,t+1,i) )^2;
     end
     Method5_Standard_Deviation1(t) = Method5_Deviation_SUM(t) / 500;
     Method5_Standard_Deviation2(t) = Method5_Deviation_SUM2(t) / 500;
     Method5_Standard_Deviation(t) = (Method5_Standard_Deviation1(t)^2+Method5_Standard_Deviation2(t)^2)^(1/2);
     Method5_STD_dB(t)             = 10 * log10( Method5_Standard_Deviation(t) ); 
end
n5=Method5_STD_dB(500);

  % 基于实际参数的Kalman 
% -------------------------------------------------------------------------- %
xn	= zeros(2, 1);                   % x-negetive（x-）先验估计
xp	= zeros(2, 1);                   % x-positive（x+）后验估计
pn	= zeros(2);                   % p-negetive（p-）先验估计方差
pp	= zeros(2);                   % p-positive（p+）后验估计方差
g	= zeros(2, 1);                   % 稳定增益矩阵
        % 稳定增益矩阵

for i=1:500
    for t=1:1001
            if t == 1
                pp = ( PI0^(-1) + bar_C' * R^(-1) * bar_C )^(-1);
                xp = pp * bar_C' * R^(-1) * Y1(:,t,i);
                
            else
                xn = Real_A{i}{t-1} * xp+u{i}{t}{1}*Real_B{i}{t-1}*[1;0.1];
                pn = Real_A{i}{t-1} * pp * Real_A{i}{t-1}' + bar_Q;
                g  = (pn * bar_C') * inv( bar_C * pn * bar_C' + R);
                xp = xn + g * (Y1(:,t,i) - bar_C * xn);
                pp = (eye(2) - g * bar_C) * pn;   
                
                
            end
            X_hat(:,t,i)	= xp';
         
    end
end
for i=1:500
    for t=1:1001
            if t == 1
                pp = ( PI0^(-1) + bar_C' * R^(-1) * bar_C )^(-1);
                xp = pp * bar_C' * R^(-1) * Y2(:,t,i);
                
            else
                xn =Real_A2{i}{t-1} * xp+u{i}{t}{2}*Real_B{i}{t-1}*[1;0.1];
                pn = Real_A2{i}{t-1} * pp * Real_A2{i}{t-1}' + bar_Q;
                g  = (pn * bar_C') * inv( bar_C * pn * bar_C' + R);
                xp = xn + g * (Y2(:,t,i) - bar_C * xn);
                pp = (eye(2) - g * bar_C) * pn;   
                
                
            end
            X_hat2(:,t,i)	= xp';
         
    end
end
Method6_Deviation_SUM       = zeros(1,1000);
Method6_Deviation_SUM2       = zeros(1,1000);
Method6_Standard_Deviation1	= zeros(1,1000);
Method6_Standard_Deviation2	= zeros(1,1000);
Method6_Standard_Deviation	= zeros(1,1000);
Method6_STD_dB              = zeros(1,1000);
for t=1:1000
     for i=1:500
       Method6_Deviation_SUM(t) = Method6_Deviation_SUM(t) + ...
                                     ( X_hat(1,t+1,i) - X1(1,t+1,i) )^2+( X_hat(2,t+1,i) - X1(2,t+1,i) )^2;
     
       Method6_Deviation_SUM2(t) = Method6_Deviation_SUM2(t) + ...
                                     ( X_hat2(1,t+1,i) - X2(1,t+1,i) )^2+( X_hat2(2,t+1,i) - X2(2,t+1,i) )^2;
     end
     Method6_Standard_Deviation1(t) = Method6_Deviation_SUM(t) / 500;
     Method6_Standard_Deviation2(t) = Method6_Deviation_SUM2(t) / 500;
     Method6_Standard_Deviation(t) = (Method6_Standard_Deviation1(t)^2+Method6_Standard_Deviation2(t)^2)^(1/2);
     Method6_STD_dB(t)             = 10 * log10( Method6_Standard_Deviation(t) ); 
end
n6=Method6_STD_dB(500);
epsilon11=cell(500,1);
for i=1:500
    epsilon11{i}=cell(1000,1);
    for t=1:1000
        
%   while 1
%          epsilon11{i}{t}=normrnd(0,1);   %%%正态分布    %%修改随机产生模型误差的方差
%     if  abs(epsilon11{i}{t})<=1;        %%%%把幅值改大，效果就好些了.5左右差不多，再大就更好，再小就差
%         break;
%     end
%   end
        epsilon11{i}{t}=unifrnd(0,0);   %%%%均匀分布
    end
    
end

% epsilon12=cell(500,1);
% for i=1:500
% %     while 1
% %          epsilon12{i}=normrnd(0,0);   %%正态分布   %%%修改随机产生模型误差的方差
% %     if  abs(epsilon12{i})<=1;        %%%%不同实验随机产生一个值，此值为0
% %         break;
% %     end
% %     end
%         epsilon12{i}=unifrnd(-1,1);   %%%%均匀分布      
% end

%%%其中，二者都等于0 ，表示无模型误差；epsilon12不等于0，epsilon11=0 ，表示不同实验模型误差随机变化，每次实验中模型误差随时间不再随机变化；
%%%%  epsilon11不等于0，epsilon12=0，表示每次实验模型误差随时间随机变化； 二者皆不等于0，表示不同实验随机产生一个值，每次实验时以此随机产生的值为均值随机变化。


epsilon1=cell(500,1);
for i=1:500
    epsilon1{i}=cell(1000,1);
    for t=1:1000
        epsilon1{i}{t}=epsilon12{i}+epsilon11{i}{t}; %%%%不同实验随机产生一个值，每次实验在此值基础上随时间随机变化
    end
end


% %矩阵A有模型误差
% A=[0.9802 0.0196;0 0.9802]+[0.0198;0.0]*epsilon*[0.0 5.0];

%矩阵A1没有模型误差，名义模型系数
z=cell(500,1);%先前变量Z(k)
for i=1:500
    z{i}=cell(1000,1);
    for t=1:1000
        z{i}{t}=0;
    end
end
z0=cell(500,1);%初始先前变量Z(k)
for i=1:500
    z0{i}=0;
end
u=cell(500,1);%隶属函数
for i=1:500
    u{i}=cell(1000,1);
    for t=1:1000
        u{i}{t}=cell(2,1);
        for k=1:2
            u{i}{t}{k}=0;
        end
    end
end
u0=cell(500,1);%初始隶属函数
for i=1:500
    u0{i}=cell(2,1);
    u0{i}{1}=1-1/(1+exp(pi/2));
    u0{i}{2}=1-u0{i}{1};
end
A1=cell(2,1);%Ai(0)
A1{1}=[0.9804 0.0196;0 0.9804];%%改变0.0196的值，相对的改变了模型误差的影响
A1{2}=[0.9902 0.0196;0 0.9902];%%改变0.0196的值，相对的改变了模型误差的影响
A11=cell(500,1);%Ak(0)
for i=1:500
    A11{i}=cell(1000,1);
    for t=1:1000
        A11{i}{t}=[0 0;0 0];
    end
end

%矩阵AA有模型误差
AA=cell(500,1);%Ai(epsilon)
AAA=cell(500,1);%Ak(epsilon)
for i=1:500
    AA{i}=cell(1000,1);
    AAA{i}=cell(1000,1);
    for t=1:1000
        AA{i}{t}=cell(2,1);
        AAA{i}{t}=[0 0;0 0];
        for k=1:2
            AA{i}{t}{k}=A1{k}+[0 0.099;0 0]*(epsilon1{i}{t});%%%%改变0.099的值，即改变模型误差的大小，后面相应的灵敏度惩罚算法等都要修改
        end
    end

end

U=1;      
B=[0.0196;0];
C=[1 0;0 1];
D=[1 -1];
R=1.0000;
Q=[1.9608 0.0195;0.0195 1.9605];
PI0=[1 0;0 1]; 


%%构造标准Kalman滤波需要的各个向量（optimal）,向量的每个元素置为0
x0=cell(500,1);%Xi(k)初始真实状态
for i=1:500
    x0{i}=cell(2,1);
    for k=1:2
        x0{i}{k}=[0;0];
    end
end

xx0=cell(500,1);%X(k)初始真实状态
for i=1:500
    xx0{i}=[0;0];
end

x=cell(500,1);%Xi(k)真实状态
for i=1:500
    x{i}=cell(1000,1);  
    for t=1:1000
        x{i}{t}=cell(2,1);
        for k=1:2
            x{i}{t}{k}=[0;0];
        end
    end
end

xx=cell(500,1);%X(k)真实状态
for i=1:500
    xx{i}=cell(1000,1);
    for t=1:1000
        xx{i}{t}=[0;0];
    end
end

%%构造需要的各个向量
y0=cell(500,1);%Yi(k)初始输出
for i=1:500
    y0{i}=cell(2,1);
    for k=1:2
        y0{i}{k}=0;
    end
end

yy0=cell(500,1);%Y(k)初始输出
for i=1:500
    yy0{i}=0;
end

y=cell(500,1);%Yi(k)输出
for i=1:500
    y{i}=cell(1000,1);    
    for t=1:1000
        y{i}{t}=cell(2,1);
        for k=1:2
            y{i}{t}{k}=0;
        end
    end
end

yy=cell(500,1);%Y(k)输出
for i=1:500
    yy{i}=cell(1000,1);
    for t=1:1000
        yy{i}{t}=0;
    end
end

e0=cell(500,1);%初始输出误差:Yi(k)-D*hatXi(k|k-1)
for i=1:500
    e0{i}=cell(2,1);
    for k=1:2
        e0{i}{k}=0;
    end
end

e=cell(500,1);%输出误差:Yi(k)-D*hatXi(k|k-1)
for i=1:500
    e{i}=cell(1000,1);  
    for t=1:1000
        e{i}{t}=cell(2,1);
        for k=1:2
            e{i}{t}{k}=0;
        end
    end
end



% m1=cell(500,1);%输出显示用
% for i=1:500
% m1{i}=cell(1000,1);   
% for t=1:1000
%    m1{i}{t}=0;
% end
% end
% 
% m2=cell(500,1);%输出显示用
% for i=1:500
% m2{i}=cell(1000,1);   
% for t=1:1000
% m2{i}{t}=0;
% end
% end

Ree0=cell(500,1);%初始输出误差:D*Pi(k|k-1)*D'+ui^2*R
for i=1:500
    Ree0{i}=cell(2,1);
    for k=1:2
        Ree0{i}{k}=0;
    end
end

Ree=cell(500,1);%输出误差:D*Pi(k|k-1)*D'+ui^2*R
for i=1:500
    Ree{i}=cell(1000,1);  
    for t=1:1000
        Ree{i}{t}=cell(2,1);
        for k=1:2
            Ree{i}{t}{k}=0;
        end
    end
end

P10=cell(500,1);%初始误差方差阵Pi(k|k-1)
for i=1:500
    P10{i}=cell(2,1);
    for k=1:2
        P10{i}{k}=[0 0;0 0];
    end
end

P1=cell(500,1);%误差方差阵Pi(k|k-1)
for i=1:500
    P1{i}=cell(1000,1); 
    for t=1:1000
        P1{i}{t}=cell(2,1);
        for k=1:2
            P1{i}{t}{k}=[0 0;0 0];
        end
    end
end

P20=cell(500,1);%初始误差方差阵Pi(k|k)
for i=1:500
    P20{i}=cell(2,1);
    for k=1:2
        P20{i}{k}=[0 0;0 0];
    end
end

P2=cell(500,1);%误差方差阵Pi(k|k)
for i=1:500
    P2{i}=cell(1000,1);   
    for t=1:1000
        P2{i}{t}=cell(2,1);
        for k=1:2
            P2{i}{t}{k}=[0 0;0 0];
        end
    end
end


P30=cell(500,1);%初始误差方差阵hat Pi(k|k)
for i=1:500
    P30{i}=cell(2,1);
    for k=1:2
        P30{i}{k}=[0 0;0 0];
    end
end

P3=cell(500,1);%误差方差阵 hat Pi(k|k)
for i=1:500
    P3{i}=cell(1000,1);   
    for t=1:1000
        P3{i}{t}=cell(2,1);
        for k=1:2
            P3{i}{t}{k}=[0 0;0 0];
        end
    end
end

Q20=cell(500,1);%初始hat Qi 
for i=1:500
    Q20{i}=cell(2,1);
    for k=1:2
        Q20{i}{k}=[0 0;0 0];
    end
end

Q2=cell(500,1);% hat Qi
for i=1:500
    Q2{i}=cell(1000,1);   
    for t=1:1000
        Q2{i}{t}=cell(2,1);
        for k=1:2
            Q2{i}{t}{k}=[0 0;0 0];
        end
    end
end

A20=cell(500,1);%初始hat Ai
for i=1:500
    A20{i}=cell(2,1);
    for k=1:2
        A20{i}{k}=[0 0;0 0];
    end
end

A2=cell(500,1);%A hat
for i=1:500
    A2{i}=cell(1000,1);   
    for t=1:1000
        A2{i}{t}=cell(2,1);
        for k=1:2
            A2{i}{t}{k}=[0 0;0 0];
        end
    end
end

B20=cell(500,1);%初始hat Bi 
for i=1:500
    B20{i}=cell(2,1);
    for k=1:2
        B20{i}{k}=[0;0];
    end
end

B2=cell(500,1);%hat Bi
for i=1:500
    B2{i}=cell(1000,1);   
    for t=1:1000
        B2{i}{t}=cell(2,1);
        for k=1:2
            B2{i}{t}{k}=[0;0];
        end
    end
end

C20=cell(500,1);%初始hat Ci 
for i=1:500
    C20{i}=cell(2,1);
    for k=1:2
        C20{i}{k}=[0 0;0 0];
    end
end

C2=cell(500,1);%hat Ci
for i=1:500
    C2{i}=cell(1000,1);   
    for t=1:1000
        C2{i}{t}=cell(2,1);
        for k=1:2
            C2{i}{t}{k}=[0 0;0 0];
        end
    end
end

R20=cell(500,1);%初始hat Ri 
for i=1:500
    R20{i}=cell(2,1);
    for k=1:2
        R20{i}{k}=0;
    end
end

R2=cell(500,1);%hat Ri
for i=1:500
    R2{i}=cell(1000,1);   
    for t=1:1000
        R2{i}{t}=cell(2,1);
        for k=1:2
            R2{i}{t}{k}=0;
        end
    end
end



cx0=cell(500,1);%初始估计状态hat Xi(k|k)
for i=1:500
    cx0{i}=cell(2,1);
    for k=1:2
        cx0{i}{k}=[0;0];
    end
end

cx=cell(500,1);%估计状态hat Xi(k|k)
for i=1:500
    cx{i}=cell(1000,1); 
    for t=1:1000
        cx{i}{t}=cell(2,1);     
        for k=1:2
            cx{i}{t}{k}=[0;0];
        end
    end
end

cxx0=cell(500,1);%初始估计状态hat X(k|k)
for i=1:500
    cxx0{i}=cell(2,1);
    for k=1:2
        cxx0{i}{k}=[0;0];
    end
end

cxx=cell(500,1);%估计状态hat X(k|k)
for i=1:500
    cxx{i}=cell(1000,1); 
    for t=1:1000
        cxx{i}{t}=[0;0];     
    end
end

cx10=cell(500,1);%初始由状态方程计算的估计状态hat Xi(k|k-1)
for i=1:500
    cx10{i}=cell(2,1);
    for k=1:2
        cx10{i}{k}=[0;0];
    end
end

cx1=cell(500,1);%由状态方程计算的估计状态hat Xi(k|k-1)
for i=1:500
    cx1{i}=cell(1000,1); 
    for t=1:1000
        cx1{i}{t}=cell(2,1);
        for k=1:2
            cx1{i}{t}{k}=[0;0];
        end
    end
end

W=cell(500,1);%W(k)
for i=1:500
    W{i}=cell(1000,1);
    for t=1:1000
       W{i}{t}=[normrnd(0,1);normrnd(0,1)]; 
    end
end

V=cell(500,1);%V(k)
for i=1:500
    V{i}=cell(1000,1);
    for t=1:1000
       V{i}{t}=normrnd(0,1); 
    end
end

%产生500次实验数据

for i=1:500 %500次实验
    for t=2:1000
        for k=1:2
            x0{i}{k}=[0;0];
            y0{i}{k}=D*x0{i}{k}+u0{i}{k}*R^(1/2)*V{i}{1};%定义初始输出
            x{i}{1}{k}=AA{i}{1}{k}*x0{i}{k}+u0{i}{k}*B*U+u0{i}{k}*C*Q^(1/2)*W{i}{1}; %用A产生的数据，后面用A估计 %%%注意，这里的t=1时刻用到的AA的系数也是第一个元素，若按照状态方程，应该为0元素
            y{i}{1}{k}=D*x{i}{1}{k}+u0{i}{k}*R^(1/2)*V{i}{1};
            xx{i}{1}=x{i}{1}{1}+x{i}{1}{2};
            %yy{i}{1}=y{i}{1}{1}+y{i}{1}{2};
            z{i}{1}=(xx{i}{1}(1)*3.1416/180-0.0098*xx{i}{1}(2)*3.1416/180);
            u{i}{1}{1}=1-1/(1+exp(3.1416/2-z{i}{1}));
            u{i}{1}{2}=1-u{i}{1}{1};
            x{i}{t}{k}=AA{i}{t}{k}*x{i}{t-1}{k}+u{i}{t-1}{k}*B*U+u{i}{t-1}{k}*C*Q^(1/2)*W{i}{t};%理想状态及输出计算 %%%注意，这里的t=1时刻用到的AA的系数也是第一个元素，若按照状态方程，应该为0元素
            y{i}{t}{k}=D*x{i}{t}{k}+u{i}{t-1}{k}*R^(1/2)*V{i}{t};
            %每个时刻的x为2个元素的列向量，每个时刻的y为一个元素的标量            
            xx{i}{t}=x{i}{t}{1}+x{i}{t}{2};
            z{i}{t}=(xx{i}{t}(1)*3.1416/180-0.0098*xx{i}{t}(2)*3.1416/180);
            u{i}{t}{1}=1-1/(1+exp(3.1416/2-z{i}{t}));
            u{i}{t}{2}=1-u{i}{t}{1};
            %yy{i}{t}=yy{i}{t}+y{i}{t}{k};
        end
    end
end

%%%%%%求期望
%%%%%%求期望
%%%%%%求期望
%%模型误差的均值必须为0，下面这些式子才对。
%%表达式中的期望
H11=cell(500,1);
for i=1:500
    H11{i}=cell(1000,1);
    for t=1:1000
        H11{i}{t}=cell(2,1);
        for k=1:2
            H11{i}{t}{k}=[A1{k}';u{i}{t}{k}*C']*D'*R^(-1)*D*[A1{k} u{i}{t}{k}*C]+...
                (1/3+0/3)*[([0 0.099;0 0])';zeros(2)]*D'*R^(-1)*D*[([0 0.099;0 0])';zeros(2)]';         %%%% 修改此参数
        %%%%%%%%%%%%%%%%%%此处有个模型误差的方差均值
        end
    end
end
H110=cell(500,1);
for i=1:500
    H110{i}=cell(2,1);
    for k=1:2
        H110{i}{k}=[A1{k}';u0{i}{k}*C']*D'*R^(-1)*D*[A1{k} u0{i}{k}*C]+...
                (1/3+0/3)*[([0 0.099;0 0])';zeros(2)]*D'*R^(-1)*D*[([0 0.099;0 0])';zeros(2)]';
    end
end
H12=cell(500,1);
for i=1:500
    H12{i}=cell(1000,1);
    for t=1:1000
        H12{i}{t}=cell(2,1);
        for k=1:2
            H12{i}{t}{k}=[A1{k}';u{i}{t}{k}*C']*D'; 
        end
    end
end
H120=cell(500,1);
for i=1:500
    H120{i}=cell(2,1);
    for k=1:2
        H120{i}{k}=[A1{k}';u0{i}{k}*C']*D';
    end
end
H13=cell(500,1);
for i=1:500
    H13{i}=cell(1000,1);
    for t=1:1000
        H13{i}{t}=cell(2,1);
        for k=1:2
            H13{i}{t}{k}=[A1{k}';u{i}{t}{k}*C']*D'*R^(-1)*D*u{i}{t}{k}*B;
        end
    end
end
H130=cell(500,1);
for i=1:500
    H130{i}=cell(2,1);
    for k=1:2
        H130{i}{k}=[A1{k}';u0{i}{k}*C']*D'*R^(-1)*D*u0{i}{k}*B;
    end
end
%%表达式中的期望
% H14=[A1';B']*C'*R^(-1)*C*A1+(1+1)*[([0.0198;0.0]*[0.0 5.0])';0 0;0 0]*C'*R^(-1)*C*[0.0198;0.0]*[0.0 5.0];        %%%%%%%%%%%%%%%%%%此处有个模型误差的方差均值 %%表达式中的期望
H14=cell(500,1);
for i=1:500
    H14{i}=cell(1000,1);
    for t=1:1000
        H14{i}{t}=cell(2,1);
        for k=1:2
            H14{i}{t}{k}=H11{i}{t}{k}*[eye(2);zeros(2)];
        end
    end
end
H140=cell(500,1);
for i=1:500
    H140{i}=cell(2,1);
    for k=1:2
        H140{i}{k}=H110{i}{k}*[eye(2);zeros(2)];
    end
end
G11=cell(500,1);
for i=1:500
    G11{i}=cell(1000,1);
    for t=1:1000
        G11{i}{t}=cell(2,1);
        for k=1:2
            G11{i}{t}{k}=H11{i}{t}{k}-[A1{k}';u{i}{t}{k}*C']*D'*R^(-1)*D*[A1{k} u{i}{t}{k}*C];
        end
    end
end
G110=cell(500,1);
for i=1:500
    G110{i}=cell(2,1);
    for k=1:2
        G110{i}{k}=H110{i}{k}-[A1{k}';u0{i}{k}*C']*D'*R^(-1)*D*[A1{k} u0{i}{k}*C];
    end
end
T11=cell(500,1);
for i=1:500
    T11{i}=cell(1000,1);
    for t=1:1000
        T11{i}{t}=cell(2,1);
        for k=1:2
            T11{i}{t}{k}=H13{i}{t}{k}-[A1{k}';u{i}{t}{k}*C']*D'*R^(-1)*D*u{i}{t}{k}*B;
        end
    end
end
T110=cell(500,1);
for i=1:500
    T110{i}=cell(2,1);
    for k=1:2
        T110{i}{k}=H130{i}{k}-[A1{k}';u0{i}{k}*C']*D'*R^(-1)*D*u0{i}{k}*B;
    end
end
for i=1:500 %对500次实验数据分别运用新算法
    for t=2:1000
        for k=1:2
            %t=1时刻的估计状态
            P20{i}{k}=(PI0^(-1)+D'*R^(-1)*D)^(-1); %P0|0 ，由于C中无模型误差才可以这样写，否则需要前面给出      
            cx0{i}{k}=P20{i}{k}*D'*R^(-1)*y0{i}{k};      %hat x0|0
            
            P30{i}{k}=(P20{i}{k}^(-1)+[G110{i}{k}(1,1) G110{i}{k}(1,2);G110{i}{k}(2,1) G110{i}{k}(2,2)])^(-1); %hat P0|0  %%%%%%%%%%%%%%%%%%此处有个模型误差的方差均值
            Q20{i}{k}=(Q^(-1)+[G110{i}{k}(3,3) G110{i}{k}(3,4);G110{i}{k}(4,3) G110{i}{k}(4,4)]-...
                [G110{i}{k}(1,3) G110{i}{k}(1,4);G110{i}{k}(2,3) G110{i}{k}(2,4)]'*P30{i}{k}*[G110{i}{k}(1,3) G110{i}{k}(1,4);G110{i}{k}(2,3) G110{i}{k}(2,4)])^(-1);%hat Q
            C20{i}{k}=u0{i}{k}*C-A1{k}*P30{i}{k}*[G110{i}{k}(1,3) G110{i}{k}(1,4);G110{i}{k}(2,3) G110{i}{k}(2,4);];  %hat C
            B20{i}{k}=u0{i}{k}*B-A1{k}*P30{i}{k}*[T110{i}{k}(1,1);T110{i}{k}(2,1)]-C20{i}{k}*Q20{i}{k}*([T110{i}{k}(3,1);T110{i}{k}(4,1)]-[G110{i}{k}(3,1) G110{i}{k}(3,2);G110{i}{k}(4,1) G110{i}{k}(4,2)]*P30{i}{k}*[T110{i}{k}(1,1);T110{i}{k}(2,1)]);            %hat B
            %     A20{i}=A1-A1*P30{i}*[1 0 0 0;0 1 0 0]*H13-B20{i}*Q20{i}*([-P30{i}*[G11(1,3) G11(1,4);G11(2,3) G11(2,4)];1 0;0 1])'*H13+...
            %         (A1*P30{i}*A1'+B20{i}*Q20{i}*B20{i}')*C'*R^(-1)*C*A1;%hat A
            A20{i}{k}=(A1{k}-C20{i}{k}*Q20{i}{k}*[G110{i}{k}(3,1) G110{i}{k}(3,2);G110{i}{k}(4,1) G110{i}{k}(4,2)])*(eye(2)-P30{i}{k}*[G110{i}{k}(1,1) G110{i}{k}(1,2);G110{i}{k}(2,1) G110{i}{k}(2,2)]);%%hat A
            P1{i}{1}{k}=A1{k}*P30{i}{k}*A1{k}'+C20{i}{k}*Q20{i}{k}*C20{i}{k}';%P1为Pi+1|i
            %     Ree{i}{1}{k}=u(k)^2*R+D*P1{i}{1}{k}*D';
            P2{i}{1}{k}=(P1{i}{1}{k}^(-1)+D'*R^(-1)*D)^(-1);%P2为Pi+1|i+1
            cx{i}{1}{k}=A20{i}{k}*cx0{i}{k}+B20{i}{k}*U+P2{i}{1}{k}*(P1{i}{1}{k}^(-1)*(A1{k}*P30{i}{k}*[eye(2) zeros(2)]+C20{i}{k}*Q20{i}{k}*[zeros(2) eye(2)])*...
                H120{i}{k}*R^(-1)*y{i}{1}{k}-D'*R^(-1)*D*(A20{i}{k}*cx0{i}{k}+B20{i}{k}*U));
            
            P3{i}{1}{k}=(P2{i}{1}{k}^(-1)+[G110{i}{k}(1,1) G110{i}{k}(1,2);G110{i}{k}(2,1) G110{i}{k}(2,2)])^(-1);%hat P1|1 %%%%%%%%%%%%%%%%%%此处有个模型误差的方差均值
            Q2{i}{1}{k}=(Q^(-1)+[G110{i}{k}(3,3) G110{i}{k}(3,4);G110{i}{k}(4,3) G110{i}{k}(4,4)]-...
                [G110{i}{k}(1,3) G110{i}{k}(1,4);G110{i}{k}(2,3) G110{i}{k}(2,4)]'*P3{i}{1}{k}*[G110{i}{k}(1,3) G110{i}{k}(1,4);G110{i}{k}(2,3) G110{i}{k}(2,4)])^(-1);%hat Q
            C2{i}{1}{k}=u0{i}{k}*C-A1{k}*P3{i}{1}{k}*[G110{i}{k}(1,3) G110{i}{k}(1,4);G110{i}{k}(2,3) G110{i}{k}(2,4)];  %hat C
            B2{i}{1}{k}=u0{i}{k}*B-A1{k}*P3{i}{1}{k}*[T110{i}{k}(1,1);T110{i}{k}(2,1)]-C2{i}{1}{k}*Q2{i}{1}{k}*([T110{i}{k}(3,1);T110{i}{k}(4,1)]-[G110{i}{k}(3,1) G110{i}{k}(3,2);G110{i}{k}(4,1) G110{i}{k}(4,2)]*P3{i}{1}{k}*[T110{i}{k}(1,1);T110{i}{k}(2,1)]); %hat B
            %     A2{i}{1}=A1-A1*P3{i}{1}*[1 0 0 0;0 1 0 0]*H13-B2{i}{1}*Q2{i}{1}*([-P3{i}{1}*[G11(1,3) G11(1,4);G11(2,3) G11(2,4)];1 0;0 1])'*H13+...
            %         (A1*P3{i}{1}*A1'+B2{i}{1}*Q2{i}{1}*B2{i}{1}')*C'*R^(-1)*C*A1;%hat A
            A2{i}{1}{k}=(A1{k}-C2{i}{1}{k}*Q2{i}{1}{k}*[G110{i}{k}(3,1) G110{i}{k}(3,2);G110{i}{k}(4,1) G110{i}{k}(4,2)])*(eye(2)-P3{i}{1}{k}*[G110{i}{k}(1,1) G110{i}{k}(1,2);G110{i}{k}(2,1) G110{i}{k}(2,2)]);%%hat A

            %误差惩罚方法
        
            P1{i}{t}{k}=A1{k}*P3{i}{t-1}{k}*A1{k}'+C2{i}{t-1}{k}*Q2{i}{t-1}{k}*C2{i}{t-1}{k}';%P1为Pi+1|i
            %     Ree{i}{t}=R+C*P1{i}{t}*C';
            P2{i}{t}{k}=(P1{i}{t}{k}^(-1)+D'*R^(-1)*D)^(-1);%P2为Pi+1|i+1
            cx{i}{t}{k}=A2{i}{t-1}{k}*cx{i}{t-1}{k}+B2{i}{t-1}{k}*U+P2{i}{t}{k}*(P1{i}{t}{k}^(-1)*(A1{k}*P3{i}{t-1}{k}*[eye(2) zeros(2)]+C2{i}{t-1}{k}*Q2{i}{t-1}{k}*[zeros(2) eye(2)])*...
                H12{i}{t-1}{k}*R^(-1)*y{i}{t}{k}-D'*R^(-1)*D*(A2{i}{t-1}{k}*cx{i}{t-1}{k}+B2{i}{t-1}{k}*U));
            P3{i}{t}{k}=(P2{i}{t}{k}^(-1)+[G11{i}{t-1}{k}(1,1) G11{i}{t-1}{k}(1,2);G11{i}{t-1}{k}(2,1) G11{i}{t-1}{k}(2,2)])^(-1);%hat Pi|i %%%%%%%%%%%%%%%%%%此处有个模型误差的方差均值
            Q2{i}{t}{k}=(Q^(-1)+[G11{i}{t-1}{k}(3,3) G11{i}{t-1}{k}(3,4);G11{i}{t-1}{k}(4,3) G11{i}{t-1}{k}(4,4)]-...
                [G11{i}{t-1}{k}(1,3) G11{i}{t-1}{k}(1,4);G11{i}{t-1}{k}(2,3) G11{i}{t-1}{k}(2,4)]'*P3{i}{t}{k}*[G11{i}{t-1}{k}(1,3) G11{i}{t-1}{k}(1,4);G11{i}{t-1}{k}(2,3) G11{i}{t-1}{k}(2,4)])^(-1);%hat Q
            C2{i}{t}{k}=u{i}{t-1}{k}*C-A1{k}*P3{i}{t}{k}*[G11{i}{t-1}{k}(1,3) G11{i}{t-1}{k}(1,4);G11{i}{t-1}{k}(2,3) G11{i}{t-1}{k}(2,4)];  %hat C
            B2{i}{t}{k}=u{i}{t-1}{k}*B-A1{k}*P3{i}{t}{k}*[T11{i}{t-1}{k}(1,1);T11{i}{t-1}{k}(2,1)]-C2{i}{t}{k}*Q2{i}{t}{k}*([T11{i}{t-1}{k}(3,1);T11{i}{t-1}{k}(4,1)]-[G11{i}{t-1}{k}(3,1) G11{i}{t-1}{k}(3,2);G11{i}{t-1}{k}(4,1) G11{i}{t-1}{k}(4,2)]*P3{i}{t}{k}*[T11{i}{t-1}{k}(1,1);T11{i}{t-1}{k}(2,1)]); %hat B
            %     A2{i}{t}=A1-A1*P3{i}{t}*[1 0 0 0;0 1 0 0]*H13-B2{i}{t}*Q2{i}{t}*([-P3{i}{t}*[G11(1,3) G11(1,4);G11(2,3) G11(2,4)];1 0;0 1])'*H13+...
            %         (A1*P3{i}{t}*A1'+B2{i}{t}*Q2{i}{t}*B2{i}{t}')*C'*R^(-1)*C*A1;%hat A
            A2{i}{t}{k}=(A1{k}-C2{i}{t}{k}*Q2{i}{t}{k}*[G11{i}{t-1}{k}(3,1) G11{i}{t-1}{k}(3,2);G11{i}{t-1}{k}(4,1) G11{i}{t-1}{k}(4,2)])*(eye(2)-P3{i}{t}{k}*[G11{i}{t-1}{k}(1,1) G11{i}{t-1}{k}(1,2);G11{i}{t-1}{k}(2,1) G11{i}{t-1}{k}(2,2)]);%%hat A
            %cxx{i}{t}=cxx{i}{t}+cx{i}{t}{k};
        end
    end
    %cxx{i}{1}=cx{i}{1}{1}+cx{i}{1}{2}+cx{i}{1}{3};
end


%通过计算总体均值或者叫统计平均值来估计实际误差方差
a1=cell(1000,1);%求500次和
for t=1:1000
    a1{t}=0;    
end

a2=cell(1000,1);%求500次和
for t=1:1000
    a2{t}=0;    
end

b1=cell(1000,1);%求500次平均值
for t=1:1000
    b1{t}=0;    
end

b2=cell(1000,1);%求500次平均值
for t=1:1000
    b2{t}=0;    
end

b=cell(1000,1);%求500次平均值
for t=1:1000
    b{t}=0;    
end

c1=cell(1000,1);%求500次平均值的db
for t=1:1000
    c1{t}=0;    
end

for t=1:1000
     for i=1:500
         a1{t}=((cx{i}{t}{1}(1))-x{i}{t}{1}(1))^2+((cx{i}{t}{1}(2))-x{i}{t}{1}(2))^2+a1{t};
         a2{t}=((cx{i}{t}{2}(1))-x{i}{t}{2}(1))^2+((cx{i}{t}{2}(2))-x{i}{t}{2}(2))^2+a2{t};
     end
     b1{t}=a1{t}/500;
     b2{t}=a2{t}/500;
     b{t}=(b1{t}^2+b2{t}^2)^(1/2);
     c1{t}=10*log10(b{t});
 
end

n7=cell2mat(c1);
figure(2);
semilogx(Method1_STD_dB','-m'); hold on;%0.035
plot(200,Method1_STD_dB(200),'b+');hold on;

semilogx(Method2_STD_dB','-b'); hold on;%0.2859
plot(200,Method2_STD_dB(200),'bD');hold on;


semilogx(Method3_STD_dB','-g'); hold on;%0.75
plot(200,Method3_STD_dB(200),'bV');hold on;


semilogx(Method4_STD_dB','-k'); hold on;%0.949
plot(200,Method4_STD_dB(200),'bo');hold on;

semilogx(Method5_STD_dB','-c'); hold on;
plot(200,Method5_STD_dB(200),'b*');hold on;

semilogx(Method6_STD_dB','-r'); hold on;
plot(200,Method6_STD_dB(200),'bs');hold on;
% plot(n,'-b')
semilogx(n7,'-y')
semilogx(log10(10^300),n7(300),'-kP')
xlabel('Sampled instant i');
ylabel('Estimation errors variance (dB)');
% legend('鲁棒状态估计','基于名义参数的Kalman','基于实际模型的Kalman');
legend('Robust state estimation gamma=0.035','Robust state estimation gamma=0.2859','Robust state estimation gamma=0.75','Robust state estimation gamma=0.9596','Kalman filter based on nominal parameters','Kalman filter based on actual model');
