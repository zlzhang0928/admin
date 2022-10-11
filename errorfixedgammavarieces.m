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
    esk =-0.8508;
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
                      

S=[0,0;0,0.099];
T1=[0,0;0,0.099];      
T2=[0,0;0,0];
    c99=cell(1,100);
for i=1:100
    c99{i}=0; 
end
for  gamma =1:100
     gamma1=(1-gamma/100)/(gamma/100);
     
    Pi1i= zeros(2,2);
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
c99{gamma}=Method1_STD_dB(500);
fprintf('预测模型1+鲁棒状态估计方法已完成\n');
end
c99=cell2mat(c99);
gaa=[0.01:0.01:1];
hh=[0.01:0.01:1];
c99=c99*(gaa/hh);
fprintf('鲁棒状态估计for循环已完成\n');
 
 
 
 
 
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

gaa=[0.01:0.01:1];
hh=[0.01:0.01:1];
n5=n5*(gaa/hh);
b=zeros(1,100);
for i=1:100
    b(i)=n5;
end
fprintf('名义参数Kalman方法已完成\n');
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
                xn = Real_A_Temp * xp+u{i}{t}{1}*Real_B_Temp*[1;0.1];
                pn = Real_A_Temp * pp * Real_A_Temp' + bar_Q;
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
                xn = Real_A2_Temp * xp+u{i}{t}{2}*Real_BB_Temp*[1;0.1];
                pn = Real_A2_Temp * pp * Real_A2_Temp' + bar_Q;
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
gaa=[0.01:0.01:1];
hh=[0.01:0.01:1];
n6=n6*(gaa/hh);
a=zeros(1,100);
for i=1:100
    a(i)=n6;
end
fprintf('实际参数Kalman方法已完成\n');
figure(1);
plot(gaa,c99, ':m',gaa, b, '-b',gaa, a, '-r');
xlabel('Design parameters \gamma');
ylabel('Estimation errors variance (dB)');hold on;
plot(0.5,b(5),'-m*');hold on;%%名义
plot(0.5,a(5),'-bs');hold on;%%实际
plot(0.5,c99(5),'-rD');hold on;%%鲁棒
