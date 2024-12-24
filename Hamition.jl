# 定义哈密顿量
#=
LinearAlgebra 是一个标准库，提供了许多线性代数操作的功能。
可以使用 LinearAlgebra 模块来执行矩阵和向量的基本运算、求解线性方程组、
计算特征值和特征向量等操作。
=#
using LinearAlgebra
#=
:: 是用于类型声明和类型注解的符号。
:: 是可选的，因为 Julia 是一种动态类型语言，可以自动推断类型。
:: 可以帮助编译器进行更好的优化，并提高代码的可读性。
=#
function Kane_Mele_Hamiltonian(kx::Float64,ky::Float64)
    
    # 键盘输入 \lambda ，然后按一下Tab键可以打出 λ
    t=1;λR=0.05;λSO=0.06;λv=0.1;
    
    # 泡利矩阵
    s0=[1 0;0 1];
    sx=[0 1;1 0];
    sy=[0 -im;im 0];
    sz=[1 0;0 -1];
    
    # 第一项：最近邻
    # 在Julia中 [kx,ky] 给出列向量，[kx ky] 给出行向量
    k=[kx,ky]; 
    # 单斜杠 / 表示通常的浮点数除法
    # 双斜杠 // 用于创建有理数（Rational）对象，表示精确的分数
    r1=[1,0];r2=[1/2,sqrt(3)/2];r3=[-1/2,sqrt(3)/2];
    r4=[-1,0];r5=[-1/2,-sqrt(3)/2];r6=[1/2,-sqrt(3)/2];
    # .* 是元素级的逐个相乘运算符，* 是矩阵乘法运算符
    # dot 表示点乘
    h1=[0 t;t 0]+
        [0 0;t 0].*exp(im*dot(k,r2))+[0 0;t 0].*exp(im*dot(k,r3))+
        [0 t;0 0].*exp(im*dot(k,r5))+[0 t;0 0].*exp(im*dot(k,r6))
    H1=kron(h1,s0)
    
    # 第二项：次近邻 SOC
    # 定义匿名函数，可以通过 ->（箭头符号）创建
    # hcat 将两个矩阵水平拼接，vcat 将两个矩阵垂直拼接
    matrix1=(x1,x2)->vcat(hcat(x1,zeros(2,2)),hcat(zeros(2,2),x2))
    h2=im*λSO*sz
    H2=matrix1(-h2,h2)*exp(im*dot(k,r1))+
        matrix1(h2,-h2)*exp(im*dot(k,r2))+
        matrix1(-h2,h2)*exp(im*dot(k,r3))+
        matrix1(h2,-h2)*exp(im*dot(k,r4))+
        matrix1(-h2,h2)*exp(im*dot(k,r5))+
        matrix1(h2,-h2)*exp(im*dot(k,r6))
    
    # 第三项：Rashba SOC
    d1=[0,1,0];d2=[-sqrt(3)/2,-1/2,0];d3=[sqrt(3)/2,-1/2,0];
    d4=[0,-1,0];d5=[sqrt(3)/2,1/2,0];d6=[-sqrt(3)/2,1/2,0];
    s=[sx,sy,zeros(2,2)]
    # dot 表示叉乘，这里对结果只取z方向
    D0=zeros(2,2);
    D1=cross(s,-d1)[3];
    D2=cross(s,-d2)[3];
    D3=cross(s,-d3)[3];
    D4=cross(s,-d4)[3];
    D5=cross(s,-d5)[3];
    D6=cross(s,-d6)[3];
    # 定义函数的另一种形式
    matrix2(x1,x2)=im*λR*vcat(hcat(zeros(2,2),x1),hcat(x2,zeros(2,2)))
    H3=matrix2(D1,D4)+
        matrix2(D0,D5)*exp(im*dot(k,r2))+
        matrix2(D0,D6)*exp(im*dot(k,r3))+
        matrix2(D2,D0)*exp(im*dot(k,r5))+
        matrix2(D3,D0)*exp(im*dot(k,r6))
    
    # 第四项：交错的子晶格势
    h4=[λv 0;0 -λv]
    H4=kron(h4,s0)
    
    # 总哈密顿量
    H=H1+H2+H3+H4
    #=
    return 语句用于从函数中返回一个值。
    是可选的，因为 Julia 具有隐式的返回语义，即函数的最后一个表达式的值将作为返回值。
    但在某些情况下，需要显式地使用 return 来提前退出函数或返回一个特定的值
    =#
    return H
end

# 正格子
a1=[1,0,0];a2=[1/2,sqrt(3)/2,0];a3=[0,0,1];
aa=[a1,a2,a3];
# 倒格子
V=dot(a1,cross(a2,a3));
b1=2pi*cross(a2,a3)/V;b2=2pi*cross(a3,a1)/V;b3=2pi*cross(a1,a2)/V;
bb=[b1,b2,b3];
# 高对称点坐标
G=[0.,0.,0.];K=[2/3,1/3,0];M=[1/2,0,0];
# 将分数坐标变换为直角坐标
ktrans=k -> k[1]*b1+k[2]*b2+k[3]*b3
GG=ktrans(G);KK=ktrans(K);MM=ktrans(M);
GKM=[GG,KK,MM]
# 使用 range 函数来生成两个向量之间的线性插值
# 前段路径中的最后一个点和后段路径中的第一个点相同，用切片操作去除
GM=range(GG,MM,length=40);
MK=range(MM,KK,length=20)[2:end];
KG=range(KK,GG,length=40)[2:end];
# 用 vcat 把三段路径连接到一起
kpath=vcat(GM,MK,KG);
# 高对称点的位置索引
kk=[1,40,59,98];
# 计算kpath中相邻两个点之间的距离
dis=[norm(kpath[i+1]-kpath[i]) for i in 1:length(kpath)-1];
pushfirst!(dis,0);
# 计算前i个dis的距离之和
dis0=[sum(dis[1:i]) for i in 1:length(dis)];
# 高对称点的距离
KK=[dis0[1],dis0[40],dis0[59],dis0[98]]

# 求解哈密顿量
vals=zeros(4,length(kpath))
for i in 1:length(kpath)
    val,vec=eigen(Kane_Mele_Hamiltonian(kpath[i][1],kpath[i][2]))
    vals[:,i]=real.(val)
end

# 画图
using Plots
plot(dis0,vals[1,:],color=:red,label="band1",lw=2)
plot!(dis0,vals[2,:],color=:blue,label="band2",lw=2)
plot!(dis0,vals[3,:],color=:green,label="band3",lw=2)
plot!(dis0,vals[4,:],color=:purple,label="band4",lw=2)
xticks!(KK,["Γ","M","K","Γ"])
xlabel!("Wave vector")
ylabel!("Energy(eV)")
title!("Kane Mele Model")