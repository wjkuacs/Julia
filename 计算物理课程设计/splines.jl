using LinearAlgebra
using Plots

# 定义数据点
data = [(1, 1), (2, 3), (4, 4), (5, 2)]
x = [d[1] for d in data]
y = [d[2] for d in data]

# 计算每个子区间的差值
h = [x[i+1] - x[i] for i in 1:length(x)-1]

# 计算方程右边的向量
d = [6((y[i+2] - y[i+1])/h[i+1] - (y[i+1] - y[i])/h[i]) for i in 1:length(x)-2]

# 构造方程组的系数矩阵并求解
n = length(x) - 1
mu = [h[i]/(h[i] + h[i+1]) for i in 1:n-1]
lambda = [1 - mu[i] for i in 1:n-1]
matrix = diagm(0 => [1; [2*(h[i] + h[i+1]) for i in 1:n-1]; 1]) +
         diagm(1 => vcat(mu, 0)) +
         diagm(-1 => vcat(0, lambda))
rhs = vcat(0, d, 0)
m = matrix \ rhs

# 计算每个子区间的系数
a = [y[i] for i in 1:n]
b = [(y[i+1] - y[i])/h[i] - (2*m[i] + m[i+1])*h[i]/6 for i in 1:n]
c = [m[i]/2 for i in 1:n]
d = [(m[i+1] - m[i])/(6*h[i]) for i in 1:n]

# 构造样条插值函数并输出显式形式
splines = []
for i in 1:n
    push!(splines, (ξ -> a[i] + b[i]*(ξ - x[i]) + c[i]*(ξ - x[i])^2 + d[i]*(ξ - x[i])^3))
    println("样条函数 S$i(ξ) = $(a[i]) + $(b[i])*(ξ - $(x[i])) + $(c[i])*(ξ - $(x[i]))^2 + $(d[i])*(ξ - $(x[i]))^3")
end

# 计算 f(3)
f3 = splines[2](3)

# 输出 f(3)
println("f(3) = ", f3)

# 可视化
scatter(x, y, label="Data Points", color=:red, markersize=8)
plot!(ξ -> splines[1](ξ), 1, 2, label="Spline Interpolation", color=:blue)
plot!(ξ -> splines[2](ξ), 2, 4, color=:blue)
plot!(ξ -> splines[3](ξ), 4, 5, color=:blue)
scatter!([3], [f3], label="f(3)", color=:green, markersize=8)
