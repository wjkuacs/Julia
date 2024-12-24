using LinearAlgebra
using Plots

# 参数设置
L = 1.0                # 位势阱宽度
N = 500                # 网格点数
x = LinRange(0, L, N)  # 空间网格
dx = x[2] - x[1]       # 网格步长
hbar2_2m = 1.0         # 设置 ħ²/(2m) = 1（自然单位制）

# 构造有限差分哈密顿矩阵
H = zeros(N-2, N-2)  # 初始化矩阵大小为 (N-2) x (N-2)

# 主对角线
for i in 1:(N-2)
    H[i, i] = -2.0
end

# 上下对角线
for i in 1:(N-3)
    H[i, i+1] = 1.0
    H[i+1, i] = 1.0
end

# 乘以 -ħ²/(2mΔx²)
H *= -hbar2_2m / dx^2

# 求解本征值和本征向量
eigvals, eigvecs = eigen(H)

# 提取本征值（能量）和本征态
E = eigvals
ψ = eigvecs

# 可视化前几个本征态
p = plot(x[2:end-1], ψ[:, 1], label="n=1, E=$(round(E[1], digits=4))", xlabel="x", ylabel="ψ_n(x)",fontfamily="Microsoft YaHei",
         title="无限深位势阱的本征态", legend=:right)
for n in 2:5
    plot!(p, x[2:end-1], ψ[:, n] .+ E[n], label="n=$n, E=$(round(E[n], digits=4))")
end

# 显示图像
display(p)


##
using LinearAlgebra
using Plots

# 参数设置
L = 1.0                # 位势阱宽度
N = 500                # 网格点数
x = LinRange(0, L, N)  # 空间网格
dx = x[2] - x[1]       # 网格步长
hbar2_2m = 1.0         # 设置 ħ²/(2m) = 1（自然单位制）

# 构造有限差分哈密顿矩阵
H = zeros(N-2, N-2)  # 初始化矩阵大小为 (N-2) x (N-2)

# 主对角线
for i in 1:(N-2)
    H[i, i] = -2.0
end

# 上下对角线
for i in 1:(N-3)
    H[i, i+1] = 1.0
    H[i+1, i] = 1.0
end

# 乘以 -ħ²/(2mΔx²)
H *= -hbar2_2m / dx^2

# 求解本征值和本征向量
eigvals, eigvecs = eigen(H)

# 提取本征值（能量）和本征态
E = eigvals
ψ = eigvecs

# 可视化前几个本征态（波函数）
p = plot(x[2:end-1], ψ[:, 1], label="n=1, E=$(round(E[1], digits=4))", xlabel="x", ylabel="ψ_n(x)",fontfamily="Microsoft YaHei",title="无限深位势阱的波函数", legend=:right)
for n in 2:5
    plot!(p, x[2:end-1], ψ[:, n], label="n=$n, E=$(round(E[n], digits=4))")
end

# 显示图像
display(p)

# 如果需要将图像保存为文件，可以使用以下命令：
# savefig(p, "wavefunction_plot.png")



##
using Roots
using Distributions
using PyCall
using PyPlot