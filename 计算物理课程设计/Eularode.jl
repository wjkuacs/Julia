using Plots

# 定义微分方程
function f(x, y)
    return y - (2 * x) / y
end

# 初始条件
y0 = 1.0
x0 = 0.0
xEnd = 1.0
h = 0.1

# 计算步骤数
n = Int((xEnd - x0) / h)

# 初始化列表
xValues = collect(LinRange(x0, xEnd, n+1))
yValues = zeros(n+1)
yValues[1] = y0

# 欧拉法迭代
for i in 1:n
    yValues[i+1] = yValues[i] + h * f(xValues[i], yValues[i])
end

# 解析解求解
function exact_solution(x)
    return sqrt(2 * x + 1)  # 替换为实际解析解
end

# 生成解析解数值
yExactArray = [exact_solution(x) for x in xValues]

# 计算误差
error = yValues - yExactArray

# 显示结果
println("数值解: ", yValues)
println("解析解: ", yExactArray)
println("误差: ", error)

# 绘制图形
plot(xValues, yValues, label="Euler Method Numerical Solution", linestyle=:solid, color=:blue)
plot!(xValues, yExactArray, label="Exact Solution", linestyle=:dash, color=:red)
plot!(xValues, error, label="Error", linestyle=:dot, color=:green)
xlabel!("x")
ylabel!("y")
title!("Euler Method Numerical Solution, Exact Solution, and Error")
