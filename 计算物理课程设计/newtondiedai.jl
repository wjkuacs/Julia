using LinearAlgebra

# 初始猜测
x = [0.0, 0.0, 0.0]

# 误差限
epsilon = 0.00000000001

# 定义方程组和雅可比矩阵
function f(x)
    return [12 * x[1] - 3 * x[2] + 3 * x[3] - 15,
            -18 * x[1] + 3 * x[2] - x[3] + 15,
            x[1] + x[2] + x[3] - 6]
end

function j(x)
    return [12 -3 3;
            -18 3 -1;
            1 1 1]
end

# 牛顿迭代法
error = 1.0
iterations = 0

while error > epsilon
    x_old = copy(x)
    x = x - j(x) \ f(x)
    error = maximum(abs.(x - x_old))
    iterations += 1
end

# 显示结果
println("解: ", x)
println("迭代次数: ", iterations)
