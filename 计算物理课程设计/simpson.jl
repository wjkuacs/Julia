# 导入必要的包
using QuadGK

# 定义函数
function f(x)
    x == 0 ? 1.0 : sin(x) / x
end

# 辛普森方法的实现
function simpson(f, a, b, n)
    h = (b - a) / n
    sum1 = sum(f(a + (2i - 1) * h) for i in 1:(n ÷ 2))
    sum2 = sum(f(a + 2i * h) for i in 1:(n ÷ 2 - 1))
    (h / 3) * (f(a) + 4 * sum1 + 2 * sum2 + f(b))
end

# 误差限
epsilon = 0.0001
n = 2
error = Inf
result = 0.0

# 迭代计算，直到误差小于 0.0001
while error > epsilon
    n *= 2
    result1 = simpson(f, 0, 1, n)
    error = abs(result1 - result)
    result = result1
end

# 显示结果
println("积分结果: ", result)
println("误差: ", error)
println("分区数: ", n)
