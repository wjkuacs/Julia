using Plots

# 定义初始条件和参数
Q0 = 1.0  # 初始电荷量
tMax = 5.0  # 最大时间
h = 0.01  # 步长

# 实现欧拉方法
function euler_method(Q0, tMax, h)
    Q = Q0
    t = 0.0
    Qlist = [(t, Q)]  # 初始化列表时记录初始条件
    while t <= tMax
        Q += h * (-Q)  # 使用欧拉方法更新Q
        t += h  # 更新时间
        push!(Qlist, (t, Q))  # 记录每一步的Q和t
    end
    return Qlist
end

# 使用欧拉方法求解
numerical_solution = euler_method(Q0, tMax, h)

# 定义解析解
analytical_solution(t) = Q0 * exp(-t)

# 定义误差函数
function error_function(numerical_solution, analytical_solution, tMax, h)
    error_list = []
    for t in 0:h:tMax
        # 查找数值解中对应时间点的值
        idx = Int(round(t / h)) + 1  # 将时间 t 转换为索引，+1 是因为 Julia 的索引从 1 开始
        Q_numerical = numerical_solution[idx][2]  # 获取对应时间点的数值解
        Q_analytical = analytical_solution(t)  # 解析解
        push!(error_list, (t, abs(Q_analytical - Q_numerical)))
    end
    return error_list
end

# 计算误差
error_list = error_function(numerical_solution, analytical_solution, tMax, h)

# 提取绘图数据
t_numerical = [p[1] for p in numerical_solution]
Q_numerical = [p[2] for p in numerical_solution]
t_error = [p[1] for p in error_list]
Q_error = [p[2] for p in error_list]

# 绘制数值解、解析解和误差曲线
plot(t_numerical, Q_numerical, label="Numerical Solution", color=:red, linewidth=2)
plot!(0:h:tMax, [analytical_solution(t) for t in 0:h:tMax], label="Analytical Solution", color=:blue, linewidth=2)
plot!(t_error, Q_error, label="Error Curve", color=:green, linewidth=2)
xlabel!("t")
ylabel!("Q")
