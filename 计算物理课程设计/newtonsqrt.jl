using Plots

# 定义牛顿迭代法计算平方根并绘制图形的函数
function NewtonSqrtPlot(x; guess = 10.0, tol = 1.0e-10, maxIter = 100)
    y = guess
    iterValues = Float64[]
    errors = Float64[]
    trueValue = sqrt(x)
    
    # 记录初始猜测值和误差
    push!(iterValues, y)
    push!(errors, abs(y - trueValue))
    
    # 开始迭代
    i = 0
    while abs(y^2 - x) > tol && i < maxIter
        # 更新猜测值
        y = y - (y^2 - x) / (2 * y)
        # 增加迭代次数
        i += 1
        # 记录每次迭代的值和误差
        push!(iterValues, y)
        push!(errors, abs(y - trueValue))
    end
    
    # 绘制迭代值和真实值曲线在同一张图中
    iterationPlot = plot(1:length(iterValues), iterValues, label = "Iterative Values", xlabel = "Iteration", ylabel = "Value", title = "Newton's Method: Iteration Values and True Value")
    plot!(iterationPlot, 1:length(iterValues), [trueValue for _ in iterValues], label = "True Value")
    
    # 绘制误差曲线
    errorPlot = plot(1:length(errors), errors, label = "Error", xlabel = "Iteration", ylabel = "Error", title = "Error in Newton's Method Iteration")
    
    return iterationPlot, errorPlot, errors, y
end

# 计算113的平方根并绘制图形，返回误差
iterationPlot, errorPlot, errors, finalValue = NewtonSqrtPlot(113)

# 分别显示两个图像
display(iterationPlot)
display(errorPlot)

# 显示最终平方根值
println("最终平方根值: ", finalValue)
