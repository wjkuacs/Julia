### 这段代码是为了学习DifferentialEquations.jl宏包而创建的
using DifferentialEquations
using Plots

## 求解标量方程

# 定义问题
f(u,p,t) = 1.01 * u
u_0 = 1/2
t_span = (0.0 , 1.0)
prob = ODEProblem(f, u_0, t_span)
# 解决问题，使用solve函数
sol = solve(prob, BS3(),reltol = 1e-8)
plot(sol, linewidth = 5, title = "Solution to the linear ODE with a thick line",xaxis = "Time (t)", yaxis = "u(t) (in μm)", label = "My Thick Line!") # legend=false
plot!(sol.t, t -> 0.5 * exp(1.01t), lw = 3, ls = :dash, label = "True Solution!")

sol[5]

sol.t[8]





###

## 求解方程组

# 定义一下 Lorenz 方程
function lorenz!(du,u,p,t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) -u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan)
sol = solve(prob)

plot(sol, idxs = (1, 2, 3))
