using Plots

# Integral of sin(x) from 0 to π
function f(x)
    return sin(x)
end

function F(x)
    return -cos(x)
end

x0 = 0
xf = pi

I_exact = F(xf) - F(x0)

function trap(f, x0, xf, n)
    h = (xf - x0) / n
    
    I = 0.5 * (f(x0) + f(xf)) # add end points
    
    for i = 1:n-1
        I += f(x0 + i * h)
    end
    return I * h
end

function simp(f, x0, xf, n)
    h = (xf - x0) / n
    
    I = (f(x0) + f(xf)) # add end points
    
    for i = 1:n-1 # n must be even
        if mod(i, 2) == 1
            I += 4 * f(x0 + i * h)
        else
            I += 2 * f(x0 + i * h)
        end
    end
    return I * h / 3
end

# Loop over h values and compute trap/simps integrals, h = (xf - x0) / n
I_trap = []
I_simp = []
hs = []

for n = 4:2:30 # only even n
    push!(I_trap, trap(f, x0, xf, n))
    push!(I_simp, simp(f, x0, xf, n))
    push!(hs, (xf - x0) / n)
end

# Scaling of Error
plot(log.(hs), log.(abs.(I_exact .- I_trap)), markershape = :circle, markersize = 2, label = "Trapezoid")
plot!(log.(hs), log.(abs.(I_exact .- I_simp)), xlabel = "log(h)", ylabel = "log(∆I)", markershape = :circle, markersize = 2, label = "Simpson's")

# Calculate slope of log-log plot
∆x = log(hs[end]) - log(hs[1])

trap_∆y = log(abs(I_exact - I_trap[end])) - log(abs(I_exact - I_trap[1]))
simp_∆y = log(abs(I_exact - I_simp[end])) - log(abs(I_exact - I_simp[1]))

println("Trapezoid slope = ", trap_∆y / ∆x)
println("Simpson's slope = ", simp_∆y / ∆x)

# Simple Romberg Integration
println("h_0 = ", hs[1], ", I_0 = ", I_trap[1])
println("h_1 = ", hs[3], ", I_1 = ", I_trap[3])
println("h_1/h_0 = ", hs[3] / hs[1])
romberg_simple = (4/3) * I_trap[3] - (1/3) * I_trap[1]
println()
println("Simple Romberg approximation (m = 1): I = ", romberg_simple)

p = plot(hs, I_trap, xlabel = "h", ylabel = "I", markershape = :circle, markersize = 2, label = "Trapezoid")
plot!(hs, I_simp, markershape = :circle, markersize = 2, label = "Simpson's")

plot!([hs[1], hs[end]], [romberg_simple, romberg_simple], linestyle = :dash, c = :black, label = "Simple Romberg", legend = :right)

plot!(
    hs[5:end], I_trap[5:end], markershape = :circle, markersize = 2, label = "", 
    inset = (1, bbox(0.1, 0.15, 0.4, 0.4, :bottom, :left)), 
    xticks = round.(hs[5:3:end], digits = 2), 
    yticks = [1.997, 2], 
    subplot = 2
)

plot!(p[2], hs[5:end], I_simp[5:end], markershape = :circle, markersize = 2, label = "")
plot!(p[2], [hs[5], hs[end]], [romberg_simple, romberg_simple], linestyle = :dash, c = :black, label = "")

# Integration using QuadGK
using QuadGK

quadgk(f, x0, xf)

function g(x)
    return log(cos(exp(x))^2)  # 1/(x-pi/2)^2
end

plot(x0:0.01:xf, g.(x0:0.01:xf))

I_QGK, σ_QGK = quadgk(g, x0, xf)

# Loop over h values and compute trap/simps integrals, h = (xf - x0) / n
I_trap = []
I_simp = []
hs = []

for n = 400:200:5000 # only even n
    push!(I_trap, trap(g, x0, xf, n))
    push!(I_simp, simp(g, x0, xf, n))
    push!(hs, (xf - x0) / n)
end

i, j = 4, 9
println("h_0 = ", hs[i], ", I_0 = ", I_trap[i])
println("h_1 = ", hs[j], ", I_1 = ", I_trap[j])
println("h_1/h_0 = ", hs[j] / hs[i])
romberg_simple = (4/3) * I_trap[j] - (1/3) * I_trap[i]
println()
println("Simple Romberg approximation (m = 1): I = ", romberg_simple)

p = plot(hs, I_trap, xlabel = "h", ylabel = "I", markershape = :circle, markersize = 2, label = "Trapezoid")
plot!(hs, I_simp, markershape = :circle, markersize = 2, label = "Simpson's")

plot!([hs[1], hs[end]], [romberg_simple, romberg_simple], linestyle = :dash, c = :black, label = "Simple Romberg")
plot!([hs[1], hs[end]], [I_QGK, I_QGK], linestyle = :dash, c = :red, label = "QuadGK", legend = :right)
