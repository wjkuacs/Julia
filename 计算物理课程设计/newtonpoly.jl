using Plots

# 定义数据点
data = [(1, 0), (2, 1), (3, 8), (5, 15), (6, 25)]

function newtonPoly(data)
    n = length(data)
    x_values = [d[1] for d in data]
    y_values = [d[2] for d in data]
    diff = zeros(n, n)
    diff[:, 1] = y_values

    for j in 2:n
        for i in j:n
            diff[i, j] = (diff[i, j-1] - diff[i-1, j-1]) / (x_values[i] - x_values[i-j+1])
        end
    end

    coeffs = [diff[i, i] for i in 1:n]
    return (x) -> begin
        poly = coeffs[1]
        for i in 2:n
            term = coeffs[i]
            for j in 1:i-1
                term *= (x - x_values[j])
            end
            poly += term
        end
        return poly
    end, coeffs, x_values
end

function display_poly(coeffs, x_values)
    n = length(coeffs)
    poly_str = string(coeffs[1])
    for i in 2:n
        term = string(coeffs[i])
        for j in 1:i-1
            term *= "* (x - " * string(x_values[j]) * ")"
        end
        poly_str *= " + " * term
    end
    return poly_str
end

# 生成插值多项式
poly, coeffs, x_values = newtonPoly(data)

# 显示插值多项式
println("插值多项式: ", display_poly(coeffs, x_values))

# 绘制数据点和插值多项式
scatter([d[1] for d in data], [d[2] for d in data], label = "Data Points", color = :red, markersize = 8)
plot!(poly, 0, 7, label = "Interpolation Polynomial", color = :blue)
