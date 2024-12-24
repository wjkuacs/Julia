using LinearAlgebra

# 定义系数矩阵 A 和常数向量 b
A = [12.0 -3.0 3.0; -18.0 3.0 -1.0; 1.0 1.0 1.0]
b = [15.0, -15.0, 6.0]

# 将矩阵 A 和 b 组合成增广矩阵
augmentedMatrix = hcat(A, b)

# 获取增广矩阵的行数
n = size(augmentedMatrix, 1)

# 前向消元：使用列主元法
for i in 1:n
    # 找到当前列中绝对值最大的元素所在的行
    maxRow = argmax(abs.(augmentedMatrix[i:end, i]))[1] + i - 1
    # 如果最大元素不是当前行，交换当前行和主元行
    if maxRow != i
        augmentedMatrix[[i, maxRow], :] = augmentedMatrix[[maxRow, i], :]
    end
    # 对当前行以下的每一行进行消元
    for j in i+1:n
        factor = augmentedMatrix[j, i] / augmentedMatrix[i, i]
        augmentedMatrix[j, :] .= augmentedMatrix[j, :] .- factor * augmentedMatrix[i, :]
    end
end

# 回代求解
x = zeros(n)
for i in n:-1:1
    x[i] = (augmentedMatrix[i, end] - sum(augmentedMatrix[i, i+1:end-1] .* x[i+1:end])) / augmentedMatrix[i, i]
end

# 输出解
x
