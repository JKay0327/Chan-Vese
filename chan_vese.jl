using MosaicViews
using Images, TestImages

function calculate_averages(img::AbstractArray{T, N}, H𝚽::AbstractArray{T, M}) where {T<:Real, N, M}
    H𝚽ⁱ = @. 1. - H𝚽
    ∫H𝚽 = sum(H𝚽)
    ∫H𝚽ⁱ = sum(H𝚽ⁱ)
    if ndims(img) == 2
        ∫u₀H𝚽 = sum(img .* H𝚽)
        ∫u₀H𝚽ⁱ = sum(img .* H𝚽ⁱ)
    elseif ndims(img) == 3
        ∫u₀H𝚽 = sum(img .* H𝚽, dims=(1, 2))
        ∫u₀H𝚽ⁱ = sum(img .* H𝚽ⁱ, dims=(1, 2))
    end
    if ∫H𝚽 != 0
        c₁ = ∫u₀H𝚽 / ∫H𝚽
    end
    if ∫H𝚽ⁱ != 0
        c₂ = ∫u₀H𝚽ⁱ / ∫H𝚽ⁱ
    end

    return c₁, c₂
end

function difference_from_average_term(img::AbstractArray{T, N}, H𝚽::AbstractArray{T, M}, λ₁::Float64, λ₂::Float64) where {T<:Real, N, M}
    c₁, c₂ = calculate_averages(img, H𝚽)

    if ndims(img) == 2
        return @. -λ₁ * (img - c₁)^2 + λ₂ * (img - c₂)^2
    elseif ndims(img) == 3
        return -λ₁ .* sum((img .- c₁).^2, dims=3) .+ λ₂ .* sum((img .- c₂).^2, dims=3)
    end
end

function Hₕ(x::AbstractArray{T,N}, h::Float64=1.0) where {T<:Real, N}
    return @. 1. / 2. * (1. + 2. / pi * atan(x / h))
end

function δₕ(x::AbstractArray{T,N}, h::Float64=1.0) where {T<:Real, N}
    return @. h / (h^2 + x^2)
end
    
function initial_level_set(shape::Tuple)
    x₀ = reshape(collect(0:shape[begin]-1), shape[begin], 1)
    y₀ = reshape(collect(0:shape[begin+1]-1), 1, shape[begin+1])
    𝚽₀ = @. sin(pi / 5 * x₀) * sin(pi / 5 * y₀)
end

function calculate_variation(img::AbstractArray{T, N}, 𝚽ⁿ::AbstractArray{T, M}, μ::Float64, λ₁::Float64, λ₂::Float64, Δt::Float64) where {T<:Real, N, M}
    ϵ = 1e-16
    𝚽⁺ = padarray(𝚽ⁿ, Pad(1, 1))
    r = axes(𝚽⁺)
    @views 𝚽ᵢ₊ = 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 2:last(r[2])] .- 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1]
    @views 𝚽ᵢ₋ = 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1] .- 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]):last(r[2]) - 2]
    @views 𝚽ᵢ  = (𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 2:last(r[2])] .- 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]):last(r[2]) - 2]) / 2.

    @views 𝚽ⱼ₊ = 𝚽⁺[first(r[1]) + 2:last(r[1]), first(r[2]) + 1:last(r[2]) - 1] .- 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1]
    @views 𝚽ⱼ₋ = 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1] .- 𝚽⁺[first(r[1]):last(r[1]) - 2, first(r[2]) + 1:last(r[2]) - 1]
    @views 𝚽ⱼ  = (𝚽⁺[first(r[1]) + 2:last(r[1]), first(r[2]) + 1:last(r[2]) - 1] .- 𝚽⁺[first(r[1]):last(r[1]) - 2, first(r[2]) + 1:last(r[2]) - 1]) / 2.

    C₁ = @. 1. / sqrt(ϵ + 𝚽ᵢ₊^2 + 𝚽ⱼ^2)
    C₂ = @. 1. / sqrt(ϵ + 𝚽ᵢ₋^2 + 𝚽ⱼ^2)
    C₃ = @. 1. / sqrt(ϵ + 𝚽ᵢ^2 + 𝚽ⱼ₊^2)
    C₄ = @. 1. / sqrt(ϵ + 𝚽ᵢ^2 + 𝚽ⱼ₋^2)

    @views K = 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 2:last(r[2])] .* C₁ .+
               𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]):last(r[2]) - 2] .* C₂ .+
               𝚽⁺[first(r[1]) + 2:last(r[1]), first(r[2]) + 1:last(r[2]) - 1] .* C₃ .+
               𝚽⁺[first(r[1]):last(r[1]) - 2, first(r[2]) + 1:last(r[2]) - 1] .* C₄
    H𝚽 = @. 1. * (𝚽ⁿ > 0)

    diff = difference_from_average_term(img, H𝚽, λ₁, λ₂)
    𝚽ⁿ⁺¹ = (𝚽ⁿ .+ Δt .* δₕ(𝚽ⁿ) .* (μ * K .+ diff)) ./ (1 .+ μ .* Δt .* δₕ(𝚽ⁿ) .* (C₁ .+ C₂ .+ C₃ .+ C₄))
end

function calculate_reinitial(𝚽::AbstractArray{T, M}, Δt::Float64) where {T<:Real, M}
    ϵ = 1e-8
    𝚽⁺ = padarray(𝚽, Pad(1, 1))
    r = axes(𝚽⁺)
    @views a = @. 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1] - 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]):last(r[2]) - 2]
    @views b = @. 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 2:last(r[2])] - 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1]
    @views c = @. 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1] - 𝚽⁺[first(r[1]):last(r[1]) - 2, first(r[2]) + 1:last(r[2]) - 1]
    @views d = @. 𝚽⁺[first(r[1]) + 2:last(r[1]), first(r[2]) + 1:last(r[2]) - 1] - 𝚽⁺[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1]

    a⁺ = max.(a, 0)
    a⁻ = min.(a, 0)
    b⁺ = max.(b, 0)
    b⁻ = min.(b, 0)
    c⁺ = max.(c, 0)
    c⁻ = min.(c, 0)
    d⁺ = max.(d, 0)
    d⁻ = min.(d, 0)

    G = zeros(size(𝚽))
    index⁺ = 𝚽 .> 0
    index⁻ = 𝚽 .< 0
    @. G = (sqrt(max(a⁺^2, b⁻^2) + max(c⁺^2, d⁻^2)) - 1) * index⁺ + (sqrt(max(a⁻^2, b⁺^2) + max(c⁻^2, d⁺^2)) - 1) * index⁻ 
    sign𝚽 = @. 𝚽 / sqrt(𝚽^2 + ϵ)
    𝚿 = @. 𝚽 - Δt * sign𝚽 * G
end

function reinitialize(𝚽::AbstractArray{T, M}, Δt::Float64, max_reiter::Int64=5) where {T<:Real, M}
    iter = 0
    while iter < max_reiter
        𝚽 .= calculate_reinitial(𝚽, Δt)
        iter += 1
    end

    return 𝚽
end

function chan_vese(img::AbstractArray{T,N};
                   μ::Float64=0.25,
                   λ₁::Float64=1.0,
                   λ₂::Float64=1.0,
                   tol::Float64=1e-3,
                   max_iter::Int64=500,
                   Δt::Float64=0.5,
                   reinitial_flag::Bool=false) where {T<:Real, N}
    iter = 0
    D = ndims(img)
    if D == 3
        img = PermutedDimsArray(img, (2, 3, 1))
    end
    𝚽ⁿ = initial_level_set(size(img))
    δ = tol + 1
    img .= img .- minimum(img)

    if maximum(img) != 0
        img .= img ./ maximum(img)
    end

    while (δ > tol) & (iter < max_iter)
        𝚽ⁿ⁺¹ = calculate_variation(img, 𝚽ⁿ, μ, λ₁, λ₂, Δt)
        δ = sqrt(meanfinite((𝚽ⁿ⁺¹ .- 𝚽ⁿ).^2, (1, 2))[1])
        if reinitial_flag
            𝚽ⁿ .= reinitialize(𝚽ⁿ⁺¹, Δt)
        else
            if D == 3
                r = axes(𝚽ⁿ⁺¹)
                @views 𝚽ⁿ .= 𝚽ⁿ⁺¹[:, :, first(r[3])]
            else
                @views 𝚽ⁿ .= 𝚽ⁿ⁺¹
            end
        end

        iter += 1
    end

    return 𝚽ⁿ, iter
end

img = float64.(channelview(testimage("cameraman")))
# img = float64.(channelview(testimage("lena_color_512")))
𝚽, iter_num = chan_vese(img, μ=0.25, λ₁=1.0, λ₂=1.0, tol=1e-3, max_iter=200, Δt=0.5, reinitial_flag=false)

segmentation = 𝚽 .> 0
print(iter_num)
𝚽 .= 𝚽 .- minimum(𝚽)

if maximum(𝚽) != 0
    𝚽 .= 𝚽 ./ maximum(𝚽)
end

output = mosaicview(img, segmentation, 𝚽; nrow=1, ncol=3, rowmajor=true)
save("demo.png", output)
