using MosaicViews
using Images, TestImages

function calculate_averages(img::AbstractArray{T, N}, Hð½::AbstractArray{T, M}) where {T<:Real, N, M}
    Hð½â± = @. 1. - Hð½
    â«Hð½ = sum(Hð½)
    â«Hð½â± = sum(Hð½â±)
    if ndims(img) == 2
        â«uâHð½ = sum(img .* Hð½)
        â«uâHð½â± = sum(img .* Hð½â±)
    elseif ndims(img) == 3
        â«uâHð½ = sum(img .* Hð½, dims=(1, 2))
        â«uâHð½â± = sum(img .* Hð½â±, dims=(1, 2))
    end
    if â«Hð½ != 0
        câ = â«uâHð½ / â«Hð½
    end
    if â«Hð½â± != 0
        câ = â«uâHð½â± / â«Hð½â±
    end

    return câ, câ
end

function difference_from_average_term(img::AbstractArray{T, N}, Hð½::AbstractArray{T, M}, Î»â::Float64, Î»â::Float64) where {T<:Real, N, M}
    câ, câ = calculate_averages(img, Hð½)

    if ndims(img) == 2
        return @. -Î»â * (img - câ)^2 + Î»â * (img - câ)^2
    elseif ndims(img) == 3
        return -Î»â .* sum((img .- câ).^2, dims=3) .+ Î»â .* sum((img .- câ).^2, dims=3)
    end
end

function Hâ(x::AbstractArray{T,N}, h::Float64=1.0) where {T<:Real, N}
    return @. 1. / 2. * (1. + 2. / pi * atan(x / h))
end

function Î´â(x::AbstractArray{T,N}, h::Float64=1.0) where {T<:Real, N}
    return @. h / (h^2 + x^2)
end
    
function initial_level_set(shape::Tuple)
    xâ = reshape(collect(0:shape[begin]-1), shape[begin], 1)
    yâ = reshape(collect(0:shape[begin+1]-1), 1, shape[begin+1])
    ð½â = @. sin(pi / 5 * xâ) * sin(pi / 5 * yâ)
end

function calculate_variation(img::AbstractArray{T, N}, ð½â¿::AbstractArray{T, M}, Î¼::Float64, Î»â::Float64, Î»â::Float64, Ît::Float64) where {T<:Real, N, M}
    Ïµ = 1e-16
    ð½âº = padarray(ð½â¿, Pad(1, 1))
    r = axes(ð½âº)
    @views ð½áµ¢â = ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 2:last(r[2])] .- ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1]
    @views ð½áµ¢â = ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1] .- ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]):last(r[2]) - 2]
    @views ð½áµ¢  = (ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 2:last(r[2])] .- ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]):last(r[2]) - 2]) / 2.

    @views ð½â±¼â = ð½âº[first(r[1]) + 2:last(r[1]), first(r[2]) + 1:last(r[2]) - 1] .- ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1]
    @views ð½â±¼â = ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1] .- ð½âº[first(r[1]):last(r[1]) - 2, first(r[2]) + 1:last(r[2]) - 1]
    @views ð½â±¼  = (ð½âº[first(r[1]) + 2:last(r[1]), first(r[2]) + 1:last(r[2]) - 1] .- ð½âº[first(r[1]):last(r[1]) - 2, first(r[2]) + 1:last(r[2]) - 1]) / 2.

    Câ = @. 1. / sqrt(Ïµ + ð½áµ¢â^2 + ð½â±¼^2)
    Câ = @. 1. / sqrt(Ïµ + ð½áµ¢â^2 + ð½â±¼^2)
    Câ = @. 1. / sqrt(Ïµ + ð½áµ¢^2 + ð½â±¼â^2)
    Câ = @. 1. / sqrt(Ïµ + ð½áµ¢^2 + ð½â±¼â^2)

    @views K = ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 2:last(r[2])] .* Câ .+
               ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]):last(r[2]) - 2] .* Câ .+
               ð½âº[first(r[1]) + 2:last(r[1]), first(r[2]) + 1:last(r[2]) - 1] .* Câ .+
               ð½âº[first(r[1]):last(r[1]) - 2, first(r[2]) + 1:last(r[2]) - 1] .* Câ
    Hð½ = @. 1. * (ð½â¿ > 0)

    diff = difference_from_average_term(img, Hð½, Î»â, Î»â)
    ð½â¿âºÂ¹ = (ð½â¿ .+ Ît .* Î´â(ð½â¿) .* (Î¼ * K .+ diff)) ./ (1 .+ Î¼ .* Ît .* Î´â(ð½â¿) .* (Câ .+ Câ .+ Câ .+ Câ))
end

function calculate_reinitial(ð½::AbstractArray{T, M}, Ît::Float64) where {T<:Real, M}
    Ïµ = 1e-8
    ð½âº = padarray(ð½, Pad(1, 1))
    r = axes(ð½âº)
    @views a = @. ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1] - ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]):last(r[2]) - 2]
    @views b = @. ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 2:last(r[2])] - ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1]
    @views c = @. ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1] - ð½âº[first(r[1]):last(r[1]) - 2, first(r[2]) + 1:last(r[2]) - 1]
    @views d = @. ð½âº[first(r[1]) + 2:last(r[1]), first(r[2]) + 1:last(r[2]) - 1] - ð½âº[first(r[1]) + 1:last(r[1]) - 1, first(r[2]) + 1:last(r[2]) - 1]

    aâº = max.(a, 0)
    aâ» = min.(a, 0)
    bâº = max.(b, 0)
    bâ» = min.(b, 0)
    câº = max.(c, 0)
    câ» = min.(c, 0)
    dâº = max.(d, 0)
    dâ» = min.(d, 0)

    G = zeros(size(ð½))
    indexâº = ð½ .> 0
    indexâ» = ð½ .< 0
    @. G = (sqrt(max(aâº^2, bâ»^2) + max(câº^2, dâ»^2)) - 1) * indexâº + (sqrt(max(aâ»^2, bâº^2) + max(câ»^2, dâº^2)) - 1) * indexâ» 
    signð½ = @. ð½ / sqrt(ð½^2 + Ïµ)
    ð¿ = @. ð½ - Ît * signð½ * G
end

function reinitialize(ð½::AbstractArray{T, M}, Ît::Float64, max_reiter::Int64=5) where {T<:Real, M}
    iter = 0
    while iter < max_reiter
        ð½ .= calculate_reinitial(ð½, Ît)
        iter += 1
    end

    return ð½
end

function chan_vese(img::AbstractArray{T,N};
                   Î¼::Float64=0.25,
                   Î»â::Float64=1.0,
                   Î»â::Float64=1.0,
                   tol::Float64=1e-3,
                   max_iter::Int64=500,
                   Ît::Float64=0.5,
                   reinitial_flag::Bool=false) where {T<:Real, N}
    iter = 0
    D = ndims(img)
    if D == 3
        img = PermutedDimsArray(img, (2, 3, 1))
    end
    ð½â¿ = initial_level_set(size(img))
    Î´ = tol + 1
    img .= img .- minimum(img)

    if maximum(img) != 0
        img .= img ./ maximum(img)
    end

    while (Î´ > tol) & (iter < max_iter)
        ð½â¿âºÂ¹ = calculate_variation(img, ð½â¿, Î¼, Î»â, Î»â, Ît)
        Î´ = sqrt(meanfinite((ð½â¿âºÂ¹ .- ð½â¿).^2, (1, 2))[1])
        if reinitial_flag
            ð½â¿ .= reinitialize(ð½â¿âºÂ¹, Ît)
        else
            if D == 3
                r = axes(ð½â¿âºÂ¹)
                @views ð½â¿ .= ð½â¿âºÂ¹[:, :, first(r[3])]
            else
                @views ð½â¿ .= ð½â¿âºÂ¹
            end
        end

        iter += 1
    end

    return ð½â¿, iter
end

img = float64.(channelview(testimage("cameraman")))
# img = float64.(channelview(testimage("lena_color_512")))
ð½, iter_num = chan_vese(img, Î¼=0.25, Î»â=1.0, Î»â=1.0, tol=1e-3, max_iter=200, Ît=0.5, reinitial_flag=false)

segmentation = ð½ .> 0
print(iter_num)
ð½ .= ð½ .- minimum(ð½)

if maximum(ð½) != 0
    ð½ .= ð½ ./ maximum(ð½)
end

output = mosaicview(img, segmentation, ð½; nrow=1, ncol=3, rowmajor=true)
save("demo.png", output)
