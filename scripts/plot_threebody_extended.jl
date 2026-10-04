# Reproduce the prescribed three-body initial-value problem and show its complete plotted interval.
using LinearAlgebra, Plots, Printf
ENV["GKSwstype"] = "100"
gr()

m = [0.3, 0.03, 0.03]
y0 = [2.0, 2.0, 0.0, 0.0, -2.0, -2.0, 0.2, -0.2, 0.0, 0.0, -0.2, 0.2]
function rhs(y)
    a = zeros(6)
    for i in 1:3, j in 1:3
        i == j && continue
        dx = y[2j-1] - y[2i-1]
        dy = y[2j] - y[2i]
        r2 = dx^2 + dy^2
        c = m[j] / (r2 * sqrt(r2))
        a[2i-1] += c * dx
        a[2i] += c * dy
    end
    [y[7:12]; a]
end
function rk4(y, h)
    k1 = rhs(y)
    k2 = rhs(y .+ (h/2) .* k1)
    k3 = rhs(y .+ (h/2) .* k2)
    k4 = rhs(y .+ h .* k3)
    y .+ (h/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end
function solve(h, T)
    n = round(Int, T/h)
    y = copy(y0)
    Y = Matrix{Float64}(undef, 12, n+1)
    Y[:,1] = y
    for k in 1:n
        y = rk4(y,h)
        Y[:,k+1] = y
    end
    Y
end

T = 200.0
Y = solve(0.001, T)
Yhalf = solve(0.0005, T)
@assert maximum(abs.(Y - Yhalf[:,1:2:end])) < 1e-8
# Verify that the original t=30 result has not changed.
expected30 = [7.073899882375271,-3.4947616887680986,6.105001404608125,-5.8841162255751644,-4.8440002283608194,4.831733113255898]
@assert maximum(abs.(Y[1:6,30001] - expected30)) < 1e-9
com = (m[1] .* Y[1:2,:] .+ m[2] .* Y[3:4,:] .+ m[3] .* Y[5:6,:]) ./ sum(m)
R = copy(Y[1:6,:])
for i in 1:3
    R[2i-1:2i,:] .-= com
end
colors = [:royalblue, :darkorange, :seagreen]
function panel(lastindex, titletext)
    p = plot(aspect_ratio=:equal, size=(1000,480), legend=:topright,
        xlabel="x - Xcm", ylabel="y - Ycm", title=titletext,
        grid=true, framestyle=:box, guidefontsize=10, tickfontsize=8,
        legendfontsize=8, titlefontsize=11)
    for i in 1:3
        x=R[2i-1,1:lastindex]; y=R[2i,1:lastindex]
        plot!(p,x,y,color=colors[i],lw=2,label="body $i")
        scatter!(p,[x[1]],[y[1]],marker=:circle,ms=5,color=colors[i],label="")
        scatter!(p,[x[end]],[y[end]],marker=:diamond,ms=6,color=colors[i],label="")
    end
    p
end
p1=panel(size(Y,2),"Full computed interval: 0 ≤ t ≤ 200")
p2=panel(30001,"Initial segment: 0 ≤ t ≤ 30")
plot(p1,p2,layout=(1,2),size=(1500,600),margin=5Plots.mm)
mkpath("figures")
savefig("figures/three_com.png")
# Pair distances over the same plotted interval; sample the dense integration only for display.
idx = 1:100:size(Y,2)
t = (idx .- 1) .* 0.001
pdist = plot(xlabel="time", ylabel="pair distance", title="Pair distances: 0 ≤ t ≤ 200",
    size=(1000,480), legend=:topleft, grid=true, framestyle=:box,
    guidefontsize=10, tickfontsize=8, legendfontsize=9)
for (k,(i,j)) in enumerate(((1,2),(1,3),(2,3)))
    dist = [norm(Y[2i-1:2i,z]-Y[2j-1:2j,z]) for z in idx]
    plot!(pdist,t,dist,color=colors[k],lw=2,label="r$(i)$(j)")
end
vline!(pdist,[30.0],color=:gray,ls=:dash,lw=1,label="original T=30")
savefig(pdist,"figures/three_distances.png")
@printf("T=200, RK4 step-half max state difference = %.6e\n", maximum(abs.(Y-Yhalf[:,1:2:end])))
for (i,j) in ((1,2),(1,3),(2,3))
    d=norm(Y[2i-1:2i,end]-Y[2j-1:2j,end])
    @printf("r%d%d(200) = %.6f\n",i,j,d)
end
