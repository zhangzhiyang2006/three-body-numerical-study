# Plot a complete relative two-body orbit for the assignment's initial conditions.
using LinearAlgebra, Plots, Printf
ENV["GKSwstype"] = "100"
gr()
m=[0.3,0.03]
y0=[2.0,2.0,0.0,0.0,0.2,-0.2,-0.01,0.01]
function rhs(y)
    dx=y[3]-y[1]; dy=y[4]-y[2]
    r2=dx^2+dy^2; q=1/(r2*sqrt(r2))
    [y[5],y[6],y[7],y[8],m[2]*q*dx,m[2]*q*dy,-m[1]*q*dx,-m[1]*q*dy]
end
function step(y,h)
    k1=rhs(y);k2=rhs(y .+ h/2 .* k1);k3=rhs(y .+ h/2 .* k2);k4=rhs(y .+ h .* k3)
    y .+ h/6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end
function solve(h,T)
    n=round(Int,T/h);Y=Matrix{Float64}(undef,8,n+1);Y[:,1]=y0
    for k in 1:n;Y[:,k+1]=step(Y[:,k],h);end
    Y
end
T=40.0
Y=solve(0.001,T);Yhalf=solve(0.0005,T)
@assert maximum(abs.(Y-Yhalf[:,1:2:end]))<1e-8
com=(m[1].*Y[1:2,:].+m[2].*Y[3:4,:])./sum(m)
R=copy(Y[1:4,:]);R[1:2,:].-=com;R[3:4,:].-=com
p=plot(aspect_ratio=:equal, size=(900,650),legend=:topright,grid=true,
    framestyle=:box,xlabel="x - Xcm",ylabel="y - Ycm",title="Two-body center-of-mass frame: 0 ≤ t ≤ 40",
    guidefontsize=11,tickfontsize=9,legendfontsize=9,titlefontsize=12)
colors=[:royalblue,:darkorange]
for i in 1:2
    x=R[2i-1,:];y=R[2i,:]
    plot!(p,x,y,color=colors[i],lw=2,label="body $i")
    scatter!(p,[x[1]],[y[1]],marker=:circle,color=colors[i],ms=6,label="")
    scatter!(p,[x[20001]],[y[20001]],marker=:star5,color=colors[i],ms=7,label="")
    scatter!(p,[x[end]],[y[end]],marker=:diamond,color=colors[i],ms=6,label="")
end
scatter!(p,[NaN],[NaN],marker=:circle,color=:black,label="t=0")
scatter!(p,[NaN],[NaN],marker=:star5,color=:black,label="t=20")
scatter!(p,[NaN],[NaN],marker=:diamond,color=:black,label="t=40")
savefig(p,"figures/two_com.png")
@printf("Two-body T=40 RK4 step-half max state difference = %.6e\n",maximum(abs.(Y-Yhalf[:,1:2:end])))
