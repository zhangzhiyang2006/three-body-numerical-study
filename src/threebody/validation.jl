# Diagnostics shared by the notebook and regression tests.
using LinearAlgebra, DelimitedFiles, Printf

function parabolic_state(t)
    D=2sinh(asinh(3t/8)/3)
    [4D, 2-2D^2, 1/(1+D^2), -D/(1+D^2)]
end

function centers(Y, masses)
    N=length(masses)
    sum(masses[i].*Y[2i-1:2i,:] for i in 1:N)/sum(masses)
end

function relative_positions(Y, masses)
    N=length(masses); C=centers(Y,masses)
    reduce(vcat,[Y[2i-1:2i,:]-C for i in 1:N])
end

function orbit_diagnostics(t,Y,masses;G=1.0)
    N=length(masses); nt=length(t)
    P=sum(masses[i].*Y[2N+2i-1:2N+2i,:] for i in 1:N)
    C=centers(Y,masses)
    L=zeros(nt);E=zeros(nt);rmin=fill(Inf,nt)
    for k in 1:nt
        for i in 1:N
            x,y=Y[2i-1:2i,k];vx,vy=Y[2N+2i-1:2N+2i,k]
            E[k]+=masses[i]*(vx^2+vy^2)/2
            L[k]+=masses[i]*(x*vy-y*vx)
            for j in i+1:N
                r=norm(Y[2i-1:2i,k]-Y[2j-1:2j,k])
                E[k]-=G*masses[i]*masses[j]/r
                rmin[k]=min(rmin[k],r)
            end
        end
    end
    predicted=C[:,1].+(P[:,1]/sum(masses))*(t.-t[1])'
    (;energy=E,angular=L,momentum=P,center=C,min_distance=rmin,
      energy_error=maximum(abs.(E.-E[1])),angular_error=maximum(abs.(L.-L[1])),
      momentum_error=maximum(abs.(P.-P[:,1])),cm_error=maximum(abs.(C.-predicted)),
      internal_energy=E[1]-sum(abs2,P[:,1])/(2sum(masses)))
end

function reversal_error(f,y0,T;h=.001)
    N=length(y0)÷4
    _,Y=integrate(f,y0,(0.,T);h=h)
    z=copy(Y[:,end]);z[2N+1:end].*=-1
    _,Z=integrate(f,z,(0.,T);h=h)
    zback=copy(Z[:,end]);zback[2N+1:end].*=-1
    maximum(abs.(zback-y0))
end

function read_reference(name)
    path=joinpath("data","threebody_reference",name*".csv")
    isfile(path) || error("缺少独立参考数据 $(path)；请在项目根目录运行笔记，或先运行 scripts/threebody_reference.py")
    a=readdlm(path,',',Float64)
    a[:,1],permutedims(a[:,2:end])
end

function aligned_solution(f,y0,tref;h)
    T=tref[end]-tref[1]
    t,Y=integrate(f,y0,(tref[1],tref[end]);h=h)
    indices=round.(Int,(tref.-tref[1])./(t[2]-t[1])).+1
    maximum(abs.(t[indices]-tref)) < 1e-9 || error("比较网格未对齐")
    t,Y,indices
end

function basic_validation()
    cases=[("single",[1.],[0.,2.,1.,0.],50.),
           ("two",[.3,.03],[2.,2.,0.,0.,.2,-.2,-.01,.01],20.),
           ("three",[.3,.03,.03],[2.,2.,0.,0.,-2.,-2.,.2,-.2,0.,0.,-.2,.2],30.)]
    results=Dict{String,Any}()
    for (name,m,y0,T) in cases
        f=name=="single" ? ((z,t)->kepler_f(z,t,1.,1.)) : ((z,t)->nbody_f(z,t,m,1.))
        tr,ref=read_reference(name)
        hs=[.1,.05,.025,.0125,.001,.0005]
        errors=Float64[]; sampled=Matrix{Float64}[]
        for h in hs
            _,Y,idx=aligned_solution(f,y0,tr;h=h)
            push!(sampled,Y[:,idx]);push!(errors,maximum(abs.(Y[:,idx]-ref)))
        end
        t,Y,idx=aligned_solution(f,y0,tr;h=.001)
        reverse=reversal_error(f,y0,T;h=.001)
        println("\n",name,": 同时刻所有状态分量最大绝对差（无量纲）")
        for (h,e) in zip(hs,errors);@printf("  h=%-8g 与独立 DOP853 差 = %.3e\n",h,e);end
        @printf("  h=0.001 减半差 = %.3e；时间反演误差 = %.3e\n",maximum(abs.(sampled[end-1]-sampled[end])),reverse)
        if name=="single"
            exact=reduce(hcat,parabolic_state.(tr))
            analytic=maximum(abs.(sampled[end-1]-exact))
            @printf("  对抛物线解析状态的最大误差 = %.3e\n",analytic)
            @assert analytic < 1e-8
            results[name]=(;hs,errors,tr,ref,sampled,analytic)
        else
            d=orbit_diagnostics(t,Y,m)
            @printf("  全步检查 max|ΔE|=%.3e, max|ΔP|=%.3e, max|ΔL|=%.3e, 质心误差=%.3e\n",d.energy_error,d.momentum_error,d.angular_error,d.cm_error)
            @printf("  初始内部能量=%.9e；全步最小两体距离=%.6f\n",d.internal_energy,minimum(d.min_distance))
            perts=NamedTuple[]
            for delta in [1e-6,1e-8]
                yp=copy(y0);yp[name=="two" ? 8 : 11]+=delta
                _,Yp,ip=aligned_solution(f,yp,tr;h=.001)
                _,Ypf,ipf=aligned_solution(f,yp,tr;h=.0005)
                rel=relative_positions(Yp[:,ip],m)-relative_positions(sampled[end-1],m)
                dr=[norm(rel[:,k]) for k in eachindex(tr)]
                signal_error=maximum(abs.((Yp[:,ip]-sampled[end-1])-(Ypf[:,ipf]-sampled[end])))
                independent_error=NaN
                if delta==1e-6
                    tp,rp=read_reference(name*"_pert")
                    @assert maximum(abs.(tp-tr))<1e-12
                    independent_error=maximum(abs.((Yp[:,ip]-sampled[end-1])-(rp-ref)))
                    @printf("  扰动信号与独立 DOP853 的最大差 = %.3e\n",independent_error)
                    @assert independent_error < 1e-8
                end
                @printf("  Δv=%.1e：末态去质心位置差范数=%.6e，扰动信号减半差=%.3e\n",delta,dr[end],signal_error)
                push!(perts,(;delta,dr,signal_error,independent_error))
                @assert signal_error < 1e-9
            end
            @assert d.cm_error < 1e-9
            results[name]=(;hs,errors,tr,ref,sampled,d,t,Y,m,perts)
        end
        @assert errors[end-1] < 1e-8
        @assert reverse < 1e-8
    end
    results
end

function figure8_validation(y0,T;periods=20)
    periods==20 || throw(ArgumentError("独立参考数据覆盖固定 20 周期"))
    tr,R=read_reference("figure8")
    @assert abs(tr[end]-periods*T)<1e-10
    f=(z,t)->nbody_f(z,t,ones(3),1.)
    perperiod=600
    ns=[300,600,1200,2400,4800,9600]
    closure=Float64[];reference_errors=Float64[]
    for n in ns
        t,Y=integrate(f,y0,(0.,T);h=T/n)
        stride=max(1,n÷perperiod);ri=n<perperiod ? (1:2:perperiod+1) : (1:perperiod+1)
        push!(closure,maximum(abs.(Y[:,end]-y0)))
        push!(reference_errors,maximum(abs.(Y[:,1:stride:end]-R[:,ri])))
        @printf("steps/period=%5d, state closure=%.3e, independent difference=%.3e\n",n,closure[end],reference_errors[end])
    end
    n=2400
    t,Y=integrate(f,y0,(0.,periods*T);h=T/n)
    tf,Yf=integrate(f,y0,(0.,periods*T);h=T/(2n))
    idx=1:(n÷perperiod):length(t);idxf=1:(2n÷perperiod):length(tf)
    @assert maximum(abs.(t[idx]-tr))<1e-10
    state_delta=Y[:,idx]-R
    errors=[maximum(abs.(state_delta[:,k])) for k in eachindex(tr)]
    @assert maximum(abs.(Y[:,idx]-Yf[:,idxf])) < 1e-5
    # Convention for this notebook's body ordering: after T/3, (1,2,3) matches (3,1,2).
    perm=[5,6,1,2,3,4,11,12,7,8,9,10]
    shift=n÷3
    choreo=maximum(abs.(Y[:,shift+1:n+1]-Y[perm,1:n+1-shift]))
    d=orbit_diagnostics(t,Y,ones(3))
    period_indices=1:n:length(t)
    position_closure=[maximum(abs.(Y[1:6,k]-y0[1:6])) for k in period_indices]
    velocity_closure=[maximum(abs.(Y[7:12,k]-y0[7:12])) for k in period_indices]
    # Tangent projection of position error onto reference velocity: a local phase-time estimate.
    phase=[dot(state_delta[1:6,k],R[7:12,k])/sum(abs2,R[7:12,k]) for k in eachindex(tr)]
    transverse=[norm(state_delta[1:6,k]-phase[k]*R[7:12,k]) for k in eachindex(tr)]
    @printf("T/3 完整状态置换残差 = %.3e\n",choreo)
    @printf("20 周期 max|ΔE|=%.3e, max|ΔP|=%.3e, max|ΔL|=%.3e, min distance=%.6f\n",d.energy_error,d.momentum_error,d.angular_error,minimum(d.min_distance))
    @printf("20 周期末位置闭合=%.3e, 速度闭合=%.3e；与同初值独立轨迹最大差=%.3e\n",position_closure[end],velocity_closure[end],maximum(errors))
    @printf("同初值独立轨迹自身一周期闭合残差=%.3e（含给定初值和周期舍入影响）\n",maximum(abs.(R[:,601]-y0)))
    @assert choreo < 1e-5
    @assert closure[end] < 1e-6
    @assert maximum(errors) < 1e-5
    @assert minimum(d.min_distance) > .1
    (;ns,closure,reference_errors,tr,R,t,Y,errors,phase,transverse,d,position_closure,velocity_closure,choreo,T)
end
