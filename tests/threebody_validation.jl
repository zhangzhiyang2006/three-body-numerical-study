@testset "Three-body scientific validation" begin
    @test isdefined(Main, :parabolic_state)
    if isdefined(Main, :parabolic_state)
        @test parabolic_state(0.0) ≈ [0.,2.,1.,0.]
        for t in [0.,1.,10.]
            z = parabolic_state(t)
            dt=1e-5
            derivative=(parabolic_state(t+dt)-parabolic_state(t-dt))/(2dt)
            @test derivative ≈ [z[3:4]; -z[1:2]/norm(z[1:2])^3] atol=1e-8
        end
        m=[.3,.03]
        z=[2.,2.,0.,0.,.2,-.2,-.01,.01]
        f=(z,t)->nbody_f(z,t,m,1.)
        t,Y=integrate(f,z,(0.,1.);h=.01)
        D=orbit_diagnostics(t,Y,m)
        @test D.cm_error < 1e-12
        @test D.energy_error < 1e-10
        shifted=copy(Y);shifted[1:2,:].+=10;shifted[3:4,:].+=10
        @test relative_positions(Y,m) ≈ relative_positions(shifted,m)
        @test reversal_error(f,z,1.;h=.01) < 1e-10
        @test maximum(abs.(m[1]*f(z,0.)[5:6]+m[2]*f(z,0.)[7:8])) < 1e-14
    end
end
@testset "RK45 component scaling and termination" begin
    @test_throws ArgumentError dp_rk45((y,t)->y,[1.],(0.,1.);h0=0.)
    @test_throws ArgumentError dp_rk45((y,t)->y,[1.],(1.,0.))
    @test_throws ErrorException dp_rk45((y,t)->fill(NaN,length(y)),[1.],(0.,1.))
    @test_throws ErrorException dp_rk45((y,t)->1000y,[1.],(0.,1.);h0=.1,hmin=.1,atol=1e-30,rtol=1e-30)
    ts,ms,_,_=dp_rk45((y,t)->y,[1e6],(0.,1.))
    tv,mv,_,_=dp_rk45((y,t)->[y[1],0.],[1e6,1e-12],(0.,1.))
    @test length(tv)==length(ts)
    @test abs(mv[1,end]/(1e6*exp(1))-1)<1e-7
    @test mv[2,end]==1e-12
end
@testset "Figure-eight periodicity and choreography" begin
    y0=[-.97000436,.24308753,0.,0.,.97000436,-.24308753,
        .466203685,.432365730,-.932407370,-.864731460,.466203685,.432365730]
    T=6.32591398
    f=(z,t)->nbody_f(z,t,ones(3),1.)
    t,Y=integrate(f,y0,(0.,T);h=T/2400)
    # Position and velocity return; three labels advance cyclically after T/3.
    @test maximum(abs.(Y[1:6,end]-y0[1:6])) < 1e-6
    @test maximum(abs.(Y[7:12,end]-y0[7:12])) < 1e-6
    perm=[5,6,1,2,3,4,11,12,7,8,9,10]
    @test maximum(abs.(Y[:,801:2401]-Y[perm,1:1601])) < 1e-6
    tr,R=read_reference("figure8")
    @test maximum(abs.(Y[:,1:4:end]-R[:,1:601])) < 1e-7
    D=orbit_diagnostics(t,Y,ones(3))
    @test minimum(D.min_distance) > .6
    @test D.energy_error < 1e-8
end
