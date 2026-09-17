using PrecompileTools: @compile_workload, @recompile_invalidations, @setup_workload

@setup_workload begin
    @parameters Ca = 0μM ROS = 0μM
    @parameters ATP = 5000μM ISO = 0μM
    @compile_workload begin
        b1arsys = get_bar_sys(ATP, ISO) |> mtkcompile
        camsys = get_camkii_sys(; Ca=Ca, ROS=ROS) |> mtkcompile
        DEFAULT_SYS = build_neonatal_ecc_sys() |> mtkcompile
        DEFAULT_PROB = ODEProblem(DEFAULT_SYS, [], (0.0, 100.0second))
    end
end
