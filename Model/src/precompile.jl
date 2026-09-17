using PrecompileTools: @compile_workload, @recompile_invalidations

@compile_workload begin
    DEFAULT_SYS = build_neonatal_ecc_sys() |> mtkcompile
    DEFAULT_PROB = ODEProblem(DEFAULT_SYS, [], (0.0, 100.0second))
end
