@compile_workload begin
    sys = build_neonatal_ecc_sys() |> mtkcompile
    prob = ODEProblem(sys, [], (0.0, 200.0second))
    sol = solve(prob, TRBDF2())
end
