import maxwell_solver as ms


def test_read_resonant_box_1D():
    adapter = ms.adapter('')

    opts = adapter.readSolverOptions()
    pd = adapter.readProblemDescription()

    solver = ms.Solver(pd, opts)

    Dx = solver.evolution.getDerivativeOperator('E', 'x')

    assert True
