import maxwell_solver as ms


def test_read_resonant_box_1D():
    adapter = ms.OpensembaAdapter('resonant_box_1D.smb.json')

    opts = adapter.readSolverOptions()
    pd = adapter.readProblem()

    solver = ms.Solver(pd, opts)

    Dx = solver.evolution.getDerivativeOperator('E', 'x')

    assert True
