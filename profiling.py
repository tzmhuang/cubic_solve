import numpy as np
import cProfile
from cubic_polynomial import CubicPolynomial 





np.random.seed(1999)

class SolverProfiler():
    def __init__(self, N): 
        random_coeff = [np.random.rand(4)*20-10 for i in range(N)]  #generate random coeff in range [-10,10)
        self.poly_list = [CubicPolynomial(c[0], c[1], c[2], c[3]) for c in random_coeff]
    
    def EigenProfiler(self):
        for p in self.poly_list:
            p.Solver.EigenSolve()
        return 

    def SJProfiler(self):
        for p in self.poly_list:
            p.Solver.SJSolve()
        return

def Spin():
    print('-- Profiling CubicPolynomial.Solver.EigenSolve()\n')
    cProfile.run('profile.EigenProfiler()')
    print('-- Profiling CubicPolynomial.Solver.SJSolve()\n')
    cProfile.run('profile.SJProfiler()')

if __name__ == "__main__":
    msg = "--- Profiler for CubicPolynomial Solvers ---"
    print(msg)
    N = input('Set N = ')
    N = int(N)
    profile = SolverProfiler(N)
    Spin()

    



