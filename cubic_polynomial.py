import numpy as np
import matplotlib.pyplot as plt


class Solver():
    def __init__(self, polynomial):
        self.polynomial = polynomial
        self.a3 = polynomial.a3
        self.a2 = polynomial.a2
        self.a1 = polynomial.a1
        self.a0 = polynomial.a0

    def EigenSolve(self):
        c0 = self.polynomial.a0/self.polynomial.a3
        c1 = self.polynomial.a1/self.polynomial.a3
        c2 = self.polynomial.a2/self.polynomial.a3
        c3 = 1
        A = np.array([[0,0,-1*c0],[1,0,-1*c1],[0,1,-1*c2]])
        eig_v = np.linalg.eig(A)
        return eig_v[0]  #eigen value of A are the roots for f(x)

    def SJSolve(self):
        roots = [0,0,0]
        A = self.a2**2 - 3 * self.a3 * self.a1
        B = self.a2 * self.a1 - 9 * self.a3 * self.a0
        C = self.a1**2 - 3 * self.a2 * self.a0
        Delta = B**2 - 4 * A * C
        #print('Delta = : {}'.format(Delta))

        if A == 0 and B == 0:   # case 1: A == B == 0
            #print('<Case 1>')
            if self.a3 != 0:
                roots[0] =  -1/3 * self.a2 * (1/self.a3)
            elif self.a2 != 0:
                roots[0] = -1 * self.a1 * (1/self.a2)
            elif self.a1 != 0:
                roots[0] = -3 * self.a0 * (1/self.a1)
            elif self.a0 != 0:
                raise ValueError('f(x) = 0, f(x) = a0, a0 != 0')
            roots[1] = roots[0]
            roots[2] = roots[0]
            return roots

        if Delta > 0:         # case 2: Delta > 0  
            #print('<Case 2>')
            Y1 = A*self.a2 + 3*self.a3*(-1*B + np.sqrt(Delta))/2
            Y2 = A*self.a2 + 3*self.a3*(-1*B - np.sqrt(Delta))/2
            if Y1 > 0:
                cbrtY1 = np.cbrt(Y1)
            else:
                cbrtY1 = -1 * np.cbrt(-1 * Y1)
            if Y2 > 0:
                cbrtY2 = np.cbrt(Y2)
            else:
                cbrtY2 = -1 * np.cbrt(-1 * Y2)
            roots[0] = (-1*self.a2 - (cbrtY1 + cbrtY2))/(3*self.a3)
            roots[1] = (-1*self.a2 + 0.5*(cbrtY1 + cbrtY2) \
                    + np.sqrt(3)*0.5*(cbrtY1 - cbrtY2)*1j)*(1/(3*self.a3))
            roots[2] = (-1*self.a2 + 0.5*(cbrtY1 + cbrtY2) \
                    - np.sqrt(3)*0.5*(cbrtY1 - cbrtY2)*1j)*(1/(3*self.a3))
            return roots


        if Delta == 0:
            #print('<Case 3>')
            roots[0] = -1*self.a2/self.a3 + B/A
            roots[1] = -1*B/(2*A)
            roots[2] = roots[1]
            return roots

        if Delta < 0:
            #print('<Case 4>')
            T = (2*A*self.a2 - 3*self.a3*B)*(1/(2*np.sqrt(A**3)))
            theta = np.arccos(T)
            roots[0] = (-1*self.a2 - 2*np.sqrt(A)*np.cos(theta*(1/3)))*(1/(3*self.a3))
            roots[1] = (-1*self.a2 + np.sqrt(A)*(np.cos(theta*(1/3)) \
                    + np.sqrt(3)*np.sin(theta*(1/3))))*(1/(3*self.a3))
            roots[2] = (-1*self.a2 + np.sqrt(A)*(np.cos(theta*(1/3)) \
                    - np.sqrt(3)*np.sin(theta*(1/3))))*(1/(3*self.a3))
            return roots



            


class CubicPolynomial():
    def __init__(self, a3, a2, a1, a0):
        self.a3 = a3
        self.a2 = a2
        self.a1 = a1
        self.a0 = a0
        self.Solver = Solver(self)
    
    def Evaluate(self, x, order):
        if order == 0:
            return ((((self.a3) * x+self.a2)*x+self.a1)*x+self.a0)
        elif order == 1:
            return (((3*self.a3*x) * x+2*self.a2)* x+self.a1)
        elif order == 3:
            return ((6*self.a3*x) * x+4*self.a2)
        else:
            print("Evaluation order out of range!")


if __name__ == "__main__":
    msg = "--- Cubic Root Solver ---"
    print(msg)
    inp = input("Input Coeff: ")
    coeff = inp.split(" ")
    coeff = [float(i) for i in coeff]
    poly = CubicPolynomial(coeff[0],coeff[1],coeff[2],coeff[3])
    print('Évaluate at -2: ', poly.Evaluate(-2,0))
    roots = poly.Solver.EigenSolve()
    sj_roots = poly.Solver.SJSolve()
    print("Roots of the polynomial are: \n")
    print("<EigenMethod>: [1] {0}, [2] {1}, [3] {2}\n".format(roots[0], roots[1], roots[2]))
    print("  <SJ_Method>: [1] {0}, [2] {1}, [3] {2}\n".format(sj_roots[0], sj_roots[1], sj_roots[2]))

    # Plotting 
    # x = np.linspace(-3,23,100)
    # y = [poly.Evaluate(i,0) for i in x]
    # fig,ax = plt.subplots()
    # ax.plot(x,y)
    # ax.plot(x,[0 for i in x], 'g-') # zero line
    # for r in roots:
    #     ax.plot(r,0,'rx')
    # for r in sj_roots:
    #     ax.plot(r,0,'gx')
    # plt.show()
    









