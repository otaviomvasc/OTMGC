
import sys
from mip import Model,xsum,minimize,CBC,OptimizationStatus,BINARY,CONTINUOUS,ConstrsGenerator,CutPool
import matplotlib.pyplot as plt
from itertools import product
import numpy as np
from time import time


"""
REFERÊNCIAS:
MÉTODOS DE DECOMPOSIÇÃO DE BENDERS DESENVOLVIDOS PELO PROF RICARDO CAMARGO.
email: Camargo.rico@gmail.com
"""


np.random.seed(0)
GRIDSIZE = 1000
EPSILON = 1e-5

class Data():
    def __init__(self,ni,nj):
         self.ni = ni
         self.nj = nj
         self.f = np.random.randint(1000, high=2500, size=nj)
         self.d = np.random.randint(5, high=30, size=ni)
         self.icoord = np.random.randint(1, high=GRIDSIZE, size=(ni,2))
         self.jcoord = np.random.randint(1, high=GRIDSIZE, size=(nj,2))
         c = self.icoord[:,np.newaxis,:] - self.jcoord[np.newaxis,:,:]
         # Euclidean distance matrix
         self.c = np.linalg.norm(c,axis=-1)
         # argument sort of c : sorted_c has indices
         self.sorted_c = np.argsort(self.c,axis=1)

class Benders():
    def __init__(self,dat : Data, is_hotstart = False, is_knapsack = False):
        self.dat = dat
        self.create_master_problem()
        self.is_hotstart = is_hotstart
        if is_hotstart == True:
            self.is_knapsack = True
        else:
            self.is_knapsack = is_knapsack


    def create_master_problem(self):
        dat = self.dat
        m = Model(name='PMR', solver_name=CBC)
        m.verbose = 0
        I,J = range(dat.ni),range(dat.nj)

        y = [m.add_var(var_type=BINARY, name='y(%d)'%j) for j in J]
        eta = m.add_var(var_type=CONTINUOUS,lb=0.0, name='eta')
        m.y = y
        m.eta = eta

        m.objective = minimize(eta + xsum(dat.f[j] * y[j] for j in J))

        m += xsum(y[j] for j in J) >= 1
        #m.write('pmr.lp')
        self.m = m

    def run(self):
        print('\n\n\n')
        dat = self.dat
        I,J = range(dat.ni),range(dat.nj)
        ub,lb = float('inf'),-float('inf')
        h = 0
        start_time = time()
        if self.is_hotstart == True:
            ittype = 'h'
        else:
            ittype = 'i'
        while (ub - lb > EPSILON):
            h += 1
            lb, _eta, _y = self.solve_master_problem()
            assert _y.all() != None, '\n\nerror running the Benders algorithm\n\n'
            if self.is_knapsack == False:
                phi, _u, _v = self.solve_dual_subproblem(_y)
            else:
                phi, _u, _v = self.solve_dual_subproblem_knapsack(_y)
            self.add_benders_cuts(_u, _v)
            sup = lb - _eta + phi
            ub = min(ub, sup)
            run_time = time() - start_time
            self.print_iteration_info(ittype, h, sup, ub, lb, run_time)

            if self.is_hotstart == True and ub - lb < EPSILON:
                self.is_hotstart = False
                ub = float('inf')
                ittype = 'i'

    def print_iteration_info(self, ittype, h, sup, ub, lb, rt):
        gap = 100.0 * (ub - lb) / ub
        print("{:s} ".format(ittype), end='')
        print("{:3d} ".format(h), end='')
        print(" | {:15,.2f}".format(sup), end='')
        print(" | {:15,.2f}".format(ub), end='')
        print(" | {:15,.2f}".format(lb), end='')
        print(" | {:6,.2f} %".format(gap), end='')
        print(" | {:10,.2f} s".format(rt))

    def add_benders_cuts(self, _u, _v):
        m = self.m
        dat = self.dat
        I, J = range(dat.ni), range(dat.nj)
        y, eta = m.y, m.eta
        m += eta >= sum(_u) - xsum(sum(_v[i][j] for i in I) * y[j] for j in J)

    def solve_dual_subproblem(self, _y):

        dat = self.dat
        I, J = range(dat.ni), range(dat.nj)
        c, d = dat.c, dat.d
        O = [j for j in J if _y[j] > 0.9]
        C = [j for j in J if _y[j] < 0.1]
        u = [min([d[i] * c[i][j] for j in O]) for i in I]
        v = [[max(u[i] - d[i] * c[i][j], 0.0) if j in C else 0.0 for j in J] for i in I]
        return sum(u), u, v

    def solve_dual_subproblem_knapsack(self, _y):
        dat = self.dat
        I, J = range(dat.ni), range(dat.nj)
        c, d = dat.c, dat.d
        u = np.zeros(ni, dtype=float)
        v = np.zeros((ni, nj), dtype=float)
        phi = 0.0
        for i in I:
            sy = _y[dat.sorted_c[i]]
            csum = np.cumsum(sy)
            k = np.argmax(csum >= 1)
            _u = dat.d[i] * dat.c[i][dat.sorted_c[i][k]]
            u[i] = _u
            phi += _u
            for p in J:
                if p < k:
                    _j = dat.sorted_c[i][p]
                    _v = max(u[i] - dat.d[i] * dat.c[i][_j], 0.0)
                    v[i][_j] = _v
                    phi -= _v * _y[_j]
                else:
                     break

        return phi, u, v

    def solve_master_problem(self):
        m = self.m
        dat = self.dat
        I, J = range(dat.ni), range(dat.nj)
        status = m.optimize(relax=self.is_hotstart)
        lb, _eta, _y = -float('inf'), -float('inf'), None
        if status == OptimizationStatus.OPTIMAL:
            lb = m.objective_value
            _y = np.array([m.y[j].x for j in J])
            _eta = m.eta.x

        return lb, _eta, _y


class BendersMulticut():
    def __init__(self,dat : Data, is_hotstart = False, is_knapsack = False):
        self.dat = dat
        self.create_master_problem()
        self.is_hotstart = is_hotstart
        if is_hotstart == True:
            self.is_knapsack = True
        else:
            self.is_knapsack = is_knapsack

    def create_master_problem(self):
        dat = self.dat
        m = Model(name='PMR', solver_name=CBC)
        m.verbose = 0
        I, J = range(dat.ni), range(dat.nj)
        y = [m.add_var(var_type=BINARY, name='y(%d)' % j) for j in J]
        eta = [m.add_var(var_type=CONTINUOUS, lb=0.0, name='eta(%d)' % i) for i in I]
        m.y = y
        m.eta = eta
        m.objective = minimize(xsum(eta[i] for i in I) + xsum(dat.f[j] * y[j] for j in J))
        m += xsum(y[j] for j in J) >= 1
        # m.write('pmr.lp
        self.m = m

    def solve_dual_subproblem_knapsack(self,_y):
        dat = self.dat
        I, J = range(dat.ni), range(dat.nj)
        c, d = dat.c, dat.d
        u = np.zeros(ni, dtype=float)
        v = np.zeros((ni, nj), dtype=float)
        phi = 0.0
        for i in I:
            sy = _y[dat.sorted_c[i]]
            csum = np.cumsum(sy)
            k = np.argmax(csum >= 1)
            # print(k,sy,csum)
            _u = dat.d[i] * dat.c[i][dat.sorted_c[i][k]]
            u[i] = _u
            phi += _u
            for p in J:
                if p < k:
                    _j = dat.sorted_c[i][p]
                    _v = max(0.0, u[i] - dat.d[i] * dat.c[i][_j])
                    v[i][_j] = _v
                    phi -= _v * _y[_j]
                else:
                    break
        return phi, u, v
    def run(self):
        print('\n\n\n')
        dat = self.dat
        I, J = range(dat.ni), range(dat.nj)
        ub, lb = float('inf'), -float('inf')
        h = 0
        start_time = time()
        if self.is_hotstart == True:
            ittype = 'h'
        else:
            ittype = 'i'
        while (ub - lb > EPSILON):
            h += 1
            if h > 100:
                break
            lb, _eta, _y = self.solve_master_problem()

            assert _y.all() != None, '\n\nerror running the Benders algorithm\n\n'

            if self.is_knapsack == False:
                phi, _u, _v = self.solve_dual_subproblem(_y)
            else:
                phi, _u, _v = self.solve_dual_subproblem_knapsack(_y)

            self.add_benders_cuts(_u, _v)
            sup = lb - _eta.sum() + phi
            ub = min(ub, sup)
            run_time = time() - start_time
            self.print_iteration_info(ittype, h, sup, ub, lb, run_time)

            if self.is_hotstart == True and ub - lb < EPSILON:
                self.is_hotstart = False
                ub = float('inf')
                ittype = 'i'

    def print_iteration_info(self, ittype, h, sup, ub, lb, rt):
        gap = 100.0 * (ub - lb) / ub
        print("{:s} ".format(ittype), end='')
        print("{:3d} ".format(h), end='')
        print(" | {:15,.2f}".format(sup), end='')
        print(" | {:15,.2f}".format(ub), end='')
        print(" | {:15,.2f}".format(lb), end='')
        print(" | {:6,.2f} %".format(gap), end='')
        print(" | {:10,.2f} s".format(rt))

    def add_benders_cuts(self, _u, _v):
        m = self.m
        dat = self.dat
        I, J = range(dat.ni), range(dat.nj)
        y, eta = m.y, m.eta
        for i in I:
            m += eta[i] >= _u[i] - xsum(_v[i][j] * y[j] for j in J)

    def solve_dual_subproblem(self, _y):
        dat = self.dat
        I, J = range(dat.ni), range(dat.nj)
        c, d = dat.c, dat.d
        O = [j for j in J if _y[j] > 0.9]
        C = [j for j in J if _y[j] < 0.1]
        u = [min([d[i] * c[i][j] for j in O]) for i in I]
        v = [[max(u[i] - d[i] * c[i][j], 0.0) if j in C else 0.0 for j in J] for i in I]
        return sum(u), u, v

    def solve_master_problem(self):
        m = self.m
        dat = self.dat
        I,J = range(dat.ni),range(dat.nj)
        status = m.optimize(relax=self.is_hotstart)
        lb,_eta,_y= -float('inf'),None,None
        if status == OptimizationStatus.OPTIMAL:
            lb = m.objective_value
            _y = np.array([m.y[j].x for j in J])
            _eta = np.array([m.eta[i].x for i in I])
        return lb,_eta,_y


class BendersCutCutGenerator(ConstrsGenerator):
    def __init__(self, dat: Data,y,eta):
        self.dat = dat
        self.y,self.eta=y,eta

    def generate_constrs(self,model: Model, depth: int=0, npass: int=0):
        dat = self.dat
        I,J = range(dat.ni),range(dat.nj)

        y,eta = self.y,self.eta
        ly,leta = model.translate(y),model.translate(eta)

        _y = np.array([ly[j].x if ly[j] else 0 for j in J])
        _eta = np.array([leta[i].x if leta[i] else 0 for i in I])

        phi,u,v = self.solve_dual_subproblem_knapsack(_y)
        for i in I:
            if phi[i] - _eta[i] > EPSILON:
                cut = leta[i] >= u[i] - xsum(v[i][j] * ly[j] for j in J)
                model += cut


    def solve_dual_subproblem_knapsack(self,_y):
        dat = self.dat
        I,J = range(dat.ni),range(dat.nj)
        c, d = dat.c, dat.d
        u = np.zeros(ni, dtype=float)
        v = np.zeros((ni, nj), dtype=float)
        phi = np.zeros(ni, dtype=float)
        for i in I:
            sy = _y[dat.sorted_c[i]]
            csum = np.cumsum(sy)
            k = np.argmax(csum >= 1)
            _u = dat.d[i] * dat.c[i][dat.sorted_c[i][k]]
            u[i] = _u
            phi[i] = _u
            for p in J:
                if p < k:
                    _j = dat.sorted_c[i][p]
                    _v = max(0.0, u[i] - dat.d[i] * dat.c[i][_j])
                    v[i][_j] = _v
                    phi[i] -= _v * _y[_j]
                else:
                    break

        return phi, u, v


class BendersBranchAndCut():
    def __init__(self,dat : Data):
        self.dat = dat
        self.create_master_problem()
        print('here')

    def create_master_problem(self):
        dat = self.dat
        m = Model(name='PMR', solver_name=CBC)
        #m.verbose = 0
        I,J = range(dat.ni),range(dat.nj)

        y = [m.add_var(var_type=BINARY, name='y(%d)'%j) for j in J]
        eta = [m.add_var(var_type=CONTINUOUS, lb=0.0,name='eta(%d)'%i) for i in I]
        m.y = y
        m.eta = eta
        m.objective = minimize(xsum(eta[i] for i in I) + xsum(dat.f[j] * y[j] for j in J))
        m += xsum(y[j] for j in J) >= 1
        #m.write('pmr.lp')
        self.m = m

    def run(self):
        m = self.m
        m.cuts_generator = BendersCutCutGenerator(self.dat, m.y, m.eta)
        m.lazy_constrs_generator = BendersCutCutGenerator(self.dat, m.y, m.eta)
        print('\n\n\n')
        m.verbose = 0
        start_time = time()
        status = m.optimize()
        self.run_time = time() - start_time
        print(status)
        if status == OptimizationStatus.OPTIMAL:
            self.is_solution = True
            self.print_solution()

    def print_solution(self):
        if self.is_solution == True:
            dat = self.dat
            m = self.m
            I, J = range(dat.ni), range(dat.nj)
            print("Custo total : {:12,.2f}.".format(m.objective_value))
            print("Tempo total : {:12,.2f}.".format(self.run_time))
            print("facilidades ")
            for j in J:
                if m.y[j].x > 1e-6:
                   print("{:5d} ".format(j + 1), end='')
            print()
        else:
            print()
            print('Nao ha solucao disponivel para impressao')
            print()


class UFLP():
    def __init__(self,dat : Data):
        self.is_solution = False
        self.dat = dat
        self.build_model()

    def run(self):
        print('\n\n\n')
        m = self.m
        m.verbose = 0
        start_time = time()
        status = m.optimize()
        self.run_time = time() - start_time
        if status == OptimizationStatus.OPTIMAL:
            self.is_solution = True

    def print_solution(self):
        if self.is_solution == True:
            dat = self.dat
            m = self.m
            I, J = range(dat.ni), range(dat.nj)
            print("Custo total de instalacao: {:12,.2f}".format(sum([m.y[j].x * dat.f[j] for j in J])))
            print("Custo total de transporte: {:12,.2f} ".format(sum([m.x[i, j].x * dat.d[i] * dat.c[i][j] for (i, j) in product(I, J)])))
            print("Custo total : {:12,.2f}.".format(m.objective_value))
            print("Tempo total : {:12,.2f}.".format(self.run_time))
            print("facilidades : demanda : clientes ")
            for j in J:
                if m.y[j].x > 1e-6:
                    print("{:11d} : {:7.0f} : ".format(j + 1, sum([m.x[i, j].x * dat.d[i] for i in I])), end='')
            for i in I:
                if m.x[i, j].x > 1e-6:
                    print(" {:d}".format(i + 1), end='')
            print()
        else:
            print()
            print('Nao ha solucao disponivel para impressao')
            print()

    def build_model(self):
        dat = self.dat
        m = Model(name='uflp', solver_name=CBC)
        I = range(dat.ni)
        J = range(dat.nj)
        y = [m.add_var(var_type=BINARY, name='y(%d)' % j) for j in J]
        x = {(i, j): m.add_var(lb=0.0, name='x(%d,%d)' % (i, j)) for j in J for i in I}
        m.x = x
        m.y = y
        m.objective = minimize(xsum(dat.f[j] * y[j] for j in J) + xsum(dat.d[i] * dat.c[i][j] * x[i, j] for
                                                                       (i, j) in product(I, J)))
        for i in I:
            m += xsum(x[i, j] for j in J) == 1
        for (i, j) in product(I, J):
            m += x[i, j] <= y[j]

        self.m = m


if __name__ == "__main__":
    args = sys.argv
    #assert len(args) > 2, f"\n\nuse: %s ni nj\n\n" % args[0]
    ni, nj = 10, 10
    dat = Data(ni, nj)
    uflp = UFLP(dat)
    uflp.run()
    uflp.print_solution()

    # bd = Benders(dat)
    # bd.run()

    # bd = Benders(dat,is_hotstart=False,is_knapsack=True)
    # bd.run()

    # bd = Benders(dat,is_hotstart=True)
    # bd.run()

    # bdmcut = BendersMulticut(dat)
    # bdmcut.run()

    # bdmcut = BendersMulticut(dat, is_hotstart=True)
    # bdmcut.run()
    bdbc = BendersBranchAndCut(dat)
    bdbc.run()


