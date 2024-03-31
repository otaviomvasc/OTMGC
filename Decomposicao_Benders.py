import sys
from mip import Model,xsum,minimize,CBC,OptimizationStatus,BINARY,CONTINUOUS,ConstrsGenerator,CutPool
import matplotlib.pyplot as plt
from itertools import product
import numpy as np
from time import time

"""
Classe para resolver o problema de alocação através
da decomposição de benders
"""
EPSILON = 1e-5
ni, nj = 30,10
class Benders():
    def __init__(self, dat, is_hotstart = True, is_knapsack = False):
        self.dat = dat
        self.create_master_problem()
        self.is_hotstart = is_hotstart
        if is_hotstart == True:
            self.is_knapsack = True
        else:
            self.is_knapsack = is_knapsack

    def create_master_problem(self):
        """
        Dado que a variável inteira "Dificil" não muda nas formulações, o PMS também se mantem apenas com:
         variável Y: Qual localização instalar o facilite
        variável eta: maior menor valor da FO dada determinada alocação Yj
        """

        dat = self.dat
        m = Model(name='PMR', solver_name=CBC)
        m.verbose = 0
        I,J = range(dat['clientes']),range(dat['plantas']) #TODO: Não é melhor criar a classe DAT?

        y = [m.add_var(var_type=BINARY, name='y(%d)'%j) for j in J]
        eta = [m.add_var(var_type=CONTINUOUS, lb=0.0, name='eta(%d)' % i) for i in I]
        m.y = y
        m.eta = eta

        m.objective = minimize(sum(7500 * y[j] for j in J) + xsum(eta[i] for i in I))

        m += xsum(y[j] for j in J) >= 1  #Corte de viabilidade, já que solução viável precisa ter pelo menos 1 facility aberto
        #m.write('pmr.lp')
        self.m = m

    def testa_calculo_dual(self, y):
        a = [7 ,6, 0, 24, 7, 0, 1, 5, 7, 7, 3, 10, 5, 0, 6, 7, 3, 8, 3, 6, 3, 6, 10, 0, 19, 10, 12, 10, 19, 0, 0, 10, 0, 16, 19, 19, 17, 5, 23, 24, 19 ,3, 7, 6, 22, 23, 23, 6, 24, 19]
        melhor_solucao = np.unique(a)
        _y_teste = [0 if p not in melhor_solucao else 1 for p in range(len(y))]
        phi, _k, _v = self.solve_dual_problem_Segunda_Formulacao(_y_teste)
        custo_total  = len(melhor_solucao) * 7500 + phi
        solucao_otima = 796648
        b=0

    def run(self):
        print('\n\n\n')
        dat = self.dat  #porque essa atribuição?
        I,J = range(dat['clientes']),range(dat['plantas'])  #porque essa atribuição?
        ub,lb = float('inf'),-float('inf')
        h = 0
        start_time = time()
        if self.is_hotstart == True:
            ittype = 'h'
        else:
            ittype = 'i'
        while (ub - lb > EPSILON):
            if h > 100:
                self.is_hotstart = False
                lb, _eta, _y = self.solve_master_problem()

            h += 1
            lb, _eta, _y = self.solve_master_problem()
            assert _y.all() != None, '\n\nerror running the Benders algorithm\n\n'
            #self.testa_calculo_dual(_y)
            phi, _k, _v = self.solve_dual_problem_Segunda_Formulacao(_y)

            #self.add_benders_cuts(_u, _v)
            self.add_benders_cuts_formulacao(_k, _v)
            sup = lb - sum(_eta) + phi
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

    def add_benders_cuts_formulacao(self, _k, _v):
        m = self.m
        y, eta = m.y, m.eta

        for c in range(self.dat['clientes']):
            if _k[c] == 0:
                m += eta[c] >= self.dat['matriz_custo'][(self.dat['D_ord'][c][0], c)]
            else:
                if _k[c] == (len(y) - 1):
                    _k[c] = 23
                lista_ord = self.dat['D_ord'][c]
                par_s = (lista_ord[_k[c] + 1], c)
                par_a = (lista_ord[_k[c]], c)
                m += (eta[c] + xsum((self.dat['matriz_custo'][par_s] - self.dat['matriz_custo'][(j, c)])
                               * y[j] for j in lista_ord[:_k[c]]) >= self.dat['matriz_custo'][par_s])


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

    def solve_dual_problem_Segunda_Formulacao(self, _y):
        """
        Método que vai calcular de forma analitica o do problema Dual de acordo com o y encontrado no PMS!
        :param _y:
        :return: valor da FO e variáveis duais para criação dos cortes!
        """

        #Cálculo do meu K, sendo que k/ é o índice!!!
        k = dict()
        for c in range(self.dat['clientes']):
            aux = 0
            result = list()
            it_aux = [_y[p] for p in self.dat['D_ord'][c]]
            if _y[self.dat['D_ord'][c][0]] >= 1:
            #Estou implementando da forma mais parecida com o artigo para melhor entendimento.
            #Após a conferência do resultado, voltar e otimizar o código.
            #Artigo também indica que se solução for binária é mais fácil calcular, porém vou deixar assim para tentar ficar genérico
                k[c] = 0
            else:
                for i in range(len(self.dat['D_ord'][0])):
                    if it_aux[i] >= aux and it_aux[i] < 1:
                        aux = it_aux[i]
                        result.append(i)
                k[c] = max(result)

        #Após o cálculo do k~ tem que se calcular a F.O otima para dado k e y
        #Mesma ideia da parte de cima. Fazer código explicito e depois refatorar!
        v = dict()
        for c in range(self.dat['clientes']):
            v[c] = list()
            for pos in range(self.dat['clientes'] - 1):
                if k[c] == self.dat['plantas'] - 1 and pos == k[c]:
                    k[c] = k[c] - 1 #A planta aberta é a mais longe do cliente e não consigo pegar a distância k+1
                if pos <= k[c]:
                    ik1 = self.dat['D_ord'][c][pos + 1]
                    ik = self.dat['D_ord'][c][pos]
                    Dk2 = self.dat['matriz_custo'][(ik1, c)]
                    Dk1 = self.dat['matriz_custo'][(ik, c)]
                    v[c].append(Dk2 - Dk1)
                else:
                    v[c].append(0)

        return round(sum(sum(v[c]) for c in range(self.dat['clientes']))), k, v

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
        I,J = range(dat['clientes']),range(dat['plantas'])
        status = m.optimize(relax=self.is_hotstart)
        lb, _eta, _y = -float('inf'), -float('inf'), None
        if status == OptimizationStatus.OPTIMAL:
            lb = m.objective_value
            _y = np.array([m.y[j].x for j in J])
            _eta = np.array([m.eta[i].x for i in I])


        return lb, _eta, _y