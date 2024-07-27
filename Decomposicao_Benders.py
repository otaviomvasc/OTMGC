import sys

from mip import Model,xsum,minimize,CBC,OptimizationStatus,BINARY,CONTINUOUS,ConstrsGenerator,CutPool, MAXIMIZE, MINIMIZE
import matplotlib.pyplot as plt
from itertools import product
import numpy as np
from time import time
from copy import deepcopy

"""
Classe para resolver o problema de alocação através
da decomposição de benders. Métodos copiados do arquivo Decomposicao_Bender_Prof
"""
EPSILON = 1
ni, nj = 30,10
class Benders():
    def __init__(self, dat, is_hotstart = True, is_knapsack = False, multicut = True, inspecao=True):

        """
        :param dat: Dados
        :param is_hotstart:
        :param is_knapsack:
        :param multicut: Se True, será resolvido por multicut e se false por single cut
        :param inpecao: Se True. será resolvido por inspecao e se false por método exato

        """
        self.dat = dat
        self.multicut = multicut
        self.inspecao = inspecao
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
        if self.multicut:
            eta = [m.add_var(var_type=CONTINUOUS, lb=0.0, name='eta(%d)' % i) for i in I]
            m.objective = minimize(sum(dat['custo_abertura'][j] * y[j] for j in J) + xsum(eta[i] for i in I))
        else:
            eta = m.add_var(var_type=CONTINUOUS, lb=0.0, name='eta')
            m.objective = minimize(sum(dat['custo_abertura'][j] * y[j] for j in J) + eta)

        m.y = y
        m.eta = eta


        m += xsum(y[j] for j in J) >= 1  #Corte de viabilidade, já que solução viável precisa ter pelo menos 1 facility aberto
        #m.write('pmr.lp')
        self.m = m

    def testa_calculo_dual(self, y):
        a = [7 ,6, 0, 24, 7, 0, 1, 5, 7, 7, 3, 10, 5, 0, 6, 7, 3, 8, 3, 6, 3, 6, 10, 0, 19, 10, 12, 10, 19, 0, 0, 10, 0, 16, 19, 19, 17, 5, 23, 24, 19 ,3, 7, 6, 22, 23, 23, 6, 24, 19]
        melhor_solucao = np.unique(a)
        _y_teste = [0 if p not in melhor_solucao else 1 for p in range(len(y))]

        phi, _k, _v = self.solve_dual_problem_Segunda_Formulacao(_y_teste)
        custo_total = (len(melhor_solucao)-1) * 7500 + phi
        solucao_otima = 796648
        b=0

    def solve_dual_problem_solver(self,_y):
        phi = list()
        v_final = list()
        for i in range(self.dat['clientes']): #Para cada cliente
            m1 = Model(name='DualProblem', sense=MAXIMIZE)
            var_v = {v: m1.add_var(name=f'v{v}',var_type=CONTINUOUS, lb=0.0) for v in range(self.dat['plantas'])}
            m1.verbose = 0
            m1.objective = self.dat['matriz_custo'][self.dat['D_ord'][i][0], i] \
                            + ( var_v[0] * (1 - sum(_y[j] for j in range(self.dat['plantas']) 
                                                  if self.dat['matriz_custo'][self.dat['D_ord'][i][0], i] == self.dat['matriz_custo'][j, i])))\
                            - xsum(sum(_y[j] for j in range(self.dat['plantas']) if self.dat['matriz_custo'][j,i] == self.dat['matriz_custo'][self.dat['D_ord'][i][k], i])
                                   * var_v[k] for k in range(len(self.dat['D_ord'][i])) if k>0)

            for k in range(len(self.dat['D_ord'][i])):
                if k == len(self.dat['D_ord'][i]) - 1:
                    continue
                #m1 += var_v[self.dat['D_ord'][i][k]] - var_v[self.dat['D_ord'][i][k + 1]] <= self.dat['matriz_custo'][i,self.dat['D_ord'][i][k+1]] - self.dat['matriz_custo'][i,self.dat['D_ord'][i][k]]
                m1 += var_v[k] - var_v[k+1] <= self.dat['matriz_custo'][self.dat['D_ord'][i][k+1],i] - self.dat['matriz_custo'][self.dat['D_ord'][i][k],i]
            status = m1.optimize()
            #if status == OptimizationStatus.OPTIMAL:
            phi.append(deepcopy(m1.objective_value))
            v_final.append([var_v[v_t].x for v_t in var_v ])
        return sum(phi), v_final

    def add_cortes_exato(self, v, _y):
        _y = self.m.y
        if self.multicut:
            for i in range(self.dat['clientes']):
                self.m += self.m.eta[i] >= self.dat['matriz_custo'][self.dat['D_ord'][i][0], i] \
                            + ( v[i][0] * (1 - sum(_y[j] for j in range(self.dat['plantas']) 
                                                  if self.dat['matriz_custo'][self.dat['D_ord'][i][0], i] == self.dat['matriz_custo'][j, i])))\
                            - xsum(sum(_y[j] for j in range(self.dat['plantas']) if self.dat['matriz_custo'][j,i] == self.dat['matriz_custo'][self.dat['D_ord'][i][k], i])
                                   * v[i][k] for k in range(len(self.dat['D_ord'][i])) if k>0)

        else:
            a=0
            self.m += self.m.eta >=  xsum([self.dat['matriz_custo'][self.dat['D_ord'][i][0], i] \
                            + (v[i][0] * (1 - sum(_y[j] for j in range(self.dat['plantas']) 
                                                  if self.dat['matriz_custo'][self.dat['D_ord'][i][0], i] == self.dat['matriz_custo'][j, i])))\
                            - xsum(sum(_y[j] for j in range(self.dat['plantas']) if self.dat['matriz_custo'][j,i] == self.dat['matriz_custo'][self.dat['D_ord'][i][k], i])
                                   * v[i][k] for k in range(len(self.dat['D_ord'][i])) if k>0) for i in range(self.dat['clientes'])])




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
            h += 1
            lb, _eta, _y = self.solve_master_problem()
            assert _y.all() != None, '\n\nerror running the Benders algorithm\n\n'
            if self.inspecao:
                phi, _k, _v = self.solve_dual_problem_Segunda_Formulacao(_y)
                self.add_benders_cuts_formulacao(_k, _v)
            else:
                phi, _v = self.solve_dual_problem_solver(_y)
                self.add_cortes_exato(_v, _y)
            
            if self.multicut:
                sup = lb - sum(_eta) + phi
            else:
                sup = lb - _eta + phi
            ub = min(ub, sup)
            run_time = time() - start_time
            self.print_iteration_info(ittype, h, sup, ub, lb, run_time)

            if self.is_hotstart == True and ub - lb < EPSILON:
                self.is_hotstart = False
                ub = float('inf')
                ittype = 'i'

        return _y, ub

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

        if self.multicut:
            for c in range(self.dat['clientes']):
                if _k[c] == 0:
                    m += eta[c] >= self.dat['matriz_custo'][(self.dat['D_ord'][c][0], c)]
                else:
                    lista_ord = self.dat['D_ord'][c]
                    par_s = (lista_ord[_k[c]], c)
                    m += (eta[c] >= self.dat['matriz_custo'][par_s] - xsum((self.dat['matriz_custo'][par_s] - self.dat['matriz_custo'][(j, c)])
                                   * y[j] for j in lista_ord[:_k[c]]))

        else:
            #single cut
            t = list()
            for c in range(self.dat['clientes']):
                if _k[c] == 0:
                    t.append(self.dat['matriz_custo'][(self.dat['D_ord'][c][0], c)])
                else:
                    lista_ord = self.dat['D_ord'][c]
                    par_s = (lista_ord[_k[c]], c)
                    t.append(
                        self.dat['matriz_custo'][par_s] - xsum(
                            (self.dat['matriz_custo'][par_s] - self.dat['matriz_custo'][(j, c)])
                            * y[j] for j in lista_ord[:_k[c]])
                    )

            m += eta >= sum(t)


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
            it_aux = [_y[p] for p in self.dat['D_ord'][c]]
            if _y[self.dat['D_ord'][c][0]] >= 1:
                k[c] = 0
            else:
                for i in range(len(self.dat['D_ord'][0])):
                    aux += it_aux[i]
                    if aux >= 1:
                        k[c] = i
                        break

        #Após o cálculo do k~ tem que se calcular a F.O otima para dado k e y
        #Mesma ideia da parte de cima. Fazer código explicito e depois refatorar!

        v = dict()
        for c in range(self.dat['clientes']):
            v[c] = list()
            if k[c] == 0:
                ik = self.dat['D_ord'][c][0]
                Dk2 = self.dat['matriz_custo'][(ik, c)]
                v[c].append(Dk2)
            else:
                ik1 = self.dat['D_ord'][c][k[c]]
                Dk2 = self.dat['matriz_custo'][(ik1, c)]
                aux = 0
                for ll in range(len(self.dat['D_ord'][c])):
                    if ll > k[c]:
                        break
                    planta = self.dat['D_ord'][c][ll]
                    aux += (Dk2 - self.dat['matriz_custo'][(planta, c)]) *_y[planta]
                vf = Dk2 - aux
                v[c].append(vf)

        return round(sum(sum(v[c]) for c in range(self.dat['clientes']))), k, v


    def solve_master_problem(self):
        m = self.m
        dat = self.dat
        I,J = range(dat['clientes']),range(dat['plantas'])
        status = m.optimize(relax=self.is_hotstart)
        lb, _eta, _y = -float('inf'), -float('inf'), None
        if status == OptimizationStatus.OPTIMAL:
            lb = m.objective_value
            _y = np.array([m.y[j].x for j in J])
            if self.multicut:
                _eta = np.array([m.eta[i].x for i in I])
            else:
                _eta = m.eta.x


        return lb, _eta, _y


class BendersCutCutGenerator(ConstrsGenerator):
    def __init__(self, dat,y,eta):
        self.dat = dat
        self.y,self.eta=y,eta

    def generate_constrs(self,model: Model, depth: int=0, npass: int=0):
        dat = self.dat
        I,J = range(dat['clientes']),range(dat['plantas'])

        y,eta = self.y,self.eta
        ly,leta = model.translate(y),model.translate(eta)

        _y = np.array([ly[j].x if ly[j] else 0 for j in J])
        _eta = np.array([leta[i].x if leta[i] else 0 for i in I])


        phi,_k,_v = self.solve_dual_problem(_y)
        self.add_benders_cuts_formulacao(_k, _v, model, ly, leta)
        # for i in I:
        #     if phi[i] - _eta[i] > EPSILON:
        #         cut = leta[i] >= u[i] - xsum(v[i][j] * ly[j] for j in J)
        #         model += cut

    def add_benders_cuts_formulacao(self, _k, _v, model, ly, leta):
        m = model
        y, eta = ly, leta

        for c in range(self.dat['clientes']):
            if _k[c] == 0:
                m += eta[c] >= self.dat['matriz_custo'][(self.dat['D_ord'][c][0], c)]
            else:
                lista_ord = self.dat['D_ord'][c]
                par_s = (lista_ord[_k[c]], c)
                m += (eta[c] >= self.dat['matriz_custo'][par_s] - xsum((self.dat['matriz_custo'][par_s] - self.dat['matriz_custo'][(j, c)])
                               * y[j] for j in lista_ord[:_k[c]]))
    def solve_dual_problem(self, _y):
        k = dict()
        for c in range(self.dat['clientes']):
            aux = 0
            result = list()
            it_aux = [_y[p] for p in self.dat['D_ord'][c]]
            if _y[self.dat['D_ord'][c][0]] >= 1:
                # Estou implementando da forma mais parecida com o artigo para melhor entendimento.
                # Após a conferência do resultado, voltar e otimizar o código.
                # Artigo também indica que se solução for binária é mais fácil calcular, porém vou deixar assim para tentar ficar genérico
                k[c] = 0
            else:
                for i in range(len(self.dat['D_ord'][0])):
                    # if it_aux[i] >= aux and it_aux[i] <= 1:
                    # aux = it_aux[i]
                    # result.append(i)
                    aux += it_aux[i]
                    if aux >= 1:
                        k[c] = i
                        break

        # Após o cálculo do k~ tem que se calcular a F.O otima para dado k e y
        # Mesma ideia da parte de cima. Fazer código explicito e depois refatorar!
        dual = False
        if dual:
            v = dict()
            for c in range(self.dat['clientes']):
                v[c] = list()
                for pos in range(self.dat['clientes'] - 1):
                    if pos <= k[c] - 1:
                        ik1 = self.dat['D_ord'][c][k[c]]
                        ik = self.dat['D_ord'][c][pos]
                        Dk2 = self.dat['matriz_custo'][(ik1, c)]
                        Dk1 = self.dat['matriz_custo'][(ik, c)]
                        v[c].append(Dk2 - Dk1)
                    else:
                        # ik1 = self.dat['D_ord'][c][k[c]]
                        # Dk2 = self.dat['matriz_custo'][(ik1, c)]
                        v[c].append(0)
                        break
        else:
            v = dict()
            for c in range(self.dat['clientes']):
                v[c] = list()
                if k[c] == 0:
                    ik = self.dat['D_ord'][c][0]
                    Dk2 = self.dat['matriz_custo'][(ik, c)]
                    v[c].append(Dk2)
                else:
                    ik1 = self.dat['D_ord'][c][k[c]]
                    Dk2 = self.dat['matriz_custo'][(ik1, c)]
                    aux = 0
                    # aux = sum((Dk2 - self.dat['matriz_custo'][(i, c)]) * _y[i] for i in self.dat['D_ord'][c] if
                    #     self.dat['D_ord'][c].index(i) <= self.dat['D_ord'][c].index(k[c]))
                    for ll in range(len(self.dat['D_ord'][c])):
                        if ll > k[c]:
                            break
                        planta = self.dat['D_ord'][c][ll]
                        aux += (Dk2 - self.dat['matriz_custo'][(planta, c)]) * _y[planta]
                    vf = Dk2 - aux
                    v[c].append(vf)

        return round(sum(sum(v[c]) for c in range(self.dat['clientes']))), k, v


class BendersBranchAndCut():
    def __init__(self, dat):
        self.dat = dat
        self.create_master_problem()

    def create_master_problem(self):
        dat = self.dat
        m = Model(name='PMR', solver_name=CBC)
        m.verbose = 0
        I,J = range(dat['clientes']),range(dat['plantas']) #TODO: Não é melhor criar a classe DAT?

        y = [m.add_var(var_type=BINARY, name='y(%d)'%j) for j in J]
        eta = [m.add_var(var_type=CONTINUOUS, lb=0.0, name='eta(%d)' % i) for i in I]
        m.y = y
        m.eta = eta

        m.objective = minimize(sum(dat['custo_abertura'][j] * y[j] for j in J) + xsum(eta[i] for i in I))

        m += xsum(y[j] for j in J) >= 1  #Corte de viabilidade, já que solução viável precisa ter pelo menos 1 facility aberto
        #m.write('pmr.lp')
        self.m = m

    def run(self):
        ...
        m = self.m
        m.cuts_generator = BendersCutCutGenerator(self.dat, m.y, m.eta)
        m.lazy_constrs_generator = BendersCutCutGenerator(self.dat, m.y, m.eta)
        print('\n\n\n')
        m.verbose = 0
        start_time = time()
        status = m.optimize()
        self.run_time = time() - start_time
        if status == OptimizationStatus.OPTIMAL:
            self.is_solution = True
            #self.print_solution()
            return [a for a in m.y if a.x > 0], m.objective_value
    def print_solution(self):
        if self.is_solution == True:
            dat = self.dat
            m = self.m
        I, J = range(dat['clientes']),range(dat['plantas'])
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
