import sys
from itertools import product
import os
import math
import numpy as np
import json
import random
from time import time
#from ilsvnd import *
#from dotmap import DotMap
from mip import Model, xsum, maximize, minimize, CBC, MAXIMIZE, MINIMIZE,OptimizationStatus, BINARY

same_limit = 10
h_limit = 500
# ------------------------------
#
# ------------------------------
def main(sys):
    dt = read_data(sys.argv[1])

    ub = mth(dt)

    lag_rel_assignment(dt,ub)

    lag_rel_knapsack(dt,ub)

    linear_gap(dt)
# ------------------------------
#
# ------------------------------
def lag_rel_assignment(dt,ub):
    print("Dualizing assignment constraints")
    J = dt.J
    I = dt.I
    c = dt.c
    b = dt.b
    t = dt.t
    stop = False
    same = 0
    lb = 0.0
    mu = 2.0
    h = 0
    u = [min([c[i][j] for j in J]) for i in I]
    # u = [0.0 for i in I]
    gap = 100.0
    while (stop == False):
        ua = u
    # solving Lagrangian problem
        x = [[0 for j in J] for i in I]
        grad = [1 for i in I]
        infimum = sum(u)
        for j in J:
            w = [(c[i][j] - u[i]) for i in I]
            s = [t[i][j] for i in I]
            of, _x = knapsack(w, s, b[j])
            infimum += of
            for i in _x:
                grad[i] -= 1
                x[i][j] = 1

        if (lb < infimum):
            lb = infimum
            same = 0
        else:
            same += 1

        if (same > same_limit):
            mu = mu / 2.0
            same = 0
    # sub gradient
        norm = sum([grad[i] ** 2 for i in I])
        gap = 100.0 * (ub - lb) / ub
        h += 1

        if (norm < 0.001) or (mu < 0.005) or (h > h_limit) or (gap < 0.001):
            stop = True
        else:
        # solving dual Lagrangian problem
            step = mu * (ub - lb) / norm
            u = [u[i] + step * grad[i] for i in I]
            ax = [abs(v1 - v2) for v1, v2 in zip(u, ua)]
        if (sum(ax) < 0.0001):
            stop = True

        print_stats(h, ub, infimum, lb, gap, mu, norm)

def knapsack(w,s,b):
    n = len(w)
    N = range(n)
    mod = Model('knapsack',sense=MINIMIZE,solver_name=CBC)
    x = [mod.add_var(var_type=BINARY, obj=w[i]) for i in N]
    mod += xsum(s[i] * x[i] for i in N) <= b
    mod.optimize()
    of = mod.objective_value
    x = [i for i, v in enumerate(x) if v.x > 0.9]
    return (of, x)
    # ------------------

def lag_rel_knapsack(dt,ub):
    print("Dualizing knapsack constraints")
    J = dt.J
    I = dt.I
    c = dt.c
    b = dt.b
    t = dt.t
    stop = False
    same = 0
    lb = 0.0
    mu = 2.0
    h = 0
    v = [0.0 for j in J]
    gap = 100.0
    while (stop == False):
        va = v
    # solving Lagrangian problem
        x = [[0 for j in J] for i in I]
        grad = [-b[j] for j in J]
        infimum = -sum([b[j] * v[j] for j in J])

        for i in I:
            of, j = min([((c[i][j] + v[j] * t[i][j]), j) for j in J], key=lambda t: t[0])

            infimum += of

            grad[j] += t[i][j]

            x[i][j] = 1
        if (lb < infimum):
            lb = infimum
            same = 0
        else:
            same += 1
        if (same > same_limit):
            mu = mu / 2.0
            same = 0

    # sub gradient
        norm = sum([grad[j] ** 2 for j in J])

        gap = 100.0 * (ub - lb) / ub
        h += 1
    if (norm < 0.0001) or (mu < 0.0005) or (h > h_limit) or (gap < 0.0001):
        stop = True
    else:
    # solving dual Lagrangian problem
        step = mu * (ub - lb) / norm
        v = [max(v[j] + step * grad[j], 0.0) for j in J]
        ax = [abs(v1 - v2) for v1, v2 in zip(v, va)]
        if (sum(ax) < 0.000001):
            stop = True
    print_stats(h, ub, infimum, lb, gap, mu, norm)


def linear_gap(dt):
    print("solving the gap using CBC")
    ni = dt.ni
    nj = dt.nj
    J = range(nj)
    I = range(ni)
    IJ = product(I, J)
    c = dt.c
    b = dt.b
    t = dt.t

    mod = Model('gap', sense=MINIMIZE, solver_name=CBC)
    mod.verbose = 0

    x = {(i, j): mod.add_var(var_type=BINARY, obj=c[i][j]) for (i, j) in IJ}
    for i in I:
        mod += xsum(x[i, j] for j in J) == 1
    for j in J:
        mod += xsum(t[i][j] * x[i, j] for i in I) <= b[j]

    mod.optimize(relax=True)
    lb = mod.objective_value
    # print(mod.solution.get_status_string())
    print("linear relaxation : {:.2f}".format(lb))
    mod.optimize(relax=False)
    ub = mod.objective_value
    print("solution          : {:.2f}".format(ub))
    print("linear gap        : {:.2f} %".format(100.0 * (ub - lb) / ub))


def mth(dt):
    ni = dt.ni
    nj = dt.nj
    J = dt.J
    I = list(dt.I)
    IJ = [(i, j) for j in J for i in I]
    c = dt.c
    b = list(dt.b)
    t = dt.t

    L = [sorted([(j, c[i][j]) for j in J], key=lambda t: t[1]) for i in I]
    x = [-1 for i in I]
    of = 0.0
    while len(I) > 0:
        i, w = max([(i, L[i][1][1] - L[i][0][1]) if len(L[i]) > 1 else (i, L[i][0][1]) for i in I], key=lambda t: t[1])
        j = L[i][0][0]
        if (b[j] - t[i][j] < 0):
            for i in I:
                for jj, v in L[i]:
                    if jj == j:
                        L[i].remove((jj, v))
        else:
            of += c[i][j]
            x[i] = j
            b[j] -= t[i][j]
            I.remove(i)

    I = dt.I
    best = (of, -1, -1, -1)
    stop = False
    while (stop == False):
        stop = True
        for i in I:
            for j in J:
                if x[i] != j:
                    dcap = b[j] - t[i][j]
                    if dcap >= 0:
                        delta = c[i][j] - c[i][x[i]]
                        if delta < 0:
                            stop = False
                            best = (of + delta, i, j, x[i])
        if stop == False:
            (of, i, j, oj) = best
            b[oj] += t[i][oj]
            b[j] -= t[i][j]
            x[i] = j
    print("ub : {:.2f}".format(best[0]))
    return best[0]


def print_stats(h, ub, infimum, lb, gap, mu, norm):
    print("{:4d} ".format(h), end='')
    print("{:12,.2f} ".format(ub), end='')
    print("{:12,.2f} ".format(infimum), end='')
    print("{:12,.2f} ".format(lb), end='')
    print("{:12,.2f} %".format(gap), end='')
    print("{:12,.8f} ".format(mu), end='')
    print("{:12,.2f} ".format(norm, end=''))
