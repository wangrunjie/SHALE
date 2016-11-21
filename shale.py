#!/usr/bin/env python

from sympy import *
from random import random
from copy import copy

class Supply:
    def __init__(self, file):
        self.pair = {}
        self.satisfy_demand = {}
        with open(file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                i, s = line.split('\t')
                self.pair[i] = int(s)

    def get_supply(self, i):
        return self.pair[i]

    def get_satisfy_demand(self, i):
        return self.satisfy_demand[i]

    def get_all_i(self):
        return self.pair.keys()


class Demand:
    def __init__(self, file, supply):
        self.demand = {}
        self.penalty = {}
        self.target_supply = {}
        with open(file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                j, d, p, ii = line.split('\t')
                self.demand[j] = int(d)
                self.penalty[j] = float(p)
                self.target_supply[j] = ii.split(',')
        self._set_supply_satisfy_demand(supply)

    def _set_supply_satisfy_demand(self, supply):
        for (j, ii) in self.target_supply.items():
            for i in ii:
                if i not in supply.satisfy_demand:
                    supply.satisfy_demand[i] = []
                supply.satisfy_demand[i].append(j)

    def get_demand(self, j):
        return self.demand[j]

    def get_penalty(self, j):
        return self.penalty[j]

    def get_target_supply(self, j):
        return self.target_supply[j]
    
    def get_all_j(self):
        return self.demand.keys()

    def get_v(self, j):
        return 1.0


class Shale:
    def __init__(self, supply, demand):
        self.supply = supply
        self.demand = demand
        self.solve = Solve()

    def initialize(self):
        self.alpha_j = {}
        self.beta_i = {}
        self.theta_ij = {}
        self.sigma_j = {}
        for j in self.demand.get_all_j():
            self.alpha_j[j] = 0.0
            sum = 0.0
            for i in self.demand.get_target_supply(j):
                sum += self.supply.get_supply(i)
            self.theta_ij[j] = self.demand.get_demand(j) / sum
        

    def stage_one(self, iters):
        while iters > 0:
            for i in self.supply.get_all_i():
                self.update_beta(i)
            for j in self.demand.get_all_j():
                self.update_alpha(j)
            iters -= 1
        print 'stage one alpha --->', sorted(self.alpha_j.iteritems(), key=lambda d:d[0])
        print 'stage one beta --->', sorted(self.beta_i.iteritems(), key=lambda d:d[0])

    def stage_two(self):
        self.s_i = {}
        for i in self.supply.get_all_i():
            self.s_i[i] = self.supply.get_supply(i)
            self.update_beta(i)
        print 'stage two beta --->', sorted(self.beta_i.iteritems(), key=lambda d:d[0])
        sigma = {}
        for j in self.demand.get_all_j():
            #self.find_sigma(j)
            self.find_sigma2(j)
            #self.update_sigma(j)
            for i in self.demand.get_target_supply(j):
                g = max(0.0, self.theta_ij[j] * (1.0+(self.sigma_j[j]-self.beta_i[i])/self.demand.get_v(j)))
                self.s_i[i] -= min(self.s_i[i], self.supply.get_supply(i) * g)

    def output(self):
        print 'output alpha_j --->', sorted(self.alpha_j.iteritems(), key=lambda d:d[0])
        print 'output sigma_j --->', sorted(self.sigma_j.iteritems(), key=lambda d:d[0])

    def update_beta(self, i):
        coef = []
        for j in self.supply.get_satisfy_demand(i): 
            a = self.theta_ij[j] * (1.0 + self.alpha_j[j] / self.demand.get_v(j))
            b = self.theta_ij[j] / self.demand.get_v(j)
            coef.append((a / b, a, b))
        result = self.solve.max(coef, 1.0)
        if len(result) == 0 or result[0] < 0.0:
            self.beta_i[i] = 0.0
        else:
            self.beta_i[i] = result[0]

    def update_alpha(self, j):
        coef = []
        for i in self.demand.get_target_supply(j):
            s = self.supply.get_supply(i)
            a = -1.0 * s * self.theta_ij[j] * (1.0 - self.beta_i[i] / self.demand.get_v(j))
            b = -1.0 * s * self.theta_ij[j] /self.demand.get_v(j)
            coef.append((a / b, a, b))
        result = self.solve.max(coef, -1.0 * self.demand.get_demand(j))
        if len(result) == 0 or -1.0 * result[0] > self.demand.get_penalty(j):
            self.alpha_j[j] = self.demand.get_penalty(j)
        #if len(result) == 0 or result[0] < 0.0:
        #    self.alpha_j[j] = 0.0
        else:
            self.alpha_j[j] = -1.0 * result[0]

    def find_sigma(self, j):
        sigma = Symbol('sigma')
        flag = True
        for i in self.demand.get_target_supply(j):
            s = self.supply.get_supply(i)
            f = self.theta_ij[j] * (1.0 + (sigma - self.beta_i[i]) / self.demand.get_v(j))
            f_max = s * Piecewise((f, f >= 0), (0, f < 0))
            f_min = Piecewise((self.s_i[i], self.s_i[i] <= f_max), (f_max, self.s_i[i] > f_max))
            if flag:
                sum = f_min
                flag = False
            else:
                sum += f_min
        result = solve(sum - self.demand.get_demand(j), sigma)
        if len(result) == 0:
            self.sigma_j[j] = float('inf')
        else:
            self.sigma_j[j] = result[0]

    def find_sigma2(self, j):
        result = self.update_sigma(j)
        if len(result) == 0:
            self.sigma_j[j] = float('inf')
        else:
            self.sigma_j[j] = -1.0 * result[0]
            #self.sigma_j[j] = result[0]

    def update_sigma(self, j):
        solutions = []
        coef = []
        for i in self.demand.get_target_supply(j):
            s = self.supply.get_supply(i)
            a = s * self.theta_ij[j] * (1.0 - self.beta_i[i] / self.demand.get_v(j))
            b = s * self.theta_ij[j] / self.demand.get_v(j)
            coef.append(((self.s_i[i] - a) / b, a / b, a, b, self.s_i[i]))
        coef = sorted(coef, key=lambda t:t[1])

        for k in range(len(coef)):
            temp = copy(coef[k:])
            temp = sorted(temp, key=lambda t:t[0])
            sum_remained = 0.0
            sum_cons = 0.0
            sum_coef = 0.0
            for t in range(len(temp)):
                sum_cons += temp[t][2]
                sum_coef += temp[t][3]
                sum_remained += temp[t][4]
            res = (sum_cons - self.demand.get_demand(j)) / sum_coef
            if sum_remained <  self.demand.get_demand(j):
                continue
            if k == 0:
                if res <=temp[0][0] and res <= coef[k][1]:
                    solutions.append(res)
            else:
                if res <=temp[0][0] and res <= coef[k][1] and res >= coef[k-1][1]:
                    solutions.append(res)
            sum_remained = 0.0
            for t in range(1, len(temp)):
                sum_cons -= temp[t-1][2]
                sum_coef -= temp[t-1][3]
                sum_remained += temp[t-1][4]
                res = (sum_cons + sum_remained - self.demand.get_demand(j)) / sum_coef
                if k == 0:
                    if res <= temp[t][0] and res >= temp[t-1][0] and res <= coef[k][1]:
                        solutions.append(res)
                else:
                    if res <= temp[t][0] and res >= temp[t-1][0] and res <= coef[k][1] and res >= coef[k-1][1]:
                        solutions.append(res)
        return solutions


class Online:
    def __init__(self, supply, demand, alpha_j, sigma_j):
        self.supply = supply
        self.demand = demand
        self.alpha_j = alpha_j
        self.sigma_j = sigma_j
        self.theta_ij = {}
        self.beta_i = {}
        self.allocation_j = {}
        self.remaind_i = {}
        for i in self.supply.get_all_i():
            self.remaind_i[i] = supply.get_supply(i)
        for j in self.demand.get_all_j():
            sum = 0.0
            for i in self.demand.get_target_supply(j):
                sum += self.supply.get_supply(i)
            self.theta_ij[j] = self.demand.get_demand(j) / sum
            self.allocation_j[j] = 0
        self.solve = Solve()


    def allocation(self, i):
        s = 1.0
        x_ij = {}
        if i not in self.beta_i:
            self.update_beta(i)
        for j in self.supply.get_satisfy_demand(i):
            g = max(0.0, self.theta_ij[j] * (1.0+(self.sigma_j[j]-self.beta_i[i])/self.demand.get_v(j)))
            x_ij[j] = min(s, g)
            s -= x_ij[j]
        
        sum = 0.0
        for (j, p) in x_ij.items():
            sum += p
        if sum < 1.0:
            print 'there is %f chance that no conract is selected' % (1.0 - sum)

        r = random()
        sum = 0.0
        for (j, p) in x_ij.items():
            sum += p
            if r < sum:
                self.allocation_j[j] += 1
                self.remaind_i[i] -= 1
                break


    def update_beta(self, i):
        coef = []
        for j in self.supply.get_satisfy_demand(i): 
            a = self.theta_ij[j] * (1.0 + self.alpha_j[j] / self.demand.get_v(j))
            b = self.theta_ij[j] / self.demand.get_v(j)
            coef.append((a / b, a, b))
        result = self.solve.max(coef, 1.0)
        if len(result) == 0 or result[0] < 0.0:
            self.beta_i[i] = 0.0
        else:
            self.beta_i[i] = result[0]


class Solve:
    # http://math.stackexchange.com/questions/145458/solve-equations-using-the-max-function
    def max(self, coef, y):
        solutions = []
        sum_cons = 0.0
        sum_coef = 0.0
        for i in range(len(coef)):
            sum_cons += coef[i][1]
            sum_coef += coef[i][2]
        coef = sorted(coef, key=lambda t:t[0])
        res = (sum_cons - y) / sum_coef
        if res <= coef[0][0]:
            solutions.append(res)
        for i in range(1, len(coef)):
            sum_cons -= coef[i-1][1]
            sum_coef -= coef[i-1][2]
            res = (sum_cons - y) / sum_coef
            if res <= coef[i]:
                solutions.append(res)
                break
        return solutions

    def max2(self):
        pass


    #http://www.bitsofpancake.com/math/minimum-and-maximum-of-two-functions/
    #http://math.stackexchange.com/questions/602553/how-to-invert-max-and-min-operators-in-equations
    def minmax(self):
        pass
                


class Debug:
    def __init__(self, shale, online):
        self.shale = shale
        self.online = online

    def print_supply(self):
        ii = self.shale.supply.get_all_i()
        ii.sort()
        print "\nsupply:"
        print "supply_node\tinventory\tsatisfy_demand"
        for i in ii:
            print '%s\t\t%d\t\t%s' % (i, self.shale.supply.get_supply(i), ','.join(self.shale.supply.get_satisfy_demand(i)))

    def print_demand(self):
        jj = self.shale.demand.get_all_j()
        jj.sort()
        print "\ndemand:"
        print "demand_node\tdemand\tpenalty\ttarget_supply"
        for j in jj:
            print '%s\t\t%d\t%f\t%s' % (j, self.shale.demand.get_demand(j), self.shale.demand.get_penalty(j), \
                                        ','.join(self.shale.demand.get_target_supply(j)))

    def print_online_allocation(self):
        jj = self.online.demand.get_all_j()
        jj.sort()
        print "\nallocation:"
        print "demand_node\tdemand\t\tallocation"
        for j in jj:
            print '%s\t\t%d\t\t%d' % (j, self.online.demand.get_demand(j), self.online.allocation_j[j])

    def print_online_remaind(self):
        ii = self.online.supply.get_all_i()
        ii.sort()
        print "\nremaind:"
        print "supply_node\tinventory\tremaind"
        for i in ii:
            print '%s\t\t%d\t\t%s' % (i, self.online.supply.get_supply(i), self.online.remaind_i[i])


def main():
    supply = Supply('./supply.txt')
    demand = Demand('./demand.txt', supply)

    shale = Shale(supply, demand)
    shale.initialize()
    shale.stage_one(5)
    shale.stage_two()
    shale.output()

    online = Online(supply, demand, shale.alpha_j, shale.sigma_j)
    for i in supply.get_all_i():
        inventory = supply.get_supply(i)
        while inventory > 0:
            online.allocation(i)
            inventory -= 1

    debug = Debug(shale, online)
    debug.print_supply()
    debug.print_demand()
    debug.print_online_allocation()
    debug.print_online_remaind()



if __name__ == '__main__' :
    main()
