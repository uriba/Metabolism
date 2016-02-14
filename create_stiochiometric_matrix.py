# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 14:25:45 2016

@author: yinonbaron
"""

import numpy as np
import sys
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible, revert_to_reversible
from cobra.core import Reaction, Metabolite, Formula
sys.setrecursionlimit(10000)

#model = create_cobra_model_from_sbml_file("iJO1366.xml")
#model.reactions.Ec_biomass_iJO1366_WT_53p95M.delete()
#model.reactions.Ec_biomass_iJO1366_core_53p95M.delete()
model = create_cobra_model_from_sbml_file("ecoli_core_model.xml")
model.reactions.Biomass_Ecoli_core_w_GAM.delete()
# add rubisco and prk:
def add_metabolite(model, id, formula, name, compartment='C'):
    try:
        model.metabolites.get_by_id(id)
    except KeyError:
        met = Metabolite(id=id, formula=Formula(formula),
                            name=name, compartment=compartment)
        model.add_metabolites([met])

def add_reaction(model, id, name, sparse,
                 lower_bound=0, upper_bound=1000):
    """
        Adds a reaction to the model
    """
    # convert the sparse representation using the metabolites in the model
    for key in sparse.keys():
        if key not in model.metabolites:
            raise Exception("cannot find the cytoplasmic metabolite %s in the model" % key)

    r = dict([(model.metabolites[model.metabolites.index(key)], val)
              for key, val in sparse.iteritems()])
    reaction = Reaction(name)
    reaction.id = id
    reaction.add_metabolites(r)
    reaction.lower_bound = lower_bound
    reaction.upper_bound = upper_bound
    model.add_reactions([reaction])
    return reaction
 

add_metabolite(model, 'rubp_D_c', 'C5H12O11P2', 'D-ribulose 1,5-bisphosphate')
sprs = {'ru5p_D_c' : -1, 'atp_c' : -1, 'rubp_D_c' : 1, 'adp_c' : 1}
add_reaction(model, 'PRK', 'phosphoribulokinase', sprs)
add_metabolite(model, 'rubp_D_c', 'C5H12O11P2', 'D-ribulose 1,5-bisphosphate')
sprs = {'rubp_D_c' : -1, 'h2o_c' : -1, 'co2_c' : -1, '3pg_c' : 2, 'h_c' : 3}
add_reaction(model, 'RBC', 'RuBisCO carboxylation', sprs)
 
convert_to_irreversible(model)
num_metabolites = len(model.metabolites)
num_reactions = len(model.reactions)
total_list = model.reactions + model.metabolites
nodes_mapping = dict()
map(lambda k, v: nodes_mapping.update({k: v}), total_list, range(0,num_metabolites+num_reactions))

network = set()

cofactors = ['CO2','O2','H2','ADP','ATP','dATP','GMP','H','H2O','AMP','UDP','GDP','GTP','CMP','dCTP','CTP','dUTP','UTP','ITP','IDP','UMP','dUMP','CDP','IMP','Hydrogen peroxide',
        'Ammonium','Phosphate','Coenzyme A','Coenzyme-A','Fe2+','H+','Fe3+','Diphosphate','Sodium',
'Nicotinamide adenine dinucleotide - reduced','Nicotinamide adenine dinucleotide', 
'Nicotinamide adenine dinucleotide phosphate', 'Nicotinamide adenine dinucleotide phosphate - reduced',
'Nicotinamide-adenine-dinucleotide-reduced','Nicotinamide-adenine-dinucleotide',
'Nicotinamide-adenine-dinucleotide-phosphate-reduced','Nicotinamide-adenine-dinucleotide-phosphate',
'Ubiquinol-8','Ubiquinone-8',
'Flavin adenine dinucleotide reduced',
'Flavin adenine dinucleotide oxidized']

cycle_classes = {   'Pentose-phosphate':['D-Erythrose-4-phosphate','Sedoheptulose-7-phosphate'],
                    'Calvin Benson':['RuBisCO carboxylation'],
                    'FBA':['FBA_reverse'],
                    'PTS':['D-glucose transport via PEP:Pyr PTS','Fructose transport via PEP:Pyr PTS (f6p generating)'],
                    'TCA':['Acetyl-CoA']}

onestep = []
for reaction in model.reactions:
    reactants = []
    products = []
    for reactant in reaction.reactants:
        if reactant.name not in cofactors:
            network.add((nodes_mapping[reactant],nodes_mapping[reaction]))
            reactants.append(nodes_mapping[reactant])
    for product in reaction.products:
        if product.name not in cofactors:
            network.add((nodes_mapping[reaction],nodes_mapping[product]))
            products.append(nodes_mapping[product])
    for reactant in reactants:
        for product in products:
            onestep.append((reactant,(product,reaction)))

        
n= num_metabolites + num_reactions
E= network
print(len(network))

def adjacency(n,E):
    adj = {}
    for i in range(n):
        adj[i]=[]
    for (s,t) in E:
        adj[s].append(t)
    return adj

linked = adjacency(n,onestep)

def strong_connect_at(s,n,E):
   index = s
   S = []
   indices = {}
   lowlinks = {}
   def strongconnect(index,v):
       indices[v]=index
       lowlinks[v]=index
       index+=1
       S.append(v)
       for u in E[v]:
           if u not in indices:
               strongconnect(index,u)
               lowlinks[v] = min(lowlinks[v],lowlinks[u])
           elif u in S:
               lowlinks[v] = min(lowlinks[v],indices[u])
       ret = [] 
       if lowlinks[v]==indices[v]:
           while S[-1] != v:
               ret.append(S.pop())
           ret.append(S.pop())
       print("ret")
       print(ret)
       if len(ret)>1:
           path = [reverse_nodes_mapping[i].name for i in ret]
           print(path)
           #y=raw_input("continue")
       return ret
   return strongconnect(s,s)

def induced(vs,E):
    rc = {}
    for i in vs:
        rc[i]=[x for x in E[i] if x in vs]
    return rc

count = 0
gcount=0
cycles = set()
autocycles = list()

def sorted_cycles():
    count = 0
    global autocycles
    global f3
    autocycles = sorted(autocycles,key=len)
    for x in autocycles:
        path = [reverse_nodes_mapping[i].name for i in x]
        reactions = set()
        metabolites = set()
        stoich = {}
        for i in x:
            i = reverse_nodes_mapping[i]
            if 'metabolites' in dir(i):
                reactions.add(i)
                for reactant in i.reactants:
                    stoich[reactant]=0
                for product in i.products:
                    stoich[product]=0
            else:
                metabolites.add(i)
        for reaction in reactions:
            for m in reaction.metabolites:
                stoich[m]+=reaction.metabolites[m]
        count+=1
        f3.write('Cycle %d:\n' % count)
        if 'metabolites' in dir(reverse_nodes_mapping[x[0]]):
            path = path[1:]
        else:
            path = path[:-1]
        for i in range(len(path)/2):
            f3.write("%-30s -->\t %-30s\n" %(path[2*i],path[2*i+1]))
        inputs = []
        outputs = []
        for m in stoich:
            if m.name not in cofactors:
                if stoich[m]>0:
                    outputs.append("%s: %d" % (m.name,stoich[m]))
                if stoich[m]<0:
                    inputs.append("%s: %d" % (m.name,stoich[m]))
        f3.write('Inputs: %s\n' % str(inputs))
        f3.write('Outputs: %s\n' % str(outputs))
        f3.write('---------------------------------\n')


def process_cycle(x):
    global count
    global gcount
    global f
    global f2
    global autocycles
    if len(x)>5:
        if len(set(x))<len(x)-1:
            print(x)
            print(set(x))
            print "bug"
            f2.close()
            f.close()
            exit
        path = [reverse_nodes_mapping[i].name for i in x]
        reactions = set()
        metabolites = set()
        stoich = {}
        for i in x:
            i = reverse_nodes_mapping[i]
            if 'metabolites' in dir(i):
                reactions.add(i)
                for reactant in i.reactants:
                    stoich[reactant]=0
                for product in i.products:
                    stoich[product]=0
            else:
                metabolites.add(i)
        if frozenset(reactions) not in cycles:
            forbidden = set()
            for reaction in reactions:
                for m in reaction.metabolites:
                    stoich[m]+=reaction.metabolites[m]
                if reaction.id.endswith('_reverse'):
                    forbidden.add(reaction.id[:-8])
                for m in reaction.products:
                    if m not in metabolites and m.name not in cofactors:
                        #for (p,r) in linked[nodes_mapping[m]]:
                        #    if reverse_nodes_mapping[p] in metabolites and r not in reactions:
                        #        print('extra')
                        #        print((m,p,r))
                        #        for n in r.metabolites:
                        #            if n not in stoich:
                        #                stoich[n]=r.metabolites[n]
                        #            else:
                        #                stoich[n]+=r.metabolites[n]
                        #        break
                                

                        for n in metabolites:
                            if nodes_mapping[n] in zip(*linked[nodes_mapping[m]])[0]:
                                stoich[m]-=1
                                stoich[n]+=1
            omit = False
            for r in reactions:
                if r.id in forbidden:
                    omit = True
            auto = False
            neg = False
            for metabolite in metabolites:
                if stoich[metabolite] < 0:
                    neg = True
                if stoich[metabolite] > 0:
                    auto = True
            if auto or neg: # check if adding more reactions can rescue the cycle
                for reaction in reactions:
                    for m in reaction.products:
                        if m not in metabolites and m.name not in cofactors:
                            for (p,r) in linked[nodes_mapping[m]]:
                                if reverse_nodes_mapping[p] in metabolites and r not in reactions:
                                    print('extra')
                                    print((m,p,r))
                                    for n in r.metabolites:
                                        if n not in stoich:
                                            stoich[n]=r.metabolites[n]
                                        else:
                                            stoich[n]+=r.metabolites[n]
                                    break
            auto = False
            neg = False
            for metabolite in metabolites:
                if stoich[metabolite] < 0:
                    neg = True
                if stoich[metabolite] > 0:
                    auto = True
            if auto and not neg and not omit:
                print("auto-catalytic cycle %d" % gcount)
                autocycles.append(x)
                count+=1
                gcount+=1
                cycles.add(frozenset(reactions))
                f2.write('Cycle %d:\n' % gcount)
                if 'metabolites' in dir(reverse_nodes_mapping[x[0]]):
                    path = path[1:]
                else:
                    path = path[:-1]
                print(str(path))
                for i in range(len(path)/2):
                    f2.write("%-30s -->\t %-30s\n" %(path[2*i],path[2*i+1]))
                inputs = []
                outputs = []
                for m in stoich:
                    if m.name not in cofactors:
                        if stoich[m]>0:
                            outputs.append("%s: %d" % (m.name,stoich[m]))
                        if stoich[m]<0:
                            inputs.append("%s: %d" % (m.name,stoich[m]))
                f2.write('Inputs: %s\n' % str(inputs))
                f2.write('Outputs: %s\n' % str(outputs))
                cat = 'New!'
                for c in cycle_classes:
                    for m in stoich:
                        if stoich[m]<0 and m.name in cycle_classes[c]:
                            cat = c
                    for r in reactions:
                        if r.name in cycle_classes[c]:
                            cat = c
                f2.write('category: %s' % cat)

                f2.write('---------------------------------\n')

        if count>100:
            f.close()
            f2.close()
            if gcount>100000:
                exit
            f = open("/home/uri/cycles.txt",'a')
            f2 = open("/home/uri/autocycles.txt",'a')
            count=0
        f.write(str(path))
        f.write("\n")

def find_cycles(n,E): # n is the number of vertices in the graph. E is the set of directed edges (i,j) in the graph.
    E = adjacency(n,E)
    def unblock(B,blocked,v):
        if v in blocked:
            blocked.remove(v)
        if v in B:
            for w in list(B[v]):
                B[v].remove(w)
                if w in blocked:
                    unblock(B,blocked,w)

    def circuit(Ak,B,stack,blocked,v):
        found = False
        stack.append(v)
        blocked.add(v)
        for w in Ak[v]:
            if w == s:
                process_cycle(stack+[s])
                found = True
            elif w not in blocked:
                if circuit(Ak,B,stack,blocked,w):
                    found = True
        if found:
            unblock(B,blocked,v)
        else:
            for w in Ak[v]:
                if v not in B[w]:
                    B[w].append(v)
        stack.pop()
        return found

    for s in range(n):
        Es = induced(range(s,n),E)
        vs = strong_connect_at(s,n,Es)
        Ak = induced(vs,Es)
        B = {}
        for m in range(s,n):
            B[m]=[]
        blocked = set()
        stack = []
        circuit(Ak,B,stack,blocked,s)

reverse_nodes_mapping = dict()
map(lambda k, v: reverse_nodes_mapping.update({k: v}), nodes_mapping.values(), nodes_mapping.keys())

    
print("find cycles:")
##################################

A = adjacency(n,E)
def print_point_stack():
    for p in point_stack:
        print p,
    print

point_stack = list()
marked = dict()
marked_stack = list()

def backtrack(v):
    f = False
    point_stack.append(v)
    marked[v] = True
    marked_stack.append(v)
    for w in A[v]:
        if w<s:
            A[w] = 0
        elif w==s:
#            print_point_stack()
            process_cycle(point_stack+[s])
            f = True
        elif not marked[w]:
            f = backtrack(w) or f
    if f:
        while marked_stack[-1] != v:
            u = marked_stack.pop()
            marked[u] = False
        marked_stack.pop()
        marked[v] = False
    point_stack.pop()
    return f

for i in range(len(A)):
    marked[i] = False

f = open("/home/uri/cycles.txt",'w')
f2 = open("/home/uri/autocycles.txt",'w')
f3 = open("/home/uri/sortedautocycles.txt",'w')
#find_cycles(n,E)
for s in range(len(A)):
    backtrack(s)
    while marked_stack:
        u = marked_stack.pop()
        marked[u] = False
sorted_cycles()
f.close()
f2.close()
f3.close()
