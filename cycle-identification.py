# read stochiometric matrix
# read mapping of rows to metabolite names.

# initialize set of metabolites as m.

# to map cycles: for each metabolite X do:
# define the set of metabolites reachable from the initial metabolite, initialized as an empty set.
## for each reaction consuming X, add all of its products (with factors corresponding to stochiometries) to the set of metabolites reachable from X.
## for every added metabolite also include the set of reactions leading to it from X.
## remove the reaction from the network.
## recuresively repeat for all metabolites reachable from X until:
### no reactions are left leaving the set of reachable metabolites OR
### X is reached (then add path to X to list of identified cycles).
## remove X and all reactions consuming it from network and repeat for next metabolite.

#demo graph to test on:
n=4
E=[(0,1),(1,2),(2,3),(3,1),(3,0)]

blocked = set([])
B={}
AK={}
root = 0

def strong_connect_at(s,n,E):
   index = s
   S = []
   indices = {}
   lowlinks = {}
   def strongconnect(index,v):
#       print("strong connect")
#       print(v)
#       print(S)
#       print(indices)
#       print(lowlinks)
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
#       print(v)
#       print(S)
#       print(indices)
#       print(lowlinks)
       return ret
   return strongconnect(s,s)

def adjacency(n,E):
    adj = {}
    for i in range(n):
        adj[i]=[]
    for (s,t) in E:
        adj[s].append(t)
    return adj

def induced(vs,E):
    rc = {}
    for i in vs:
        rc[i]=[x for x in E[i] if x in vs]
    return rc

def write_cycle(x):
    print("cycle")
    print(x)
    f.write(str(x))
    f.write("\n")

def find_cycles(n,E): # n is the number of vertices in the graph. E is the set of directed edges (i,j) in the graph.
    E = adjacency(n,E)
    for s in range(n):
        Es = induced(range(s,n),E)
        vs = strong_connect_at(s,n,Es)
        Ak = induced(vs,Es)
        B = {}
        blocked = set()
        stack = []
        def unblock(v):
            if v in blocked:
                blocked.remove(v)
            if v in B:
                for w in B[v]:
                    B[v].remove(w)
                    if w in blocked:
                        unblock(w)

        def circuit(v):
            found = False
            stack.append(v)
            blocked.add(v)
            for w in Ak[v]:
                if w == s:
                    write_cycle(stack+[s])
                elif w not in blocked:
                    if circuit(w):
                        found = True
            if found:
                unblock(v)
            else:
                for w in Ak[v]:
                    if w not in B:
                        B[w]=[v]
                    elif v not in B[w]:
                        B[w].append(v)
            stack.pop()
            return found
        circuit(s)
    
print("find cycles:")
f = open("cycles.txt",'w')
find_cycles(n,E)
f.close()
