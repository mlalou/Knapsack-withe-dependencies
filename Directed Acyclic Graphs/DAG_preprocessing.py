import networkx as nx

print('Give the number of items:')
nbr=input()

my_coef = [0.1, 0.2, 0.5, 0.8]
n = [x * float(nbr) for x in my_coef]

for x in range(int(4)):
    # print("-------------------- Generate random DAG --------------------")
    #DAG=nx.read_edgelist("graph_dependencies_20_1+.txt")
    DAG = nx.read_edgelist('dependence_graph_20.txt', create_using=nx.DiGraph)
    DAG = nx.relabel_nodes(DAG, {a:int(a) for a in DAG.nodes()})
    print(list(DAG.edges))
    
    # #// Get file name
    # name_g="graph_dependencies_"+str(int(n[x]))+"_"+str(x)+".txt"
    
    # #// Generate random DAG       
    # import random
    # G=nx.gnp_random_graph(int(n[x]),0.5,directed=True)
    # DAG = nx.DiGraph([(u,v) for (u,v) in G.edges() if u<v])
    # nx.write_edgelist(DAG, name_g)
    
    if nx.is_directed_acyclic_graph(DAG): print('true')
    else: print ('false')

    print("-------------------- Transisitve_reduction --------------------")
    #// Get file name
    name_tr="transitive_reduction_"+str(int(n[x]))+"_"+str(x)+".txt"
    
    TR = nx.transitive_reduction(DAG)
    name="_".join([name_tr, str(x)])
    print(name)
    nx.write_edgelist(TR, name_tr)
    xx=list(TR.edges)
    print(xx)
    
    print("-------------------- Topological sort -------------------------")
    #// Get file name
    name_ts="topological_sort_"+str(int(n[x]))+"_"+str(x)+".txt"
    
    TS = nx.topological_sort(TR)
    #nx.write_edgelist(TS, "test.topological_sort")
    TopS=list(TS)
    print(TopS)
    with open(name_ts, 'w') as fp:
        fp.write("\n".join(str(item) for item in TopS))
    
 
     
    print("-------------------- Subtrees --------------------------------")
    #// Get file name
    name_st="subtrees_"+str(int(n[x]))+"_"+str(x)+".txt"
    
    lst_subtrees = []
    
    i = 1
    while i < int(n[x]-2):
        nx.descendants(TR,i)
        i += 1
        
    nbr=input()
    
    i = 1
    while i < int(n[x]-1):
        lst_subtrees.append(sorted(list(nx.descendants(TR, i))))
        i += 1

    print(lst_subtrees)

    with open(name_st, 'w') as fst:
        fst.write("\n".join(str(item) for item in lst_subtrees))
    #add a "\n" at the end according to the code in C++
    f = open(name_st, "a+")
    f.write("\n")
    
    #li=[]
    #with open(r'test.subtrees.txt', 'r') as fst:
    #    li=fst_subtrees.read()
    #    print(li)

    print("-------------------- Child nodes ------------------------------")
    #// Get file name
    name_ch="child_nodes_"+str(int(n[x]))+"_"+str(x)+".txt"
    
    lst_children = []
    i = 0
    while i < int(n[x]):
      lst_children.append(list(TR.successors(i)))
      i += 1

    with open(name_ch, 'w') as fch:
        fch.write("\n".join(str(item) for item in lst_children))

    f = open(name_ch, "a+")
    f.write("\n")
    print(lst_children)
            
    


print('Done')


