'''  IPA algorithm [1].
Scalable and Parallelizable Processing of Influence Maximization for Large-Scale Social Networks
'''

from __future__ import division
import networkx as nx
import math, time
from copy import deepcopy
from runIAC import avgIAC
import multiprocessing, json
from runIAC import avgIAC, runIAC

def updateAP(ap, S, IPAv, IPA_MIPv, Ep):

    # going from leaves to root
    sorted_MIPs = sorted(IPA_MIPv.iteritems(), key = lambda (_, MIP): len(MIP), reverse = True)
    for u, _ in sorted_MIPs:
        if u in S:
            ap[(u, IPAv)] = 1
        elif not IPAv.in_edges([u]):
            ap[(u, IPAv)] = 0
        else:
            in_edges = IPAv.in_edges([u], data=True)
            prod = 1
            for w, _, edata in in_edges:
                # p = (1 - (1 - Ep[(w, u)])**edata["weight"])
                p = Ep[(w,u)]
                prod *= 1 - ap[(w, IPAv)]*p
            ap[(u, IPAv)] = 1 - prod

def updateAlpha(alpha, v, S, IPAv, IPA_MIPv, Ep, ap):
    # going from root to leaves
    sorted_MIPs =  sorted(IPA_MIPv.iteritems(), key = lambda (_, MIP): len(MIP))
    for u, mip in sorted_MIPs:
        if u == v:
            alpha[(IPAv, u)] = 1
        else:
            out_edges = IPAv.out_edges([u])
            assert len(out_edges) == 1, "node u=%s must have exactly one neighbor, got %s instead" %(u, len(out_edges))
            w = out_edges[0][1]
            if w in S:
                alpha[(IPAv, u)] = 0
            else:
                in_edges = IPAv.in_edges([w], data=True)
                prod = 1
                for up, _, edata in in_edges:
                    if up != u:
                        # pp_upw = 1 - (1 - Ep[(up, w)])**edata["weight"]
                        pp_upw = Ep[(up, w)]
                        prod *= (1 - ap[(up, IPAv)]*pp_upw)
                alpha[(IPAv, u)] = alpha[(IPAv, w)]*(Ep[(u,w)])*prod

def computeIPAOA(G, u, theta, S, Ep):

    IPAOA = nx.DiGraph()
    IPAOA.add_node(u)
    IPAOA_MIP = {u: [u]}

    crossing_edges = set([out_edge for out_edge in G.out_edges([u]) if out_edge[1] not in S + [u]])
    edge_weights = dict()
    dist = {u: 0} # shortest paths from the root u


    while crossing_edges:
        # Dijkstra's greedy criteria
        min_dist = float("Inf")
        sorted_crossing_edges = sorted(crossing_edges) # to break ties consistently
        for edge in sorted_crossing_edges:
            if edge not in edge_weights:
                # edge_weights[edge] = -math.log(1 - (1 - Ep[edge])**G[edge[0]][edge[1]]["weight"])
                edge_weights[edge] = -math.log(Ep[edge])
            edge_weight = edge_weights[edge]
            if dist[edge[0]] + edge_weight < min_dist:
                min_dist = dist[edge[0]] + edge_weight
                min_edge = edge
        # check stopping criteria
        if min_dist < -math.log(theta):
            dist[min_edge[1]] = min_dist

            IPAOA.add_edge(min_edge[0], min_edge[1])
            IPAOA_MIP[min_edge[1]] = IPAOA_MIP[min_edge[0]] + [min_edge[1]]
            # update crossing edges
            crossing_edges.difference_update(G.in_edges(min_edge[1]))
            crossing_edges.update([out_edge for out_edge in G.out_edges(min_edge[1])
                                   if (out_edge[1] not in IPAOA) and (out_edge[1] not in S)])
        else:
            break
    return IPAOA, IPAOA_MIP

def updateIS(IS, S, u, IPAOA, IPA):
    for v in IPAOA[u]:
        for si in S:
            # if seed node is effective and it's blocked by u
            # then it becomes ineffective
            if (si in IPA[v]) and (si not in IS[v]) and (u in IPA[v][si]):
                    IS[v].append(si)

def computeIPA(G, ISv, v, theta, S, Ep):


    IPA = nx.DiGraph()
    IPA.add_node(v)
    IPA_MIP = {v: [v]}

    crossing_edges = set([in_edge for in_edge in G.in_edges([v]) if in_edge[0] not in ISv + [v]])
    edge_weights = dict()
    dist = {v: 0} # shortest paths from the root u


    while crossing_edges:

        min_dist = float("Inf")
        sorted_crossing_edges = sorted(crossing_edges) # to break ties consistently
        for edge in sorted_crossing_edges:
            if edge not in edge_weights:
                # edge_weights[edge] = -math.log(1 - (1 - Ep[edge])**G[edge[0]][edge[1]]["weight"])
                edge_weights[edge] = -math.log(Ep[edge])
            edge_weight = edge_weights[edge]
            if dist[edge[1]] + edge_weight < min_dist:
                min_dist = dist[edge[1]] + edge_weight
                min_edge = edge

        if min_dist < -math.log(theta):
            dist[min_edge[0]] = min_dist

            IPA.add_edge(min_edge[0], min_edge[1])
            IPA_MIP[min_edge[0]] = IPA_MIP[min_edge[1]] + [min_edge[0]]

            crossing_edges.difference_update(G.out_edges(min_edge[0]))
            if min_edge[0] not in S:
                crossing_edges.update([in_edge for in_edge in G.in_edges(min_edge[0])
                                       if (in_edge[0] not in IPA) and (in_edge[0] not in ISv)])
        else:
            break
    return IPA, IPA_MIP

def IPA(G, k, theta, Ep):
    start = time.time()

    S = []
    IncInf = dict(zip(G.nodes(), [0]*len(G)))
    IPA = dict() # node to tree
    IPAOA = dict()
    IPA_MIP = dict() # node to MIPs (dict)
    IPAOA_MIP = dict()
    ap = dict()
    alpha = dict()
    IS = dict()
    for v in G:
        IS[v] = []
        IPA[v], IPA_MIP[v] = computeIPA(G, IS[v], v, theta, S, Ep)
        for u in IPA[v]:
            ap[(u, IPA[v])] = 0
        updateAlpha(alpha, v, S, IPA[v], IPA_MIP[v], Ep, ap)
        for u in IPA[v]:
            IncInf[u] += alpha[(IPA[v], u)]*(1 - ap[(u, IPA[v])])
    print 'Finished initialization'
    print time.time() - start

    # main loop
    for i in range(k):
        u, _ = max(IncInf.iteritems(), key = lambda (dk, dv): dv)
        # print i+1, "node:", u, "-->", IncInf[u]
        IncInf.pop(u) # exclude node u for next iterations
        IPAOA[u], IPAOA_MIP[u] = computeIPAOA(G, u, theta, S, Ep)
        for v in IPAOA[u]:
            for w in IPA[v]:
                if w not in S + [u]:
                    IncInf[w] -= alpha[(IPA[v],w)]*(1 - ap[(w, IPA[v])])

        updateIS(IS, S, u, IPAOA_MIP, IPA_MIP)

        S.append(u)

        for v in IPAOA[u]:
            if v != u:
                IPA[v], IPA_MIP[v] = computeIPA(G, IS[v], v, theta, S, Ep)
                updateAP(ap, S, IPA[v], IPA_MIP[v], Ep)
                updateAlpha(alpha, v, S, IPA[v], IPA_MIP[v], Ep, ap)
                # add new incremental influence
                for w in IPA[v]:
                    if w not in S:
                        IncInf[w] += alpha[(IPA[v], w)]*(1 - ap[(w, IPA[v])])

    return S

def getCoverage((G, S, Ep)):
    return len(runIAC(G, S, Ep))

if __name__ == "__main__":
    import time
    start = time.time()

    model = "Categories"

    if model == "MultiValency":
        ep_model = "range"
    elif model == "Random":
        ep_model = "random"
    elif model == "Categories":
        ep_model = "degree"

    dataset = "gnu09"

    G = nx.read_gpickle("../../graphs/%s.gpickle" %dataset)
    print 'Read graph G'
    print time.time() - start

    Ep = dict()
    with open("Ep_%s_%s1.txt" %(dataset, ep_model)) as f:
        for line in f:
            data = line.split()
            Ep[(int(data[0]), int(data[1]))] = float(data[2])

    ALGO_NAME = "IPA"
    FOLDER = "Data4InfMax"
    SEEDS_FOLDER = "Seeds"
    TIME_FOLDER = "Time"
    DROPBOX_FOLDER = "/home/sergey/Dropbox/Influence Maximization"
    seeds_filename = SEEDS_FOLDER + "/%s_%s_%s_%s.txt" %(SEEDS_FOLDER, ALGO_NAME, dataset, model)
    time_filename = TIME_FOLDER + "/%s_%s_%s_%s.txt" %(TIME_FOLDER, ALGO_NAME, dataset, model)

    theta = 1.0/20
    pool = None
    I = 1000
    l2c = [[0, 0]]
    # open file for writing output
    seeds_file = open(seeds_filename, "a+")
    time_file = open(time_filename, "a+")
    dbox_seeds_file = open("%/%", DROPBOX_FOLDER, seeds_filename, "a+")
    dbox_time_file = open("%/%", DROPBOX_FOLDER, time_filename, "a+")
    for length in range(1, 250, 5):
        time2length = time.time()
        print "Start finding solution for length = %s" %length

        time2S = time.time()
        S = IPA(G, length, theta, Ep)
        time2complete = time.time() - time2S
        print >>time_file, (time2complete)
        print >>dbox_time_file, (time2complete)
        print 'Finish finding S in %s sec...' %(time2complete)

        print 'Writing S to files...'
        print >>seeds_filename, json.dumps(S)
        print >>dbox_seeds_file, json.dumps(S)


    seeds_file.close()
    dbox_seeds_file.close()
    time_file.close()
    dbox_time_file.close()
    print 'Total time: %s' %(time.time() - start)