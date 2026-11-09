"""Brute-force checks of the CEN/UV cube family (finite depths only; not a proof).
For every construction word up to a given depth: Barnette class (simple, cubic,
bipartite, planar, 3-connected) and global H^{+-} by exhaustive Hamilton-cycle
enumeration; isomorphism-class counts, rooted (old-core-fixed) type counts, and
comparison with the lower bound ceil(2^{max(k-2,0)}/(n*3^7))."""
import itertools, math, sys, json, gc, platform
import networkx as nx
from networkx.algorithms import isomorphism as iso

CUBE=[(0,1),(0,3),(0,4),(1,2),(1,7),(2,3),(2,6),(3,5),(4,5),(4,7),(5,6),(6,7)]

def expand(G, site, nxt):
    u,v,x,y=site
    assert G.has_edge(u,v) and G.has_edge(x,y)
    a,b,c,d=nxt,nxt+1,nxt+2,nxt+3
    G.remove_edge(u,v); G.remove_edge(x,y)
    G.add_edges_from([(u,a),(a,b),(b,v),(x,c),(c,d),(d,y),(a,d),(b,c)])
    return {'CEN':(b,c,d,a),'UV':(u,a,b,v)}

def build(word, refl=None):
    """word: letters for successor choice after each expansion (len k-1 => depth k)."""
    G=nx.Graph(CUBE); site=(0,1,2,3); nxt=8
    refl=refl or [False]*len(word)
    succ=expand(G,site,nxt); nxt+=4
    for L,r in zip(word,refl):
        s=succ[L]
        if r: s=(s[2],s[3],s[0],s[1])
        succ=expand(G,s,nxt); nxt+=4
    return G

def barnette(G):
    return (all(d==3 for _,d in G.degree()) and nx.is_bipartite(G)
            and nx.check_planarity(G)[0] and nx.node_connectivity(G)>=3)

def ham_cycles(G):
    nodes=sorted(G); idx={e:i for i,e in enumerate(frozenset(e) for e in G.edges())}
    adj={v:list(G[v]) for v in nodes}; n=len(nodes); s=nodes[0]
    out=[]; path=[s]; used={s}
    def eid(p,q): return idx[frozenset((p,q))]
    def dfs(v,mask):
        if len(path)==n:
            if G.has_edge(v,s):
                m=mask|(1<<eid(v,s))
                out.append(m)
            return
        for w in adj[v]:
            if w not in used:
                used.add(w); path.append(w); dfs(w,mask|(1<<eid(v,w))); path.pop(); used.remove(w)
    sys.setrecursionlimit(10000); dfs(s,0)
    # each cycle found twice (two directions); dedupe
    return list(set(out)), len(idx)

def global_hpm(G):
    cyc,m=ham_cycles(G); full=(1<<m)-1
    for e in range(m):
        acc=0
        for c in cyc:
            if c>>e&1: acc|=full&~c
        if acc!=full&~(1<<e): return False,len(cyc)
    return True,len(cyc)

def rooted_classes(gs):
    """classes up to isomorphism fixing the 8 old cube vertices pointwise"""
    reps=[]
    for G in gs:
        for v in G: G.nodes[v]['lab']=v if v<8 else -1
        nm=iso.categorical_node_match('lab',-1)
        if not any(nx.is_isomorphic(G,H,node_match=nm) for H in reps): reps.append(G)
    return len(reps)

def classes(gs):
    reps={}
    for G in gs:
        h=nx.weisfeiler_lehman_graph_hash(G,iterations=5)
        L=reps.setdefault(h,[])
        if not any(nx.is_isomorphic(G,H) for H in L): L.append(G)
    return sum(len(v) for v in reps.values())

if __name__=='__main__':
    maxk=int(sys.argv[1]) if len(sys.argv)>1 else 6
    res=[dict(python=platform.python_version(),networkx=nx.__version__)]
    print(res[0],flush=True)
    Q=nx.Graph(CUBE); ok,nc=global_hpm(Q)
    assert ok and barnette(Q)
    res.append(dict(k=0,n=8,words=1,barnette=barnette(Q),hpm_all=ok,classes=1,rooted=1,bound=1))
    for k in range(1,maxk+1):
        words=list(itertools.product(['CEN','UV'],repeat=k-1)) or [()]
        gs=[build(w) for w in words]
        b=all(barnette(G) for G in gs)
        h=all(global_hpm(G)[0] for G in gs)
        n=8+4*k; bound=math.ceil(2**max(k-2,0)/(n*3**7))
        r=dict(k=k,n=n,words=len(words),barnette=b,hpm_all=h,classes=classes(gs),
               rooted=rooted_classes(gs),rooted_pred=2**max(k-2,0),bound=bound)
        print(r,flush=True); res.append(r)
        assert b and h and r['rooted']==r['rooted_pred'] and r['classes']>=bound, r
        del gs; gc.collect()
    # reflected words (all 4^(k-1) choices) for small depth
    for k in range(1,min(maxk,5)+1):
        allok=True; cnt=0
        for w in itertools.product(['CEN','UV'],repeat=k-1):
            for rf in itertools.product([False,True],repeat=k-1):
                G=build(w,list(rf)); cnt+=1
                allok&=barnette(G) and global_hpm(G)[0]
                del G; gc.collect()
        print(dict(reflected_k=k,graphs=cnt,all_barnette_and_hpm=allok),flush=True)
        res.append(dict(reflected_k=k,graphs=cnt,ok=allok))
        assert allok and cnt==4**(k-1), (k,cnt)
    json.dump(res,open('check_family_results.json','w'),indent=1)
    print('FAMILY_CHECKS_PASSED')
