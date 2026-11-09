"""Check the cube seed conditions of the transfer theorem using the explicit
witnesses of the manuscript: role of each witness, signed states (P,N) of the
ten retained edges lie in the 40-state list, and joint separation J(r,s)."""
import itertools
Q=[(0,1),(0,3),(0,4),(1,2),(1,7),(2,3),(2,6),(3,5),(4,5),(4,7),(5,6),(6,7)]
E={frozenset(e) for e in Q}; p,q=frozenset((0,1)),frozenset((2,3)); R=E-{p,q}
def F(s): return {frozenset(map(int,(a,b))) for a,b in (x.split('-') for x in s)}
def path(seq,closed): 
    es=list(zip(seq,seq[1:]))+([(seq[-1],seq[0])] if closed else []); S={frozenset(e) for e in es}; assert S<=E and len(set(seq))==8; return S
W={'A':[path([0,1,2,3,5,6,7,4],1)],'B':[path([0,1,2,6,7,4,5,3],1)],'C':[path([0,3,2,1,7,6,5,4],1)],
   'D':[path([1,7,4,0,3,5,6,2],0)],'E':[path([0,4,7,1,2,6,5,3],0)]}
extra={'A1':[(0,3),(1,7),(2,6),(4,5),(4,7),(5,6)],'A2':[(0,4),(1,2),(3,5),(4,7),(5,6),(6,7)],'A3':[(0,4),(1,7),(2,6),(3,5),(4,5),(6,7)],'L':[(0,4),(1,7),(2,6),(3,5),(4,7),(5,6)]}
import networkx as nx
for k,es in extra.items():
    S={frozenset(e) for e in es}; G=nx.Graph(es); G.add_nodes_from(range(8))
    comps=[sorted(c) for c in nx.connected_components(G)]; ends=sorted(v for v in G if G.degree(v)==1)
    if k.startswith('A'):
        # with p,q added must be one Hamilton cycle; pairing after deletion (03)(12)
        H=nx.Graph(es+[(0,1),(2,3)]); assert all(d==2 for _,d in H.degree()) and nx.is_connected(H)
        pr=sorted(tuple(sorted(v for v in comp if v in (0,1,2,3))) for comp in comps)
        print(k,'C11, pairing',pr); assert pr==[(0,3),(1,2)]
        W['A'].append(S)
    else:
        pr=sorted(tuple(sorted(v for v in comp if v in (0,1,2,3))) for comp in comps)
        assert all(G.degree(v)<=2 for v in G) and len(comps)==2 and ends==[0,1,2,3]
        print('L pairing',pr); assert pr==[(0,1),(2,3)]; W['L']=[S]
roles='ABCDEL'
PN={}
for r in R:
    P=''.join(t for t in roles if any(r in w for w in W.get(t,[]))); N=''.join(t for t in roles if any(r not in w for w in W.get(t,[])))
    PN[r]=P+'/'+N
S40="""ABC/ADEL ABCD/AEL ABCE/ADL ABDEL/AC ACDEL/AB ADEL/ABCDL ADEL/ABCEL ABCDL/ADEL ABCEL/ADEL ABDEL/ACDL ABDEL/ACEL ACDEL/ABDL ACDEL/ABEL ABCDEL/ADL ABCDEL/AEL ADEL/ABCDEL ABCDL/ABDEL ABCDL/ACDEL ABCEL/ABDEL ABCEL/ACDEL ABDEL/ABCDL ABDEL/ABCEL ABDEL/ACDEL ACDEL/ABCDL ACDEL/ABCEL ACDEL/ABDEL ABCDEL/ABDL ABCDEL/ABEL ABCDEL/ACDL ABCDEL/ACEL ABCDEL/ADEL ABCDL/ABCDEL ABCEL/ABCDEL ABDEL/ABCDEL ACDEL/ABCDEL ABCDEL/ABCDL ABCDEL/ABCEL ABCDEL/ABDEL ABCDEL/ACDEL ABCDEL/ABCDEL""".split()
assert len(set(S40))==40
for r,s in PN.items(): print(sorted(r),s, s in S40)
allw=[w for t in W for w in W[t]]
Jok=all(any(r in w and s not in w for w in allw) for r in R for s in R if r!=s)
print('witnesses',len(allw),'J nonempty for all ordered pairs:',Jok)
assert len(allw)==9 and len({frozenset(w&R) for w in allw})==8 and Jok
assert all(PN[r] in S40 for r in R)
print('SEED_CHECKS_PASSED')
