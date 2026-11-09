"""Finite checks of the local tables: gadget traversal patterns (Table 1),
seven minimal covers, cube witnesses and their roles, the cap isomorphism."""
import itertools, networkx as nx
E=['ua','ab','bv','xc','cd','dy','ad','bc']
pats=[]
for bits in itertools.product([0,1],repeat=8):
    S=[e for e,b in zip(E,bits) if b]
    if all(sum(z in e for e in S)==2 for z in 'abcd'): pats.append(''.join(map(str,bits)))
print('degree-2 patterns:',pats)
table={'D11':'11111100','D10':'10101011','D01':'01010111','Xvx':'01111010','Xuy':'11001101','L':'10110111'}
assert sorted(pats)==sorted(list(table.values())+['01101011'[:0] or ''.join('1' if e in('ab','cd','ad','bc') else '0' for e in E)])
miss={k:{E[i] for i,c in enumerate(v) if c=='0'} for k,v in table.items()}
pairs=[p for p in itertools.combinations(table,2) if not miss[p[0]]&miss[p[1]]]
print('covering pairs:',pairs, len(pairs))
assert len(pairs)==7
inter=[p for p in itertools.combinations(table,2) if miss[p[0]]&miss[p[1]]]
Gc=nx.Graph(inter); assert not any(len(c)==3 for c in nx.enumerate_all_cliques(Gc)); print('intersection graph triangle-free')
# cube witnesses
Q=nx.Graph([(0,1),(0,3),(0,4),(1,2),(1,7),(2,3),(2,6),(3,5),(4,5),(4,7),(5,6),(6,7)])
def es(seq,closed):
    p=list(zip(seq,seq[1:]))+([(seq[-1],seq[0])] if closed else []); assert all(Q.has_edge(*e) for e in p); assert len(set(seq))==8; return {frozenset(e) for e in p}
W={'A0':es([0,1,2,3,5,6,7,4],1),'B':es([0,1,2,6,7,4,5,3],1),'C':es([0,3,2,1,7,6,5,4],1),'D':es([1,7,4,0,3,5,6,2],0),'E':es([0,4,7,1,2,6,5,3],0)}
p,q=frozenset((0,1)),frozenset((2,3))
assert p in W['A0'] and q in W['A0'] and p in W['B'] and q not in W['B'] and q in W['C'] and p not in W['C']
assert p not in W['D'] and q not in W['D'] and p not in W['E'] and q not in W['E']
print('cube supports:',{''.join(sorted(map(str,e))):''.join(k for k in W if e in W[k]) for e in map(frozenset,Q.edges())})
# root cube global H+- from six listed cycles
cyc=[[0,1,2,3,5,6,7,4],[0,1,2,6,7,4,5,3],[0,1,7,4,5,6,2,3],[0,1,7,6,2,3,5,4],[0,3,2,1,7,6,5,4],[0,3,5,6,2,1,7,4]]
C=[es(c,1) for c in cyc]
assert all(any(e in c and f not in c for c in C) for e in map(frozenset,Q.edges()) for f in map(frozenset,Q.edges()) if e!=f); print('cube global H+- from six cycles: OK')
# cap isomorphism (8.10)
F0={('U','a'),('a','b'),('b','V'),('X','c'),('c','d'),('d','Y'),('b','e'),('c','f'),('d','g'),('a','h'),('e','f'),('f','g'),('g','h'),('h','e')}
# F1(B): inner basic square e f g h with U-e,V-f,X-g,Y-h; outer a-b-c-d-a, X-c,Y-d, a-n_V=f, b-n_X=g; new U=n_U(e) side: U-e, V-h
F1={('a','b'),('b','c'),('c','d'),('d','a'),('X','c'),('Y','d'),('a','f'),('b','g'),('e','f'),('f','g'),('g','h'),('h','e'),('U','e'),('V','h')}
m={'a':'e','b':'h','c':'c','d':'d','e':'g','f':'b','g':'a','h':'f','U':'U','V':'V','X':'X','Y':'Y'}
img={frozenset((m[x],m[y])) for x,y in F0}; capok=img=={frozenset(e) for e in F1}
print('cap map is boundary-fixed isomorphism:', capok)
assert capok and all(m[v]==v for v in 'UVXY') and len(set(m.values()))==len(m)
print('LOCAL_CHECKS_PASSED')
