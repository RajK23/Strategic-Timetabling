import xml.etree.ElementTree as ET
import os
from gurobipy import * 
import operator
from collections import defaultdict

Days = 5
PPD = -1
Rooms = []
Courses = []
P = []
D = []
Curricula = {}
T = defaultdict(list)
L  = {} 
MND = {} 
DEM = {}
CAP = {} 
S = {0}

PeriodConstraints = defaultdict(list)
delta = 25
deltaMode = False


script_dir = os.path.dirname(__file__)
relpath = 'CCT/ITC2007/comp14.xml'
absFilePath = os.path.join(script_dir, relpath)

tree = ET.parse(absFilePath)
root = tree.getroot()
#DATA ANALYSIS
for child in root:
    if child.tag == "descriptor":
        for gchild in child:
            if (gchild.tag == "periods_per_day"):
                PPD = int(gchild.attrib["value"])
            elif (gchild.tag == 'days'):
                Days = int(gchild.attrib["value"])
    
    for gchild in child:
        attribs = gchild.attrib
        tag = gchild.tag

        if tag == "course":
            id = attribs['id']
            tc = attribs['teacher']

            Courses.append(id)
            MND[id] = int(attribs['min_days'])
            DEM[id] = int(attribs['students'])
            T[tc].append(id)
            L[id] = int(attribs['lectures'])

        elif tag == "room":
            id = attribs['id']

            Rooms.append(id)
            if not deltaMode:
                CAP[id] = int(attribs['size'])
                S.add(int(attribs['size']))
            else:
                CAP[id] = math.ceil(int(attribs['size']) / 25) * 25
                S.add(CAP[id])

        elif tag == "curriculum":
            id = attribs['id']

            CurriculaCourses = []
            for ggchild in gchild:
                CurriculaCourses.append(ggchild.attrib['ref'])

            Curricula[id] = CurriculaCourses

        elif tag == "constraint":
            course = attribs['course']
            if (attribs['type'] == 'period'):
                for x in gchild:
                    PeriodConstraints[course].append((int(x.attrib['day']), int(x.attrib['period'])))

# S = [math.ceil(DEM[c] / 25) * 25 for c in DEM]



def TimeslotConstraint(model, X):
    print(PeriodConstraints)
    for course in PeriodConstraints:
        model.addConstr( quicksum(X[course, day*Days + period + 1]  for day, period in PeriodConstraints[course]) == 1)


print(list(S))
Obj = {}
for i in range(len(list(S)[:-1])): 
    for c in Courses:
        for p in P:
            Obj[(S[i], c, p)] = min(DEM[c] - S[i], S[i + 1] - S[i]) 
# print(Obj)
P = range(1, Days * PPD + 1)
D = [[(j+1) + i*(Days) for i in range(PPD)] for j in range(Days)]
SGT = {s : [ss for ss in S if ss > s] for s in S}
CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}



m = Model()
X = {(c, p) : m.addVar(vtype = GRB.BINARY) for c in Courses for p in P}
Y = {(s, c, p) : m.addVar(vtype = GRB.BINARY) for s in S for c in Courses for p in P}
W = {(c) : m.addVar() for c in Courses}
Z = {(c,d): m.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
Q = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
V = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}

TimeslotConstraint(m, X)
# Implement Room Constraint 

#b
# RoomsAccomodate = {(s, p) : m.addConstr(quicksum(X[c, p] for c in CGT[s]) <= len(RGT[s])) for s in S for p in P}
#c
AllLecPlanned = {c : m.addConstr(quicksum(X[c, p] for p in P) == L[c]) for c in Courses}
{p : m.addConstr(quicksum(X[c, p] for c in Courses) <= len(Rooms)) for p in P}
{(s,c,p) : m.addConstr(X[c, p] - Y[s, c, p] >= 0) for s in S for c in CGT[s] for p in P}
{(s, p) : m.addConstr(quicksum(X[c, p] - Y[s, c, p] for c in CGT[s])  <= len(RGT[s])) for s in S for p in P}
{(c, d) : m.addConstr(quicksum(X[c, p] for p in D[d]) - Z[c, d] >= 0) for c in Courses for d in range(len(D))}
{(c) : m.addConstr(quicksum(Z[c, d] for d in range(len(D))) + W[c] >= MND[c]) for c in Courses}
{(cu, p) : m.addConstr(quicksum(X[c, p] for c in Curricula[cu]) - Q[cu, p]  == 0) for cu in Curricula for p in P}

g1 = {(cu,p):  m.addConstr(-Q[cu, p-1] + Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1) in Q and (cu, p+1) in Q}
g2 = {(cu,p):  m.addConstr(Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1)  not in Q}
g3 = {(cu,p):  m.addConstr(-Q[cu, p-1] + Q[cu,p]  - V[cu, p] <= 0) for (cu, p) in Q if (cu, p+1) not in Q}

{(t, p): m.addConstr(quicksum(X[c, p] for c in T[t]) <= 1) for t in T for p in P}

# Must set z to one to be able to plannign a lecture in this class

# #Compactness
# quicksum(Obj[s, c, p] * )
# #Teacher can
m.setObjective(quicksum(5 * W[c] for c in Courses) + quicksum(2 * V[cu,p] for cu in Curricula for p in P)+ quicksum(Y[s, c, p] * Obj[s, c, p] for (s, c, p) in Obj))

m.optimize()    