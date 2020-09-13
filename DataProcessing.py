import xml.etree.ElementTree as ET
import os
from gurobipy import * 

script_dir = os.path.dirname(__file__)
relpath = 'CCT/ITC2007/comp17.xml'
absFilePath = os.path.join(script_dir, relpath)

tree = ET.parse(absFilePath)
root = tree.getroot()


Days = 5
PPD = -1
Rooms = []
Courses = []
Curricula = []
P = []
D = []

CCOURSES = {}
L  = {} 
MND = {} 
DEM = {}
CAP = {} 
S = {0}

#DATA ANALYSIS
for child in root:
    if child.tag == "descriptor":
        for gchild in child:
            if (gchild.tag == "periods_per_day"):
                PPD = int(gchild.attrib["value"]) 
    
    for gchild in child:
        attribs = gchild.attrib
        tag = gchild.tag

        if tag == "course":
            id = attribs['id']
            Courses.append(id)
            MND[id] = int(attribs['min_days'])
            DEM[id] = int(attribs['students'])
            L[id] = int(attribs['lectures'])

        elif tag == "room":
            id = attribs['id']

            Rooms.append(id)
            CAP[id] = int(attribs['size'])
            S.add(int(attribs['size']))
        
        elif tag == "curriculum":
            id = attribs['id']

            Curricula.append(id)
            CurriculaCourses = []

            for ggchild in gchild:
                CurriculaCourses.append(ggchild.attrib['ref'])

            CCOURSES[id] = CurriculaCourses


P = range(1, Days * PPD + 1)
D = [[(j+1) + i*(PPD-1) for i in range(PPD)] for j in range(Days)]
SGT = {s : [ss for ss in S if ss >= s] for s in S}
CGT = {s : [c for c in Courses if DEM[c] >= s] for s in S}
RGT = {s : [c for c in Rooms if CAP[c] >= s] for s in S}


# need to implement the linear version of this algorithm maybe in the future
# RPP  = Model()
# R = {s : RPP.addVar(vtype= GRB.INTEGER) for s in S}
# RPP.setObjective(quicksum(s * R[s] for s in S))
# SetCovering = {s :  RPP.addConstr(len(P) * quicksum(R[ss] * ss for ss in SGT[s]) >= quicksum(L[c] for c in CGT[s])) for s in S}
# RPP.optimize();

# for s in S:
#     #if R[s].x != 0:
#     print("Allocated  " +  str(int(R[s].x)) + " of size " + str(int(s)))

# # Teaching Periods Problem
# TPP = Model()
# T = {p : TPP.addVar(vtype = GRB.BINARY) for p in P}
# TPP.setObjective(quicksum(T[p] for p in P), GRB.MINIMIZE)
# SetCovering = {s: TPP.addConstr(len(RGT[s])*quicksum(T[p] for p in P) >= quicksum(L[c] for c in CGT[s])) for s in S}
# ConsecPeriods = {p: TPP.addConstr((T[p] - T[p-1]) <= 0) for p in P if p != 1}
# TPP.optimize()



# for p in P:
#     if T[p].x > 0.9:
#         print("Allocated: " +  str(p))