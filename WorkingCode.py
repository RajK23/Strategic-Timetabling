import xml.etree.ElementTree as ET
import os
from gurobipy import * 
import operator
from collections import defaultdict
import math




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


def ProcessData(dataset=1):    
    global P
    global D 
    global Days 
    global PPD 
    global Rooms 
    global Courses 
    global P 
    global D 
    global Curricula 
    global T 
    global L  
    global MND  
    global DEM 
    global CAP  
    global S 



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
    relpath = f'CCT/ITC2007/comp{dataset:0>2}.xml'
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
    P = range(1, Days * PPD + 1)
    D = [[(j+1) + i*(Days) for i in range(PPD)] for j in range(Days)]
    # SGT = {s : [ss for ss in S if ss >= s] for s in S}
    # CGT = {s : [c for c in Courses if DEM[c] >= s] for s in S}
    # RGT = {s : [r for r in Rooms if CAP[r] >= s] for s in S}

def SimpleRoomPlanningProblem(output=True) :
    # need to implement the linear version of this algorithm maybe in the future
    # S = [math.ceil(DEM[c]/25)*25 for c in DEM]
    S = [DEM[c] for c in Courses]
    S = set(S).union({0})
    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}

    RPP  = Model()
    if (output == False) :
        RPP.setParam('OutputFlag', 0)
    

    R = {s : RPP.addVar(vtype= GRB.INTEGER) for s in S}
    RPP.setObjective(quicksum(s * R[s] for s in S))
    SetCovering = {s :  RPP.addConstr(len(P) * quicksum(R[ss]  for ss in SGT[s]) >= quicksum(L[c] for c in CGT[s])) for s in S}
    
    RPP.optimize();
    print(RPP.objVal)#, sorted([(s, round(R[s].x))] for s in S )[1:])

def SimpleTimeSlot(output=True):    
    # S = [DEM[c] for c in Courses]
    # S = set(S).union({0})
    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
    RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}

    m = Model()
    TP = {p : m.addVar(vtype = GRB.BINARY) for p in P}
    
    if (output == False) :
        m.setParam('OutputFlag', 0)


    m.setObjective(quicksum(TP[p] for p in P), GRB.MINIMIZE)
    
    ConsecPeriods = {p: m.addConstr((TP[p+1] - TP[p]) <= 0) for p in P if p+1 in P}
    m.addConstr(len(Rooms) * quicksum(TP[p] for p in TP) >= sum(L[c] for c in L))

    for s in S:
        m.addConstr(len(RGT[s]) * len(P) >= sum(L[c] for c in CGT[s]))

    m.optimize();

    if (m.Status == 2):
        print(m.objVal)      
    else:
        print("Infeasible or not solved to optimality")   

def ActualTimeSlotPresented(output=True):    
    # S = [DEM[c] for c in Courses]
    # S = set(S).union({0})
    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
    RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}

    m = Model()
    TP = {p : m.addVar(vtype = GRB.BINARY) for p in P}


    if (output == False) :
        m.setParam('OutputFlag', 0)


    m.setObjective(quicksum(TP[p] for p in P), GRB.MINIMIZE)
    SetCovering = {s: m.addConstr(len(RGT[s])*quicksum(TP[p] for p in P) >= quicksum(L[c] for c in CGT[s])) for s in S}
    ConsecPeriods = {p: m.addConstr((TP[p+1] - TP[p]) <= 0) for p in P if p+1 in P}
    m.optimize();

    if (m.Status == 2):
        print(m.objVal)      
    else:
        print("Infeasible or not solved to optimality")


if __name__ == "__main__":
    # ProcessData(2)
    for i in range(1, 22):
        ProcessData(i)
        print(i, end= " ")
        ActualTimeSlotPresented(output=False)
        # SimpleTimeSlot(output=False)
