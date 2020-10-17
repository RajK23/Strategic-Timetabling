import xml.etree.ElementTree as ET
import os
from gurobipy import * 
import operator
from collections import defaultdict
import math
import itertools

Days = 5
MinPPD = -1
MaxPPD = -1
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
    global MinPPD
    global MaxPPD
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
    global PeriodConstraints

    
    Days = 5
    MinPPD = -1
    MaxPPD = -1
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
                elif (gchild.tag == 'daily_lectures'):
                    MinPPD = int(gchild.attrib["min"])
                    MaxPPD = int(gchild.attrib["max"])
        
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
                CAP[id] = int(attribs['size'])
                S.add(int(attribs['size']))


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
                        day = int(x.attrib['day'])
                        period = int(x.attrib['period'])
                        PeriodConstraints[course].append(day + period*Days + 1)

    P = range(1, Days * PPD + 1)
    D = [[(j+1) + i*(Days) for i in range(PPD)] for j in range(Days)]

def SimpleRoomPlanningProblem(output=True) :
    # need to implement the linear version of this algorithm maybe in the future
    # S = [math.ceil(DEM[c]/25)*25 for c in DEM]
    # S = [DEM[c] for c in Courses]
    # S = set(S).union({0})
    step = 25
    maxS = max([math.ceil(DEM[c] / step) * step for c in DEM])
    S = [i for i in range(step, maxS + 1, step)]

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

UseHallsConditions=True
UseRoomHallConditions=True 
UseStageIandII=True 
TuneGurobi=False
NSols=1
ConstraintPenalty=None
AbortSolverOnZeroPenalty=False
UseRoomsAsTypes=False
SaveRelaxations=False
Seed=60
SaveBoundInformation=False
RoomBoundForEachTimeSlot=False
FocusOnlyOnBounds=False

def CCTModelForRoomPlanning(output=True, ParetoFront=False, Bounds=False, WarmStart=False, TuneGurobiExp = False, TimeLimit=GRB.INFINITY) :
    UseHallsConditions = True
    UseRoomHallConditions = True
    TuneGurobi = False

    step = 25
    maxS = max([math.ceil(DEM[c] / step) * step for c in DEM])
    S = [i for i in range(step, maxS + 1, step)]
    # S.remove(0)
    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
    RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}


    m = Model()
    if (not output):
        m.setParam('OutputFlag', 0)

    if (TuneGurobi) :
        m.setParam("BranchDir", 1)
        m.setParam('Heuristics', 0.3)
        m.setParam('MIPFocus', 1)


    R = {(s) : m.addVar(vtype=GRB.INTEGER) for s in S}
    RPlus = {(s) : m.addVar(vtype=GRB.INTEGER) for s in S}
    TimeslotUsed = {p : m.addVar(vtype = GRB.BINARY) for p in P}

    Y = {(c, p): m.addVar(vtype = GRB.BINARY) for c in Courses for p in P if p not in PeriodConstraints[c]}

    DayInUse = {(c, d) : m.addVar(ub=1) for c in Courses for d in range(len(D))}
    DaysBelow = {(c) : m.addVar() for c in Courses}
    CurrAlone = {(cu, p) : m.addVar(ub=1) for cu in Curricula for p in P}
    CurrAssigned = {(cu, p) : m.addVar(ub=1) for cu in Curricula for p in P}
    
    MinimumWorkingDaysBelow = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)
    VarCurrCompactViol = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)

    # RoomsAccomodate = {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y) <= len(RGT[s])) for s in S for p in P}
    
    AllLecPlanned = {c : m.addConstr(quicksum(Y[c, p] for p in P if (c, p) in Y) == L[c]) for c in Courses}
    #Additional Room Cuts
    m.addConstr(len(P) * quicksum(R[s] for s in S) >= sum(L[c] for c in Courses))

    {(s) :  m.addConstr(quicksum(R[s] for s in SGT[s]) == RPlus[s]) for s in S}

    if (UseRoomHallConditions) :
        {(s) : m.addConstr(RPlus[s] * len(P) >= sum(L[c] for c in CGT[s])) for s in S}

    if (UseHallsConditions):
        {p : m.addConstr(quicksum(Y[c, p] for c in Courses if (c, p) in Y) <= quicksum(R[s] for s in S)) for p in P}
        {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y) <= RPlus[s]) for p in P for s in S}

    # Lecturers
    TeacherOneCourse = {(t, p): m.addConstr(quicksum(Y[c, p] for c in T[t] if (c, p) in Y) <= 1) for t in T for p in P}
    
    {(cu, p) : m.addConstr(quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)  <= 1) for cu in Curricula for p in P}
    
    {(c, d) : m.addConstr(quicksum(Y[c, p] for p in D[d] if (c, p) in Y) - DayInUse[c, d] >= 0) for c in Courses for d in range(len(D))}
    
    CalculationOfMNDViolation = {(c) : m.addConstr(quicksum(DayInUse[c, d] for d in range(len(D))) + DaysBelow[c] >= MND[c]) for c in Courses}

    {(cu, p) : m.addConstr(CurrAssigned[cu, p] ==  quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)) for cu in Curricula for p in P}

    Compactness = {(cu, D[d][pi]) : m.addConstr(CurrAssigned[cu, D[d][pi]] - quicksum(CurrAssigned[cu, pp] for pp in  TimeslotConseq(cu, pi, d)) <= CurrAlone[cu, D[d][pi]]) for cu in Curricula for d in range(len(D)) for pi in range(len(D[d]))}

    {(p) : m.addConstr(Y[c, p]  - TimeslotUsed[p] <= 0)  for p in P for c in Courses if (c, p) in Y}
    ConsecPeriods = {p: m.addConstr((TimeslotUsed[p] - TimeslotUsed[p - 1]) <= 0) for p in P if p != 1}

    m.addConstr(MinimumWorkingDaysBelow == quicksum(DaysBelow[c] for c in Courses))
    m.addConstr(VarCurrCompactViol == quicksum(CurrAlone[x] for x in CurrAlone))
    
    fSeats = quicksum(s * R[s] for s in S)
    fQual = 5 * MinimumWorkingDaysBelow + 2 * VarCurrCompactViol
    
    m.setParam('TimeLimit', TimeLimit)
    # m.setObjective(fSeats)
    # m.optimize()
    for x in R:
        R[x].BranchPriority = 10

    if (WarmStart):
        print("Conducting warm start")
        # for s in S:
        #     R[s].vtype = GRB.CONTINUOUS
        #     RPlus[s].vtype = GRB.CONTINUOUS

        # for y in Y:
        #     Y[y].vtype = GRB.CONTINUOUS

        # m.setObjective(fQual)
        # m.optimize()

        # for s in S:
        #     R[s].vtype = GRB.INTEGER
        #     RPlus[s].vtype = GRB.INTEGER

        # for y in Y:
        #     Y[y].vtype = GRB.BINARY
        
        # m.optimize()

    if (Bounds):
        OptimalObjectiveBoundsSolver(fQual, fSeats, m, 'Quality', 'Seats')
        return;
    
    if (ParetoFront):
        BiObjectiveSolver()
        return;
    
    if (TuneGurobiExp):
        TuneParameterExperimentation(fSeats, fQual, m)
        return
# Currently Broken
def CCTModelForTimeSlots(output=True) :
    global PPD;
    global P;
    global D;

    PPD = PPD + 5 ;
    P = range(1, Days * PPD + 1)
    D = [[(j+1) + i*(Days) for i in range(PPD)] for j in range(Days)]

    UseHallsConditions = True
    UseRoomHallConditions = True
    TuneGurobi = True

    # step = 25
    # maxS = max([math.ceil(DEM[c] / step) * step for c in DEM])
    # S = [i for i in range(step, maxS + 1, step)]
    # S.remove(0)
    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
    RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}


    m = Model()
    if (not output):
        m.setParam('OutputFlag', 0)

    if (TuneGurobi) :
        m.setParam("BranchDir", 1)
        m.setParam('Heuristics', 0.5)
        m.setParam('MIPFocus', 1)


    R = {(s) : m.addVar(vtype=GRB.INTEGER, ub = GRB.INFINITY) for s in S}
    RPlus = {(s) : m.addVar(vtype=GRB.INTEGER) for s in S}
    TimeslotUsed = {p : m.addVar(vtype = GRB.BINARY) for p in P}

    Y = {(c, p): m.addVar(vtype = GRB.BINARY, ub=CalcUB(c, p)) for c in Courses for p in P}

    DayInUse = {(c, d) : m.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    DaysBelow = {(c) : m.addVar() for c in Courses}
    CurrAlone = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    CurrAssigned = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    
    MinimumWorkingDaysBelow = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)
    VarCurrCompactViol = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)

    # RoomsAccomodate = {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y) <= len(RGT[s])) for s in S for p in P}
    
    AllLecPlanned = {c : m.addConstr(quicksum(Y[c, p] for p in P if (c, p) in Y) == L[c]) for c in Courses}
    #Additional Room Cuts
    m.addConstr(len(P) * quicksum(R[s] for s in S) >= sum(L[c] for c in Courses))

    {(s) :  m.addConstr(quicksum(R[s] for s in SGT[s]) == RPlus[s]) for s in S}

    if (UseRoomHallConditions) :
        # {(s) : m.addConstr(RPlus[s] * len(P) >= sum(L[c] for c in CGT[s])) for s in S}
        pass
    if (UseHallsConditions):
        pass
        {p : m.addConstr(quicksum(Y[c, p] for c in Courses if (c, p) in Y) <= quicksum(R[s] for s in S)) for p in P}
        {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y) <= RPlus[s]) for p in P for s in S}

    # Lecturers
    TeacherOneCourse = {(t, p): m.addConstr(quicksum(Y[c, p] for c in T[t] if (c, p) in Y) <= 1) for t in T for p in P}
    
    {(cu, p) : m.addConstr(quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)  <= 1) for cu in Curricula for p in P}
    
    {(c, d) : m.addConstr(quicksum(Y[c, p] for p in D[d] if (c, p) in Y) - DayInUse[c, d] >= 0) for c in Courses for d in range(len(D))}
    
    CalculationOfMNDViolation = {(c) : m.addConstr(quicksum(DayInUse[c, d] for d in range(len(D))) + DaysBelow[c] >= MND[c]) for c in Courses}

    {(cu, p) : m.addConstr(CurrAssigned[cu, p] ==  quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)) for cu in Curricula for p in P}

    Compactness = {(cu, D[d][pi]) : m.addConstr(CurrAssigned[cu, D[d][pi]] - quicksum(CurrAssigned[cu, pp] for pp in  TimeslotConseq(cu, pi, d)) <= CurrAlone[cu, D[d][pi]]) for cu in Curricula for d in range(len(D)) for pi in range(len(D[d]))}

    {(p) : m.addConstr(Y[c, p]  - TimeslotUsed[p] <= 0)  for p in P for c in Courses if (c, p) in Y}
    ConsecPeriods = {p: m.addConstr((TimeslotUsed[p] - TimeslotUsed[p - 1]) <= 0) for p in P if p != 1}

    m.addConstr(MinimumWorkingDaysBelow == quicksum(DaysBelow[c] for c in Courses))
    m.addConstr(VarCurrCompactViol == quicksum(CurrAlone[x] for x in CurrAlone))
    
    # fSeats = quicksum(s * R[s] for s in S)
    fQual = 5 * MinimumWorkingDaysBelow + 2 * VarCurrCompactViol
    fTime = quicksum(TimeslotUsed[p] for p in P)
    # for x in RPlus:
    #     RPlus[x].BranchPriority = 10
    m.setObjective(fTime)
    m.optimize()

    if (m.Status == 2):
        print(m.objVal)      
    else:
        print("Infeasible or not solved to optimality")       # # COMP 4

    # y = m.objVal
    
    # m.addConstr(quicksum(s * R[s] for s in S) == y)

    # m.setObjective(5 * MinimumWorkingDaysBelow + 2 * VarCurrCompactViol)
    # m.optimize()

def CalcUB(c, p):
    if p in PeriodConstraints[c]:
        return 1
    else:
        return 0

def TimeslotConseq(cu, pi , d) :
    ret = []

    if (pi + 1) < PPD:
        ret.append(D[d][pi + 1])
    
    if (pi - 1) >= 0:
        ret.append(D[d][pi - 1])

    return ret


def OptimalObjectiveBoundsSolver(fx, fy, m, fxString, fyString):    
    print("Solving For Best " + fxString)
    m.setObjective(fx, GRB.MINIMIZE)
    m.optimize()
    
    bestFx = m.objVal
    temp = m.addConstr(fx == bestFx)


    m.setObjective(fy, GRB.MINIMIZE)
    m.optimize();

    fyWithBestFx = m.objVal
    m.remove(temp)


    print("Solving For Best " + fyString)

    m.setObjective(fy, GRB.MINIMIZE)
    m.optimize()
    
    bestFy = m.objVal
    temp = m.addConstr(fy == bestFy)


    m.setObjective(fx, GRB.MINIMIZE)
    m.optimize();

    fxWithBestFy = m.objVal
    m.remove(temp)
    # print("Original " + str(int(originalProblem)))
    print("Best " + fxString , bestFx, "| ", fyString + " with Best " + fxString, fyWithBestFx)
    print("Best " + fyString , bestFy, "|", fxString + " with Best " + fyString, fxWithBestFy)

def BiObjectiveSolver(fx, fy, m, delta) :
    m.setObjective(fy, GRB.MINIMIZE);
    m.optimize();
    m.setParam('OutputFlag', 0)
    minY = m.objVal

    TEMP = m.addConstr(fy == minY)    
    m.setObjective(fx, GRB.MINIMIZE)
    m.optimize()

    #max x ← Minimize (f_x)
    maxX = m.objVal;
    m.remove(TEMP)

    # Minimize f(x) again without the constraint
    m.optimize()
    epsilon = m.objVal
    epsilonConstraint = m.addConstr(fx == epsilon)
    m.setObjective(fy)

    print("MIN Y", minY)
    print("MAX X", maxX)
    print("EPSILON" , epsilon)
    while True:
        m.optimize()
        fx = epsilon
        fyHat = m.objVal
        epsilon = epsilon + delta
        epsilonConstraint.RHS = epsilon
        
        if (epsilon > maxX or fyHat <= minY):
            print("BREAKING")
            print(epsilon , "   ", maxX)
            print(fyHat , "   ", minY)
            break;

def TuneParameterExperimentation(fx, fy, m):
    Params = [(1, 3), (0.3, 0.5, 0.8)]
    Names = ['MIPFocus' , 'Heuristics']
    m.setParam('TimeLimit', 420)
    m.setParam('OutputFlag', 0)

    m.addConstr(fy == 35)
    m.setObjective(fx)
    for i in itertools.product(*Params):
        print("Optimizing : ")
        for p in range(len(i)):
            m.setParam(Names[p], i[p])
            print(Names[p], i[p])
        m.optimize()

        print("After A max of 10 mins", m.objVal, m.Runtime)
        m.reset()


if __name__ == "__main__":
    for i in [2]:
        print(i, end= " ")
        ProcessData(i)
        CCTModelForRoomPlanning(output=True, Bounds=True, WarmStart = False, TuneGurobiExp=False, TimeLimit=600)
        print()

# Reducing Room Approximation
# 