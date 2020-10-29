import xml.etree.ElementTree as ET
import os
from gurobipy import * 
import operator
from collections import defaultdict
import math
import itertools
import pylab as pyl

# SETS AND DATA USED BY THE ALL OF THE FORMULATIONS
# Number of days in the current dataset
Days = 5

# Minimum Periods Per Day (Not used)
MinPPD = -1
# Maximum Periods Per Day(Not used)
MaxPPD = -1
# Periods per day
PPD = -1

# Set of Rooms
Rooms = []
# Set of Courses
Courses = []

# Set of periods
P = []

# Data of days,  array of arrays, each index contains the time periods associated with that day
D = []

# Curricula and the courses they map to
Curricula = {}
# Not Used
CoursesToCurricula = {}
# Teachers and the courses they teach
T = defaultdict(list)

# Maps courses to number of lectures required for that course
L  = {} 
# Minimum number of working days required for each course
MND = {} 
# Demand for each course
DEM = {}
# Capacity of each room
CAP = {} 
# Unique room capacities
S = {}
# Name of the dataset
Name = ""
# For a course, the timeslots it cannot be in
PeriodConstraints = defaultdict(list)

# Parses the Dataset
# Pick one of ITC, Udine, Erlangen
# ITC specify the number 1-21
# Udine specifiy the number 1-7
# Erlangen provide the name of the dataset (find the data directory)
def ProcessData(dataset=1, ITC=True, Udine=False, Erlangen=False):    
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
    global CoursesToCurricula
    global T 
    global L  
    global MND  
    global DEM 
    global CAP  
    global S 
    global PeriodConstraints
    global Name
    
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
    if (ITC):
        relpath = f'CCT/ITC2007/comp{dataset:0>2}.xml'
        Name = f"comp{dataset:0>2}"
    elif (Udine):
        relpath = f'CCT/Udine/Udine{dataset}.xml'
        Name = f"Udine{dataset}"
    elif (Erlangen):
        relpath = f'CCT/Erlangen/{dataset}.xml'
        Name = f"{dataset}"


    absFilePath = os.path.join(script_dir, relpath)
    tree = ET.parse(absFilePath)
    root = tree.getroot()
    
    # Parses the data
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
                    
                    CoursesToCurricula[ggchild.attrib['ref']] = id

                Curricula[id] = CurriculaCourses

            elif tag == "constraint":
                course = attribs['course']
                if (attribs['type'] == 'period'):
                    for x in gchild:
                        day = int(x.attrib['day'])
                        period = int(x.attrib['period'])
                        PeriodConstraints[course].append(day + period*Days + 1)

    # List of periods
    P = range(1, Days * PPD + 1)
    # Partitions periods into their associated periods
    D = [[(j+1) + i*(Days) for i in range(PPD)] for j in range(Days)]


#---------------------- SINGLE OBJECTIVE FORMUALATIONS ---------------------------------------
# Single objective implementation not really used except for verification purposes
# Singular Room Planning Problem
def SimpleRoomPlanningProblem(output=True) :
    # need to implement the linear version of this algorithm maybe in the future
    # S = [math.ceil(DEM[c]/25)*25 for c in DEM]
    # S = [DEM[c] for c in Courses]
    # S = set(S).union({0})
    step = 25
    maxS = max([math.ceil(DEM[c] / step) * step for c in DEM])
    S = [i for i in range(step, maxS + 1, step)]
    # Calculate Helper Sets
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

# Actual Formulation Presented in the Paper
def SimpleTimeSlot(output=True):    
    # S = [DEM[c] for c in Courses]
    # Calculate Helper Sets    # S = set(S).union({0})
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

# Transcribed from the Author's Code
def ActualTimeSlotPresented(output=True):    
    # S = [DEM[c] for c in Courses]
    # Calculate Helper Sets    # S = set(S).union({0})
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


# ------------------------------ MULTI (BI) OBJECTIVE FORMULATIONS-------------------------------

# Room Planning vs Quality Replicates the Paper
# output : Output optimization
# UseStep : Use simplifying constraint 
# TimeLimit : Set timelimt for objectives
#  Pick One Of
#       Pareto Front
#           Pick Save Graph
#       Bounds
#       TuneGurobiExp (tunng gurobi experiments)
#       WarmStart(Does nothing)
def CCTModelForRoomPlanning(output=True, ParetoFront=False, Bounds=False, SaveGraph=False, 
        TimeLimit=GRB.INFINITY, UseStep=True, WarmStart=True, TuneGurobiExp = False) :
    global Courses;
    global S;

    if (UseStep) :
        step = 25 
        maxS = max([math.ceil(DEM[c] / step) * step for c in DEM])
        S = [i for i in range(step, maxS + 1, step)]
    
    
    # Calculate Helper Sets    
    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
    RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}

    m = Model()
    if (not output):
        m.setParam('OutputFlag', 0)

    m.setParam("BranchDir", -1)
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

    for x in RPlus:
        RPlus[x].BranchPriority = 10;

    # for x in Y:
    #     Y[x].BranchPriority = 15;
    
    AllLecPlanned = {c : m.addConstr(quicksum(Y[c, p] for p in P if (c, p) in Y) == L[c]) for c in Courses}
    #Additional Room Cuts
    AdditionalRoomCuts = m.addConstr(len(P) * quicksum(R[s] for s in S) >= sum(L[c] for c in Courses))

    CalculatedBranchVar = {(s) :  m.addConstr(quicksum(R[s] for s in SGT[s]) == RPlus[s]) for s in S}

    RoomAllocation = {(s) : m.addConstr(RPlus[s] * len(P) >= sum(L[c] for c in CGT[s])) for s in S}

    {pConstr : m.addConstr(quicksum(Y[c, p] for c in Courses if (c, p) in Y) <= quicksum(R[s] for s in S)) for p in P}
    {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y) <= RPlus[s]) for p in P for s in S}

    # Lecturers
    TeacherOneCourse = {(t, p): m.addConstr(quicksum(Y[c, p] for c in T[t] if (c, p) in Y) <= 1) for t in T for p in P}
    
    OneCurriculaPerPeriod = {(cu, p) : m.addConstr(quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)  <= 1) for cu in Curricula for p in P}
    
    DayInUseCalculator = {(c, d) : m.addConstr(quicksum(Y[c, p] for p in D[d] if (c, p) in Y) - DayInUse[c, d] >= 0) for c in Courses for d in range(len(D))}
    
    CalculationOfMNDViolation = {(c) : m.addConstr(quicksum(DayInUse[c, d] for d in range(len(D))) + DaysBelow[c] >= MND[c]) for c in Courses}

    CurriculaInPeriodCalculator = {(cu, p) : m.addConstr(CurrAssigned[cu, p] ==  quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)) for cu in Curricula for p in P}


    # the formulation to take for foreover to get decent results
    Compactness = {(cu, D[d][pi]) : m.addConstr(CurrAssigned[cu, D[d][pi]] -  quicksum(CurrAssigned[cu, pp] for pp in  TimeslotConseq(cu, pi, d)) <= CurrAlone[cu, D[d][pi]]) for cu in Curricula for d in range(len(D)) for pi in range(len(D[d]))}

    # for x in BVar:
    #     BVar[x].BranchPriority = 20

    ConseqTimePeriodsAlloc = {(p) : m.addConstr(Y[c, p]  - TimeslotUsed[p] <= 0)  for p in P for c in Courses if (c, p) in Y}
    ConsecPeriods = {p: m.addConstr((TimeslotUsed[p] - TimeslotUsed[p - 1]) <= 0) for p in P if p != 1}

    CalculateOverallMWD = m.addConstr(MinimumWorkingDaysBelow == quicksum(DaysBelow[c] for c in Courses))
    CalculateOverallCurrViol = m.addConstr(VarCurrCompactViol == quicksum(CurrAlone[x] for x in CurrAlone))
    
    fSeats = quicksum(s * R[s] for s in S)
    fQual = 5 * MinimumWorkingDaysBelow + 2 * VarCurrCompactViol
    
    m.setParam('TimeLimit', TimeLimit)

    #Didnt really work
    # The model itself has a really weak lp relaxation
    if (WarmStart):
        # m.setObjective(fSeats);
        # m.setParam('OutputFlag', False)
        # m.optimize();

        # m.addConstr(fSeats == m.objVal)  
        # m.setObjective(fQual);

        # WarmStartModel(m, [Y, R ,RPlus], [GRB.BINARY, GRB.INTEGER, GRB.INTEGER])
        return

    if (Bounds):
        OptimalObjectiveBoundsSolver(fSeats, fQual, m, 'Seats', 'Quality')
        return;
    
    if (ParetoFront):
        BiObjectiveSolver(fSeats, fQual, m, 25, "Seats", "Quality")
        return;

    if (TuneGurobiExp):
        TuneParameterExperimentation(fSeats, fQual, m)
        return
    
    m.setObjective(fQual)
    m.optimize()
    return m.objVal

# Time Periods Vs Quality Replicates the Paper
# output : Output optimization
# UseStep : Use simplifying constraint 
# TimeLimit : Set timelimt for objectives
#  Pick One Of
#       Pareto Front
#           Pick Save Graph
#       Bounds
def CCTModelForTimeSlots(output=True, ParetoFront=False, Bounds=False, TimeLimit=GRB.INFINITY, SaveGraph=False, TuneGurobiExp = False) :
    global PPD;
    global P;
    global D;
    global S;


    # Adds 5 teaching periods to each day
    PPD += 5
    P = range(1, Days * PPD + 1)
    D = [[(j+1) + i*(Days) for i in range(PPD)] for j in range(Days)]
    # Calculate Helper Sets
    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
    RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}
    AmtS = {s : sum(1 for r in Rooms if CAP[r] == s) for s in S }


    m = Model()
    if (not output):
        m.setParam('OutputFlag', 0)

    m.setParam("BranchDir", 1)
    m.setParam('Heuristics', 0.5)
    m.setParam('MIPFocus', 1)

    R = {(s) : m.addVar(vtype=GRB.INTEGER, ub=AmtS[s]) for s in S}
    RPlus = {(s) : m.addVar(vtype=GRB.INTEGER) for s in S}
    TimeslotUsed = {p : m.addVar(vtype = GRB.BINARY) for p in P}

    Y = {(c, p): m.addVar(vtype = GRB.BINARY) for c in Courses for p in P if p not in PeriodConstraints[c]}
    
    DayInUse = {(c, d) : m.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    DaysBelow = {(c) : m.addVar() for c in Courses}
    CurrAlone = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    CurrAssigned = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    
    MinimumWorkingDaysBelow = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)
    VarCurrCompactViol = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)

    RoomsAccomodate = {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y) <= len(RGT[s])) for s in S for p in P}
    
    AllLecPlanned = {c : m.addConstr(quicksum(Y[c, p] for p in P if (c, p) in Y) == L[c]) for c in Courses}
    AdditionalRoomCuts = m.addConstr(len(P) * quicksum(R[s] for s in S) >= sum(L[c] for c in Courses))

    CalculatedBranchVar = {(s) :  m.addConstr(quicksum(R[s] for s in SGT[s]) == RPlus[s]) for s in S}

    RoomAllocationConstr = {(s) : m.addConstr(RPlus[s] * len(P) >= sum(L[c] for c in CGT[s])) for s in S}
    EnsureEnoughTimeslotsAndRooms = {p : m.addConstr(quicksum(Y[c, p] for c in Courses if (c, p) in Y) <= quicksum(R[s] for s in S)) for p in P}
    EnsureRoomSizeIsBigEnough = {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y) <= RPlus[s]) for p in P for s in S}

    # Lecturers
    TeacherOneCourse = {(t, p): m.addConstr(quicksum(Y[c, p] for c in T[t] if (c, p) in Y) <= 1) for t in T for p in P}
    
    OneCurriculaPerPeriod = {(cu, p) : m.addConstr(quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)  <= 1) for cu in Curricula for p in P}
    
    DayInUseCalculator = {(c, d) : m.addConstr(quicksum(Y[c, p] for p in D[d] if (c, p) in Y) - DayInUse[c, d] >= 0) for c in Courses for d in range(len(D))}
    
    CalculationOfMNDViolation = {(c) : m.addConstr(quicksum(DayInUse[c, d] for d in range(len(D))) + DaysBelow[c] >= MND[c]) for c in Courses}

    CurriculaInPeriodCalculator = {(cu, p) : m.addConstr(CurrAssigned[cu, p] ==  quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)) for cu in Curricula for p in P}

    Compactness = {(cu, D[d][pi]) : m.addConstr(CurrAssigned[cu, D[d][pi]] - quicksum(CurrAssigned[cu, pp] for pp in  TimeslotConseq(cu, pi, d)) <= CurrAlone[cu, D[d][pi]]) for cu in Curricula for d in range(len(D)) for pi in range(len(D[d]))}

    ConseqTimePeriodsAlloc = {(p) : m.addConstr(Y[c, p]  - TimeslotUsed[p] <= 0)  for p in P for c in Courses if (c, p) in Y}
    ConsecPeriods = {p: m.addConstr((TimeslotUsed[p] - TimeslotUsed[p - 1]) <= 0) for p in P if p != 1}

    CalculateOverallMWD = m.addConstr(MinimumWorkingDaysBelow == quicksum(DaysBelow[c] for c in Courses))
    CalculateOverallCurrViol = m.addConstr(VarCurrCompactViol == quicksum(CurrAlone[x] for x in CurrAlone))
    
    # fSeats = quicksum(s * R[s] for s in S)
    fQual = 5 * MinimumWorkingDaysBelow + 2 * VarCurrCompactViol
    fTime = quicksum(TimeslotUsed[p] for p in P)
    for x in RPlus:
        RPlus[x].BranchPriority = 10

    
    m.setParam("TimeLimit", TimeLimit)

    if (ParetoFront) :
        BiObjectiveSolver(fTime, fQual, m, 1, "Time", "Quality", SaveGraph=SaveGraph)
        return

    if (Bounds):
        OptimalObjectiveBoundsSolver(fTime, fQual, m , 'Time', "Quality")
        return

    if (TuneGurobiExp):
        TuneParameterExperimentation(fTime, fQual, m)
        return

    m.setObjective(fTime)
    m.optimize()

# Time Periods vs Room Planning Replicates the Paper
def CCTModelForTimeSlotsAndRP(output=True, ParetoFront=False, Bounds=False, TimeLimit=GRB.INFINITY, 
        SaveGraph=False, UseStep=True, WarmStart=False, TuneGurobiExp=False):
    
    global PPD;
    global P;
    global D;
    global S;
    # Add 5 time slots to each day
    # re calculate relevant sets
    PPD += 5
    P = range(1, Days * PPD + 1)
    D = [[(j+1) + i*(Days) for i in range(PPD)] for j in range(Days)]
    if (UseStep) :
        step = 25 
        maxS = max([math.ceil(DEM[c] / step) * step for c in DEM])
        S = [i for i in range(step, maxS + 1, step)]

    # Calculate Helper Sets    
    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
    RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}


    m = Model()
    if (not output):
        m.setParam('OutputFlag', 0)


    m.setParam("BranchDir", 1)
    m.setParam('Heuristics', 0.5)
    m.setParam('MIPFocus', 1)

    AmtS = {s : sum(1 for r in Rooms if CAP[r] == s) for s in S }


    R = {(s) : m.addVar(vtype=GRB.INTEGER, ub=GRB.INFINITY) for s in S}
    RPlus = {(s) : m.addVar(vtype=GRB.INTEGER) for s in S}
    TimeslotUsed = {p : m.addVar(vtype = GRB.BINARY) for p in P}

    Y = {(c, p): m.addVar(vtype = GRB.BINARY) for c in Courses for p in P if p not in PeriodConstraints[c]}

    DayInUse = {(c, d) : m.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    DaysBelow = {(c) : m.addVar() for c in Courses}
    CurrAlone = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    CurrAssigned = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    
    MinimumWorkingDaysBelow = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)
    VarCurrCompactViol = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)

    # RoomsAccomodate = {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y) <= len(RGT[s])) for s in S for p in P}
    
    AllLecPlanned = {c : m.addConstr(quicksum(Y[c, p] for p in P if (c, p) in Y) == L[c]) for c in Courses}
    #Additional Room Cuts
    AdditionalRoomCuts = m.addConstr(len(P) * quicksum(R[s] for s in S) >= sum(L[c] for c in Courses))

    CalculatedBranchVar = {(s) :  m.addConstr(quicksum(R[s] for s in SGT[s]) == RPlus[s]) for s in S}

    RoomAllocationConstr = {(s) : m.addConstr(RPlus[s] * len(P) >= sum(L[c] for c in CGT[s])) for s in S}
    EnsureEnoughTimeslotsAndRooms = {p : m.addConstr(quicksum(Y[c, p] for c in Courses if (c, p) in Y) <= quicksum(R[s] for s in S)) for p in P}
    EnsureRoomSizeIsBigEnough = {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y) <= RPlus[s]) for p in P for s in S}

    # Lecturers
    TeacherOneCourse = {(t, p): m.addConstr(quicksum(Y[c, p] for c in T[t] if (c, p) in Y) <= 1) for t in T for p in P}
    
    OneCurriculaPerPeriod = {(cu, p) : m.addConstr(quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)  <= 1) for cu in Curricula for p in P}
    
    DayInUseCalculator = {(c, d) : m.addConstr(quicksum(Y[c, p] for p in D[d] if (c, p) in Y) - DayInUse[c, d] >= 0) for c in Courses for d in range(len(D))}
    
    CalculationOfMNDViolation = {(c) : m.addConstr(quicksum(DayInUse[c, d] for d in range(len(D))) + DaysBelow[c] >= MND[c]) for c in Courses}

    CurriculaInPeriodCalculator = {(cu, p) : m.addConstr(CurrAssigned[cu, p] ==  quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)) for cu in Curricula for p in P}

    Compactness = {(cu, D[d][pi]) : m.addConstr(CurrAssigned[cu, D[d][pi]] - quicksum(CurrAssigned[cu, pp] for pp in  TimeslotConseq(cu, pi, d)) <= CurrAlone[cu, D[d][pi]]) for cu in Curricula for d in range(len(D)) for pi in range(len(D[d]))}

    ConseqTimePeriodsAlloc = {(p) : m.addConstr(Y[c, p]  - TimeslotUsed[p] <= 0)  for p in P for c in Courses if (c, p) in Y}
    ConsecPeriods = {p: m.addConstr((TimeslotUsed[p] - TimeslotUsed[p - 1]) <= 0) for p in P if p != 1}

    CalculateOverallMWD = m.addConstr(MinimumWorkingDaysBelow == quicksum(DaysBelow[c] for c in Courses))
    CalculateOverallCurrViol = m.addConstr(VarCurrCompactViol == quicksum(CurrAlone[x] for x in CurrAlone))



    fSeats = quicksum(s * R[s] for s in S)
    fQual = 5 * MinimumWorkingDaysBelow + 2 * VarCurrCompactViol
    fTime = quicksum(TimeslotUsed[p] for p in P)

    m.setObjective(fTime)
    for x in R:
        R[x].BranchPriority = 10

    m.setParam('TimeLimit', TimeLimit)

    if (ParetoFront) :
        BiObjectiveSolver(fTime, fSeats, m, 1, "Periods", "Seats")
        return

    if (Bounds) :
        OptimalObjectiveBoundsSolver(fTime, fSeats, m , 1, "Periods", "Seats")
        return

    if (TuneGurobiExp):
        TuneParameterExperimentation(fTime, fSeats, m)
        return

# Returns consequctives for a particular time period in a particular day
def TimeslotConseq(cu, pi , d) :
    ret = []

    if (pi + 1) < PPD:
        ret.append(D[d][pi + 1])
    
    if (pi - 1) >= 0:
        ret.append(D[d][pi - 1])

    return ret


# --------------------------- METHODS TRIED -------------------------------------
#  Didnt really work
def WarmStartModel(m, varsR, setBackTo):
    print("-------------------")
    print("Warm Starting Model")

    oldOutputFlag = m.getParamInfo('OutputFlag')[2]
    
    m.setParam('OutputFlag', False)
    for X in varsR:
        for k in X:
            X[k].vtype = GRB.CONTINUOUS
    
    m.optimize()    
    print("Warm Start Yield: ", m.objVal)
    for i in range(len(setBackTo)):
        for X in varsR[i].values():
            X.vtype = setBackTo[i]

    m.setParam('OutputFlag', oldOutputFlag)
    print("Ending Warm Start")
    print("-------------------")

#  Gave up on doing this early as it 
#  Generates way too many arrangements to do
#  Didn't have time to do branch - cut/ dantzig - wolfe :(
def APrioriColGenQuality() :
    MultipleCourse = []
    for c in Courses:
        MultipleCourse += L[c] * [c] 
    
    roomArrangements = set()
    A = {}

    for i in itertools.permutations(MultipleCourse, PPD):
        roomArrangements.add(i)

        A[i] = {}

        for course in i:
            if course not in A[i]:
                A[i][course] = 1;
            else:
                A[i][course] += 1;

    
    # X = {p : m.addVar(vtype = GRB.BINARY) for p in roomArrangements}   
    m = Model()
    Y = {(p, d) : m.addVar(vtype=GRB.BINARY) for p in roomArrangements for d in range(Days)}

    OneArrangmentPerDay = {d : m.addConstr(quicksum(Y[p,d] for p in A) == 1) for d in range(Days)}

    for d in range(Days):
        for p in A:
            if (Y[p, d].x > 0.9):
                print("For day " + d + ":" + p)

#  Didnt really work
def PartitioningTheModel(PartitionsSize=10, modelType = 0):
    global Courses;
    global Curricula;

    retu = CCTModelForRoomPlanning(output=False, ParetoFront=False, WarmStart = False,  Bounds=False, TimeLimit=360)
    print(retu)

    if (PartitionsSize >= 0):
        Courses2 = Courses
        split = [Courses[i:i + PartitionsSize] for i in range(0, len(Courses), PartitionsSize)]
        Curricula2 = Curricula
        

        for i in split :
            Curricula = defaultdict(list)

            Courses = i;
            
            for k in i:
                Curricula[CoursesToCurricula[k]].append(k)

            Curricula3 = {}
            for k in Curricula:
                if (len(Curricula[k])) > 1 and (Curricula2[k] != 1) :
                    Curricula3[k] = Curricula[k]

            Curricula = Curricula3
            print(Curricula3)
            if (modelType == 0):
                retu = CCTModelForRoomPlanning(output=False, ParetoFront=False, WarmStart = False,  Bounds=False, TimeLimit=360)
                print(retu)
            elif (modelType == 1):
                pass
            else:
                pass

        Curricula = Curricula2
        Courses = Courses2



# ------------------------------- SOLUTION GENERATION METHODS  ----------------

# Solves the bounds of the pareto front
# Parameters:
# fx: LinExp of First Objective
# fy: LinExp of Second Objective
# m : Gurobi Model with all variables constraints added
# fxString : String Name of the First Objective
# fyString : String Name of the Second Objective
def OptimalObjectiveBoundsSolver(fx, fy, m, fxString, fyString):    
    print("Solving For Best " + fxString)
    m.setObjective(fx, GRB.MINIMIZE)
    m.optimize()
    
    bestFx = m.objVal
    bestFxT = m.Runtime;

    temp = m.addConstr(fx <= bestFx)


    m.setObjective(fy, GRB.MINIMIZE)
    m.optimize();

    fyWithBestFx = m.objVal
    fyWithBestFxT = m.Runtime;

    m.remove(temp)


    print("Solving For Best " + fyString)

    m.setObjective(fy, GRB.MINIMIZE)
    m.optimize()
    
    bestFy = m.objVal
    bestFyT = m.Runtime

    temp = m.addConstr(fy <= bestFy)


    m.setObjective(fx, GRB.MINIMIZE)
    m.optimize();

    fxWithBestFy = m.objVal
    fxWithBestFyT = m.Runtime

    m.remove(temp)
    # print("Original " + str(int(originalProblem)))
    print("Best ", fxString , bestFx, str(bestFxT), "| ", fyString , " with Best " , fxString, fyWithBestFx, str(fyWithBestFxT))
    print("Best " + fyString , bestFy, str(bestFyT),"|", fxString + " with Best " + fyString, fxWithBestFy , str(fxWithBestFyT))

# Epsilon Constraint method
# Parameters:
# fx: LinExp of First Objective
# fy: LinExp of Second Objective
# m : Gurobi Model with all variables constraints added
# fxName : String Name of the First Objective
# fyName : String Name of the Second Objective
# SaveGraph : True, automatically saves the graph, False automatically shows you the graph
def BiObjectiveSolver(fx, fy, m, delta, fxName, fyName, SaveGraph=True):
    m.setObjective(fy, GRB.MINIMIZE);
    m.optimize();
    minY = m.objVal
    TEMP = m.addConstr(fy == minY)    
    m.setParam('Method', 3)#Barrier Method
    m.setObjective(fx, GRB.MINIMIZE)
    m.optimize()

    #max x ← Minimize (f_x)
    maxX = m.objVal;
    m.remove(TEMP)

    # Minimize f(x) again without the constraint
    m.optimize()
    epsilon = m.objVal
    epsilonConstraint = m.addConstr(fx <= epsilon)
    m.setObjective(fy)

    print("MIN X", minY)
    print("MAX Y", maxX)
    print("epsilon or MIN Y" , epsilon)
    i = 0
    fx = []
    fyObj = []
    fyLB = []

    while True and i < 200:
        m.optimize()
        fxHat = epsilon
        fyHat = m.objVal
        fyHatLB = m.ObjBound

        
        if (len(fyObj) == 0 or  fyObj[-1] != fyHat):
            fx.append(fxHat)
            fyObj.append(fyHat)
            fyLB.append(fyHatLB)


        epsilon = epsilon + delta
        epsilonConstraint.RHS = epsilon

        i += 1;
        if (epsilon > maxX or fyHat <= minY):
            break;

    GraphParetoFront(fx, fyObj, fyLB, fxName, fyName, SaveGraph = SaveGraph)


# ---------------------------------- AUXILLIARY METHODS
# Calculates the bounds for each combination of the tuning parameters specified in the method
# Parameters:
# fx: LinExp of First Objective
# fy: LinExp of Second Objective
# m : Gurobi Model with all variables constraints added
# fxName : String Name of the First Objective
# fyName : String Name of the Second Objective
def TuneParameterExperimentation(fx, fy, m):
    Names = ['MIPFocus' , 'Heuristics', 'BranchDir']
    # The values you want to install, must be in the same order as the names
    # and must be valid values for each of the names
    Params = [(1, 3), (0.3, 0.5, 0.8), (-1, 0, 1)]
    
    # Set a timelimit
    m.setParam('TimeLimit', 420)
    # Turn off output
    m.setParam('OutputFlag', 0)
    # Gets the combination of all the values in the params
    for i in itertools.product(*Params):
        print("Optimizing : ")
        for p in range(len(i)):
            # Sets the parameters for e
            m.setParam(Names[p], i[p])
            print(Names[p], i[p])

        # Function prints out bounds and the time it took
        OptimalObjectiveBoundsSolver(fx, fy, m , "room", "quality")        
        # print("After A max of 10 mins", m.objVal, m.Runtime)
        # Don't want the previous runs to influence the current run
        m.reset()


# Graphs a particular pareto front generated from the bi - objective solver
# fx: LinExp of First Objective
# fyObj : The objective values at each step
# fyLB: The lower bound for each value in the same
# fxName : String Name of the First Objective
# fyName : String Name of the Second Objective
def GraphParetoFront(fx, fyObj, fyLB, fxName, fyName, SaveGraph=True):
    print(fx)
    print(fyObj)
    print(fyLB) 
    pyl.plot(fx, fyObj, 'rx-', markersize=10)
    pyl.plot(fx, fyLB, 'ro--',mfc='none', markersize=10)
    pyl.xlabel(fxName)
    pyl.ylabel(fyName)
    pyl.title(Name)

    if (SaveGraph):
        pyl.savefig(f"{Name}{fxName}{fyName}.png")
    else:
        pyl.show()

    pyl.close()


# CODE START HERE
if __name__ == "__main__":
    # Code can be run from here

    for i in range(1, 21): # Specify Datasets here
        
        print(i, end= " ")
        #Processes the current dataset
        ProcessData(i, ITC=True, Udine=False,Erlangen=False)
        # UNCOMMENT TO RUN
        # Teaching Periods vs Room Planning
        # CCTModelForTimeSlotsAndRP(output=True, ParetoFront=True, Bounds=False, SaveGraph=True, TimeLimit=420)
        
        # Teaching Periods vs Quality
        # CCTModelForTimeSlots(output=True, ParetoFront=True, Bounds=False, SaveGraph=True, TimeLimit=420)
        
        # Room Planning vs Quality
        # CCTModelForRoomPlanning(output=True, ParetoFront=True, Bounds=False, SaveGraph=True, TimeLimit=420)
        
        

        
