''''
This File Just Serves as Place where I dump code just in case I need a backup.
'''

import xml.etree.ElementTree as ET
import os
from gurobipy import * 
import operator
from collections import defaultdict

script_dir = os.path.dirname(__file__)
relpath = 'CCT/ITC2007/comp13.xml'
absFilePath = os.path.join(script_dir, relpath)

tree = ET.parse(absFilePath)
root = tree.getroot()


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

            dbl = attribs['double_lectures']

            Courses.append(id)
            MND[id] = int(attribs['min_days'])
            DEM[id] = int(attribs['students'])
            T[tc].append(id)
            
            # if dbl == "yes":
            #     L[id] = int(attribs['lectures']) * 1.2
            # else:
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
PPD = PPD + 5
P = range(1, Days * PPD + 1)
D = [[(j+1) + i*(Days) for i in range(PPD)] for j in range(Days)]
SGT = {s : [ss for ss in S if ss > s] for s in S}
CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}

# print("Output of Processed Data")
# print(Days)
# print(PPD)
# print(Rooms)
# print(Courses)
# print(len(P))
# print(D)
# print(L)
# print()
# print(MND)
# print()
# print(DEM)
# print()
# print(CAP)
# print(S)

# print(SGT)
# print()
# print(CGT)
# print()


def RoomPlanningAndQual():
    
    # 4b
    S = [math.ceil(DEM[c]/25)*25 for c in DEM]
    S = set(S).union({0})
    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
    RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}

    m  = Model()
    # Room Planning Variables
    R = {s : m.addVar(vtype= GRB.INTEGER) for s in S}
    # Quality Variables
    X = {(c, p) : m.addVar(vtype = GRB.BINARY) for c in Courses for p in P}
    W = {(c) : m.addVar() for c in Courses}
    Z = {(c, d) : m.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    Q = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    V = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}


    # 4a
    fSeats  = quicksum(s * R[s] for s in S)
    # 4b
    fQual = quicksum(5 * W[c] for c in Courses) + quicksum(2 * V[cu,p] for cu in Curricula for p in P)
    

    #4c 4d
    Connecting = {(s, p): m.addConstr(quicksum(X[c, p] for c in CGT[s])  <= quicksum(R[ss] for ss in SGT[s])) for s in S for p in P}

    # 4e
    # 1c
    AllLecPlanned = {c : m.addConstr(quicksum(X[c, p] for p in P) == L[c]) for c in Courses}

    # 1d
    SetOnDayVariable = {(c, d) : m.addConstr(quicksum(X[c, p] for p in D[d]) - Z[c, d] >= 0) for c in Courses for d in range(len(D))}

    # 1e
    MinimumDayPenalty  = {(c) : m.addConstr(quicksum(Z[c, d] for d in range(len(D))) + W[c] >= MND[c]) for c in Courses}

    # 1f
    CurriculaPressure = {(cu, p) : m.addConstr(quicksum(X[c, p] for c in Curricula[cu]) - Q[cu, p]  == 0) for cu in Curricula for p in P}

    # 1g Calculates Compactness
    g1 = {(cu,p):  m.addConstr(-Q[cu, p-1] + Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1) in Q and (cu, p+1) in Q}
    g2 = {(cu,p):  m.addConstr(Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1)  not in Q}
    g3= {(cu,p):  m.addConstr(-Q[cu, p-1] + Q[cu,p]  - V[cu, p] <= 0) for (cu, p) in Q if (cu, p+1) not in Q}

    # 1h
    TeacherPeriods = {(t, p): m.addConstr(quicksum(X[c, p] for c in T[t]) <= 1) for t in T for p in P}

    # 4f
    SetCovering = {s :  m.addConstr(len(P) * quicksum(R[ss]  for ss in SGT[s]) >= quicksum(L[c] for c in CGT[s])) for s in S}

    m.setObjective(fSeats)
    m.optimize()
    
    m.setObjective(fQual)
    m.optimize()


def LinearTimeRoomPlanningProblem():
    RoomSizes = [math.ceil(DEM[c]/25)*25 for c in DEM]
    RoomSizes = sorted(set(RoomSizes), reverse=True)
    print(RoomSizes)
    B = sorted([(DEM[c], L[c], c) for c in L], reverse=True)
    print(len(P))
    i = 0;
    remainingTimeslots = 0;
    
    rooms = []

    for dem, lectures, c in B:
        # print(dem, lectures, c)

        if (remainingTimeslots >= lectures):
            remainingTimeslots -= lectures
        else:
            remaining =  lectures - remainingTimeslots

            while (i < len(RoomSizes) - 1 and dem <= RoomSizes[i + 1]):
                i = i+1;    
                
            remainingTimeslots = len(P) - remaining 
            rooms.append(RoomSizes[i])
    print(rooms)
    print(sum(rooms))
    # 150 = (130, 6, 'c0001'), (117, 7, 'c0004'), (75, 6, 'c0002'), (75, 3, 'c0005'), (65, 8, 'c0015')
    # 75 = (65, 7, 'c0016'), (65, 2, 'c0017'), (65, 1, 'c0014'), (55, 8, 'c0025'), (55, 5, 'c0078'), (55, 4, 'c0024'), (31, 3, 'c0033')
    # 50 = (31, 3, 'c0033') (31, 1, 'c0032'), (20, 5, 'c0030'), (14, 6, 'c0066'), (11, 5, 'c0031'),(10, 6, 'c0071'), (10, 4, 'c0062')
    # 25 = (10, 1, 'c0062') (9, 6, 'c0072'), (9, 6, 'c0068'), (8, 6, 'c0063'), (7, 6, 'c0069'), (7, 6, 'c0059')
# print(Curricula)

def QualityProblem() :
    Quality = Model()
    X = {(c, p) : Quality.addVar(vtype = GRB.BINARY) for c in Courses for p in P}
    W = {(c) : Quality.addVar() for c in Courses}
    Z = {(c, d) : Quality.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    Q = {(cu, p) : Quality.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    V = {(cu, p) : Quality.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}

    Quality.setObjective(quicksum(5 * W[c] for c in Courses) + quicksum(2 * V[cu,p] for cu in Curricula for p in P));

    TimeslotConstraint(Quality, X)
    #b
    RoomsAccomodate = {(s, p) : Quality.addConstr(quicksum(X[c, p] for c in CGT[s]) <= len(RGT[s])) for s in S for p in P}
    #c
    AllLecPlanned = {c : Quality.addConstr(quicksum(X[c, p] for p in P) == L[c]) for c in Courses}
    #Must set z to one to be able to plannign a lecture in this class
    d = {(c, d) : Quality.addConstr(quicksum(X[c, p] for p in D[d]) - Z[c, d] >= 0) for c in Courses for d in range(len(D))}
    e = {(c) : Quality.addConstr(quicksum(Z[c, d] for d in range(len(D))) + W[c] >= MND[c]) for c in Courses}
    f = {(cu, p) : Quality.addConstr(quicksum(X[c, p] for c in Curricula[cu]) - Q[cu, p]  == 0) for cu in Curricula for p in P}
    
    #Compactness
    g1 = {(cu,p):  Quality.addConstr(-Q[cu, p-1] + Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1) in Q and (cu, p+1) in Q}
    g2 = {(cu,p):  Quality.addConstr(Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1)  not in Q}
    g3 = {(cu,p):  Quality.addConstr(-Q[cu, p-1] + Q[cu,p]  - V[cu, p] <= 0) for (cu, p) in Q if (cu, p+1) not in Q}
    #Teacher can
    h = {(t, p): Quality.addConstr(quicksum(X[c, p] for c in T[t]) <= 1) for t in T for p in P}



    Quality.optimize();

    for (c, d) in Z:
        if Z[c, d].x > 0.9:
            print(f'Course {c} is on day {d}')

    # for (c, p) in X:
    #     if X[c, p].x > 0.9:
    #         print(f'Course {c} has been scheduled in Timeslot {p}')

def TimeslotConstraint(model, X):
    # print(PeriodConstraints)
    print(D)
    for course in PeriodConstraints:
        print(course)
        # model.addConstr( quicksum(X[course, day*Days + period + 1]  for day, period in PeriodConstraints[course]) ==1)
        for day, period in PeriodConstraints[course]:
            con = day + period * Days + 1
            print(day, period, con)
            model.addConstr(X[course, con] == 0) 


def RoomPlanningProblem() :
    # need to implement the linear version of this algorithm maybe in the future
    RPP  = Model()
    R = {s : RPP.addVar(vtype= GRB.INTEGER) for s in S}
    RPP.setObjective(quicksum(s * R[s] for s in S))
    SetCovering = {s :  RPP.addConstr(len(P) * quicksum(R[ss]  for ss in SGT[s]) >= quicksum(L[c] for c in CGT[s])) for s in S}
    
    A = {s : RPP.addConstr(R[s] == sum(1 for r in Rooms if CAP[r] == s)) for s in S}
    RPP.optimize();

    for s in S:
        print(s, R[s].x)

#Can be solved if the Dem is pre sorted
def TeachingPeriodsProblem() :
    # Teaching Periods Problem
    TPP = Model()
    TP = {p : TPP.addVar(vtype = GRB.BINARY) for p in P}
    print(sum(L[c] for c in L))

    TPP.setObjective(quicksum(TP[p] for p in P), GRB.MINIMIZE)

    # SetCovering = {s: TPP.addConstr(len(RGT[s])*quicksum(TP[p] for p in P) 
    #         >= quicksum(L[c] for c in CGT[s])) for s in S}
    
    TPP.addConstr(len(Rooms) * quicksum(TP[p] for p in TP) >= sum(L[c] for c in L))

    ConsecPeriods = {p: TPP.addConstr((TP[p+1] - TP[p]) <= 0) for p in P if p+1 in P}
    TPP.optimize();

    print(sum(L[c] for c in L))

def RoomsPlanningvsQuality():
    # S = [math.ceil(DEM[c] / 25) * 25 for c in DEM]
    # S = set(S).union({0})

    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] >= s] for s in S}
    RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}

    RPQ = Model() 

    X = {(c, p) : RPQ.addVar(vtype = GRB.BINARY) for c in Courses for p in P}
    Z = {(c, d) : RPQ.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    W = {(c) : RPQ.addVar(vtype = GRB.INTEGER) for c in Courses}
    Q = {(cu, p) : RPQ.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    V = {(cu, p) : RPQ.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    R = {s : RPQ.addVar(vtype= GRB.INTEGER) for s in S}
    #R_PLUS = {(s) : RPQ.addVar(vtype = GRB.BINARY) for s in S for p in P}

    fQual = quicksum(5 * W[c] for c in Courses) + quicksum(2 * V[cu,p] for cu in Curricula for p in P)
    fSeats = quicksum(s * R[s] for s in S)

    TimeslotConstraint(RPQ, X)
    c4 = {(s, p): RPQ.addConstr(quicksum(X[c, p] for c in CGT[s])  <= quicksum(R[ss] for ss in SGT[s])) for s in S for p in P}
    #4e
    AllLecPlanned = {c : RPQ.addConstr(quicksum(X[c, p] for p in P) == L[c]) for c in Courses}
    #Must set z to one to be able to plannign a lecture in this class
    d = {(c, d) : RPQ.addConstr(quicksum(X[c, p] for p in D[d]) - Z[c, d] >= 0) for c in Courses for d in range(len(D))}
    e = {(c) : RPQ.addConstr(quicksum(Z[c, d] for d in range(len(D))) + W[c] >= MND[c]) for c in Courses}
    f = {(cu, p) : RPQ.addConstr(quicksum(X[c, p] for c in Curricula[cu]) - Q[cu, p]  == 0) for cu in Curricula for p in P}
    
    #Compactness
    compactness = {(cu, D[d][pi]) : RPQ.addConstr(Q[cu, D[d][pi]] - 
            quicksum(Q[cu, p] for p in TimeslotConseq(cu, pi, d)) <=  
            V[cu, D[d][pi]]) for cu in Curricula for d in range(len(D)) for pi in range(len(D[d]))}
    #Teacher can
    h = {(t, p): RPQ.addConstr(quicksum(X[c, p] for c in T[t]) <= 1) for t in T for p in P}

    #4f
    SetCovering = {s :  RPQ.addConstr(len(P) * quicksum(R[ss] * ss for ss in SGT[s]) >= quicksum(L[c] for c in CGT[s])) for s in S}

    # BoundsSolver(fSeats, fQual, RPQ, "Room Planning", "Quality")
    # BiObjectiveSolver( fQual, fSeats, RPQ, 1)
    RPQ.setObjective(fQual);

    RPQ.optimize()
    print(sum(W[c].x for c in Courses))
    print(sum(V[cu,p].x for cu in Curricula for p in P))


def TimeslotConseq(cu, pi , d) :
    # print(cu, pi, d, D[d][pi])
    ret = []

    if (pi + 1) < PPD:
        ret.append(D[d][pi + 1])
    
    if (pi - 1) >= 0:
        ret.append(D[d][pi - 1])

    return ret
def RoomsPlanningvsTeachingPeriods():
    RPTP = Model ()

    X = {(c, p) : RPTP.addVar(vtype = GRB.BINARY) for c in Courses for p in P}
    Z = {(c, d) : RPTP.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    W = {(c) : RPTP.addVar(vtype = GRB.INTEGER) for c in Courses}
    Q = {(cu, p) : RPTP.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    V = {(cu, p) : RPTP.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    R = {s : RPTP.addVar(vtype= GRB.INTEGER) for s in S}
    TP = {p : RPTP.addVar(vtype = GRB.BINARY) for p in P}


    fTime = quicksum(TP[p] for p in P)
    fSeats = quicksum(s * R[s] for s in S)

    AllLecPlanned = {c : RPTP.addConstr(quicksum(X[c, p] for p in P) == L[c]) for c in Courses}
    #Must set z to one to be able to plannign a lecture in this class
    d = {(c, d) : RPTP.addConstr(quicksum(X[c, p] for p in D[d]) - Z[c, d] >= 0) for c in Courses for d in range(len(D))}
    e = {(c) : RPTP.addConstr(quicksum(Z[c, d] for d in range(len(D))) + W[c] >= MND[c]) for c in Courses}
    f = {(cu, p) : RPTP.addConstr(quicksum(X[c, p] for c in Curricula[cu]) - Q[cu, p]  == 0) for cu in Curricula for p in P}
    
    #Compactness
    g1 = {(cu,p):  RPTP.addConstr(-Q[cu, p-1] + Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1) in Q and (cu, p+1) in Q}
    g2 = {(cu,p):  RPTP.addConstr(Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1)  not in Q}
    g3 = {(cu,p):  RPTP.addConstr(-Q[cu, p-1] + Q[cu,p]  - V[cu, p] <= 0) for (cu, p) in Q if (cu, p+1) not in Q}
    #Teacher can
    h = {(t, p): RPTP.addConstr(quicksum(X[c, p] for c in T[t]) <= 1) for t in T for p in P}

    SetCoveringR = {s :  RPTP.addConstr(len(P) * quicksum(R[ss] * ss for ss in SGT[s]) >= quicksum(L[c] for c in CGT[s])) for s in S}

    SetCoveringTP = {s: RPTP.addConstr(len(RGT[s])*quicksum(TP[p] for p in P) >= quicksum(L[c] for c in CGT[s])) for s in S}
    ConsecPeriods = {p: RPTP.addConstr((TP[p] - TP[p-1]) <= 0) for p in P if p != 1}

    c4 = {(s, p): RPTP.addConstr(quicksum(X[c, p] for c in CGT[s])  <= quicksum(R[ss] for ss in SGT[s])) for s in S for p in P}
    c5 = {(c, p): RPTP.addConstr(X[c, p] - TP[p] <= 0) for c in Courses for p in P}

    BiObjectiveSolver(fTime, fSeats, RPTP, 1);

def TeachingPeriodsvsQuality():
    TPQ = Model()
    #Teaching Period Variable
    TP = {p : TPQ.addVar(vtype = GRB.BINARY) for p in P}
    # Quality Variables
    X = {(c, p) : TPQ.addVar(vtype = GRB.BINARY) for c in Courses for p in P}
    Z = {(c, d) : TPQ.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    W = {(c) : TPQ.addVar(vtype = GRB.INTEGER) for c in Courses}
    Q = {(cu, p) : TPQ.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    V = {(cu, p) : TPQ.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}

    TimeslotConstraint(TPQ, X)
    fTime = quicksum(TP[p] for p in P)
    fQual = quicksum(5 * W[c] for c in Courses) + quicksum(2 * V[cu,p] for cu in Curricula for p in P)

    ConnectingConstraint = {(c, p): TPQ.addConstr(X[c, p] - TP[p] <= 0) for c in Courses for p in P}
    # Quality Problem
    RoomsAccomodate = {(s, p) : TPQ.addConstr(quicksum(X[c, p] for c in CGT[s]) <= len(RGT[s])) for s in S for p in P}
    AllLecPlanned = {c : TPQ.addConstr(quicksum(X[c, p] for p in P) == L[c]) for c in Courses}
    #Must set z to one to be able to plannign a lecture in this class
    
    ObjectivePressure = {(c, d) : TPQ.addConstr(quicksum(X[c, p] for p in D[d]) - Z[c, d] >= 0) for c in Courses for d in range(len(D))}
    CalculationOfMNDViolation = {(c) : TPQ.addConstr(quicksum(Z[c, d] for d in range(len(D))) + W[c] >= MND[c]) for c in Courses}
    CalculationOfCurriculumCompactness = {(cu, p) : TPQ.addConstr(quicksum(X[c, p] for c in Curricula[cu]) - Q[cu, p]  == 0) for cu in Curricula for p in P}
    
    #Curriculum Compactness
    g1 = {(cu,p):  TPQ.addConstr(-Q[cu, p-1] + Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1) in Q and (cu, p+1) in Q}
        #Boundary Constraints
    g2 = {(cu,p):  TPQ.addConstr(Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1)  not in Q}
    g3 = {(cu,p):  TPQ.addConstr(-Q[cu, p-1] + Q[cu,p]  - V[cu, p] <= 0) for (cu, p) in Q if (cu, p+1) not in Q}
    
    #Teacher can only teaching in one course at one time
    TeacherOneCourse = {(t, p): TPQ.addConstr(quicksum(X[c, p] for c in T[t]) <= 1) for t in T for p in P}

    #Teaching Periods
    # SetCovering = {s: TPQ.addConstr(len(RGT[s])*quicksum(TP[p] for p in P) >= quicksum(L[c] for c in CGT[s])) for s in S}
    ConsecPeriods = {p: TPQ.addConstr((TP[p] - TP[p-1]) <= 0) for p in P if p != 1}
    TPQ.setObjective(quicksum(TP[p] for p in TP))
    TPQ.optimize()
    # TPQ.setParam('OutputFlag', 0)
    # Look at adding in class constraints.
    # BiObjectiveSolver(fTime, fQual, TPQ, 1)
    #BoundsSolver(fTime, fQual, TPQ, "TimePeriods", "Quality")

def BiObjectiveSolver(fx, fy, m, delta) :
    m.setObjective(fy, GRB.MINIMIZE);
    m.optimize();
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

def BoundsSolver(fx, fy, m, fxString, fyString):    

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


def CCTModel(output=True) :
    UseHallsConditions = True
    UseRoomHallConditions = True
    TuneGurobi = True

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


    if (TuneGurobi):
        m.setParam('Heuristics', 0.5)
        m.setParam('MIPFocus', 1)

    Y = {(c, p): m.addVar(vtype = GRB.BINARY) for c in Courses for p in P if p not in PeriodConstraints[c]}

    
    R = {(s) : m.addVar(vtype=GRB.INTEGER) for s in S}
    RPlus = {(s) : m.addVar(vtype=GRB.INTEGER) for s in S}
    DayInUse = {(c, d) : m.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    DaysBelow = {(c) : m.addVar() for c in Courses}
    TimeslotUsed = {p : m.addVar(vtype = GRB.BINARY) for p in P}
    RoomStability = m.addVar(lb=0,ub=len(Courses)*(len(Rooms) - 1), obj=0, vtype=GRB.CONTINUOUS);

    UnsuitableRooms = m.addVar(0,sum(L[c] for c in Courses), 0,    GRB.CONTINUOUS);

    CurrAlone = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    CurrAssigned = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}

    SoftConstraint = m.addVar()
    Cost = m.addVar() #  Room Cost Weight in the Objective

    StudentLoad = {(cu, d): m.addVar() for cu in Curricula for d in range(len(D))}
    StudentLoadViol = {(cu, d): m.addVar() for cu in Curricula for d in range(len(D))}
    CurrDayInUse = {(cu, d): m.addVar() for cu in Curricula for d in range(len(D))}

    StudentMinMaxLoad = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)
    MinimumWorkingDaysBelow = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)
    VarCurrCompactViol = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)
    VarBadTimeSlotsViol = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)

    ConsPenalty = m.addVar(0, 0, 0, GRB.CONTINUOUS)

    Proximity = m.addVar(lb=0, ub=GRB.INFINITY, obj=0, vtype=GRB.CONTINUOUS, name="proximity");
    
    AllLecPlanned = {c : m.addConstr(quicksum(Y[c, p] for p in P if (c, p) in Y) == L[c]) for c in Courses}

    #Additional Room Cuts
    m.addConstr(len(P) * quicksum(R[s] for s in S) >= sum(L[c] for c in Courses))

    {(s) :  m.addConstr(quicksum(R[s] for s in SGT[s]) == RPlus[s]) for s in S}


    if (UseRoomHallConditions) :
        {(s) : m.addConstr(RPlus[s] * len(P) >= sum(L[c] for c in CGT[s])) for s in S}


    if (UseHallsConditions):
        {p : m.addConstr(quicksum(Y[c, p] for c in Courses if (c, p) in Y) <= quicksum(R[s] for s in S)) for p in P}

        {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y)<= RPlus[s]) for p in P for s in S}

    

    # Lecturers
    TeacherOneCourse = {(t, p): m.addConstr(quicksum(Y[c, p] for c in T[t] if (c, p) in Y) <= 1) for t in T for p in P}

    {(cu, p) : m.addConstr(quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)  <= 1) for cu in Curricula for p in P}


    {(c, d) : m.addConstr(quicksum(Y[c, p] for p in D[d] if (c, p) in Y) - DayInUse[c, d] >= 0) for c in Courses for d in range(len(D))}
    CalculationOfMNDViolation = {(c) : m.addConstr(quicksum(DayInUse[c, d] for d in range(len(D))) + DaysBelow[c] >= MND[c]) for c in Courses}

    {(cu, p) : m.addConstr(CurrAssigned[cu, p] ==  quicksum(Y[c, p] for c in Curricula[cu] if (c,p ) in Y)) for cu in Curricula for p in P}

    compactness = {(cu, D[d][pi]) : m.addConstr(CurrAssigned[cu, D[d][pi]] - 
            quicksum(CurrAssigned[cu, p] for p in TimeslotConseq(cu, pi, d)) <=  
            CurrAlone[cu, D[d][pi]]) for cu in Curricula for d in range(len(D)) for pi in range(len(D[d]))}
    

    {(cu, d, p) : m.addConstr(CurrDayInUse[cu, d] >= CurrAssigned[cu, p]) for cu in Curricula for d in range(len(D)) for p in D[d]}
    {(cu, d) : m.addConstr(quicksum(CurrAssigned[cu, p] for p in D[d]) == StudentLoad[cu, d]) for cu in Curricula for d in range(len(D))}
    # {(cu, d) : m.addConstr(StudentLoad[cu, d] + StudentLoadViol[cu, d] >= MinPPD * CurrDayInUse[cu, d]) for cu in Curricula for d in range(len(D))}
    # {(cu, d) : m.addConstr(StudentLoad[cu, d] - StudentLoadViol[cu, d] <= MaxPPD) for cu in Curricula for d in range(len(D))}

    {(p, c) : m.addConstr(Y[p, c] <= TimeslotUsed[p]) for p in P for c in Courses if (p, c) in Y}
       
    ConsecPeriods = {p: m.addConstr((TimeslotUsed[p] - TimeslotUsed[p - 1]) <= 0) for p in P if p != 1}

    # m.addConstr(StudentMinMaxLoad == quicksum(StudentLoadViol[x] for x in StudentLoadViol))
    m.addConstr(MinimumWorkingDaysBelow == quicksum(DaysBelow[c] for c in Courses))
    m.addConstr(VarCurrCompactViol == quicksum(CurrAlone[x] for x in CurrAlone))

    # m.setObjective(quicksum(s * R[s] for s in S))
    # 5 MWD Violation
    # 2 Curriculum Compact 

    m.setObjective(5 * MinimumWorkingDaysBelow + 2 * VarCurrCompactViol)
    m.optimize()

    # print([(c, CurrAlone[c].x) for c in CurrAlone if CurrAlone[c].x > 0.9])
    
    # print(MinimumWorkingDaysBelow.x)
    # print(VarCurrCompactViol.x)
    # for c in Courses:
    #     for d in range(len(D)):
    #         for p in range(len(D[d])):
    #             if ((c,D[d][p]) in Y and Y[c, D[d][p]].x > 0.9):
    #                 print(c, d, p)
    # print(len(Y)) # Right
    # print(len(R)) # Right
    # print(len(RPlus)) # Right
    # print(len(DayInUse)) # Right

   # print(len(TimeslotUsed)) Right
    # print(len(DaysBelow))

    # print(len(CurrAlone)) # Should be 57
    # print(len(CurrAssigned)) # Should be 57 

    # print(len(StudentLoad))
    # print(len(StudentLoadViol))
    # for s in S:
    #     print(s, R[s].x)
    # print(m.objVal, end="")

    
def CCTModelForRoomPlanning(output=True) :
    UseHallsConditions = True
    UseRoomHallConditions = True


    step = 25
    maxS = max([math.ceil(DEM[c] / step) * step for c in DEM])
    S = [i for i in range(step, maxS + 1, step)]


    SGT = {s : [ss for ss in S if ss > s] for s in S}
    CGT = {s : [c for c in Courses if DEM[c] > s] for s in S}
    RGT = {s : [r for r in Rooms if CAP[r] > s] for s in S}


    m = Model()
    if (not output):
        m.setParam('OutputFlag', 0)

    Y = {(c, p): m.addVar(vtype = GRB.BINARY) for c in Courses for p in P if p not in PeriodConstraints[c]}
    R = {(s) : m.addVar(vtype=GRB.INTEGER) for s in S}
    RPlus = {(s) : m.addVar(vtype=GRB.INTEGER) for s in S}
    
    DayInUse = {(c, d) : m.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    DaysBelow = {(c) : m.addVar() for c in Courses}

    TimeslotUsed = {p : m.addVar(vtype = GRB.BINARY) for p in P}

    RoomStability = m.addVar(lb=0,ub=len(Courses)*(len(Rooms) - 1), obj=0, vtype=GRB.CONTINUOUS);

    UnsuitableRooms = m.addVar(0,sum(L[c] for c in Courses), 0,    GRB.CONTINUOUS);

    CurrAlone = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    CurrAssigned = {(cu, p) : m.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}

    SoftConstraint = m.addVar()
    Cost = m.addVar() #  Room Cost Weight in the Objective

    StudentLoad = {(cu, d): m.addVar() for cu in Curricula for d in range(len(D))}
    StudentLoadViol = {(cu, d): m.addVar() for cu in Curricula for d in range(len(D))}
    CurrDayInUse = {(cu, d): m.addVar() for cu in Curricula for d in range(len(D))}

    StudentMinMaxLoad = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)
    MinimumWorkingDaysBelow = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)
    VarCurrCompactViol = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)
    VarBadTimeSlotsViol = m.addVar(lb=0, ub=GRB.INFINITY, obj=0)

    ConsPenalty = m.addVar(0, 0, 0, GRB.CONTINUOUS)

    Proximity = m.addVar(lb=0, ub=GRB.INFINITY, obj=0, vtype=GRB.CONTINUOUS, name="proximity");
    
    AllLecPlanned = {c : m.addConstr(quicksum(Y[c, p] for p in P if (c, p) in Y) == L[c]) for c in Courses}

    #Additional Room Cuts
    m.addConstr(len(P) * quicksum(R[s] for s in S) >= sum(L[c] for c in Courses))

    {(s) :  m.addConstr(quicksum(R[s] for s in SGT[s]) == RPlus[s]) for s in S}


    if (UseRoomHallConditions) :
        {(s) : m.addConstr(RPlus[s] * len(P) >= sum(L[c] for c in CGT[s])) for s in S}


    if (UseHallsConditions):
        {p : m.addConstr(quicksum(Y[c, p] for c in Courses if (c, p) in Y) <= quicksum(R[s] for s in S)) for p in P}

        {(s, p) : m.addConstr(quicksum(Y[c, p] for c in CGT[s] if (c, p) in Y)<= RPlus[s]) for p in P for s in S}

    

    # Lecturers
    TeacherOneCourse = {(t, p): m.addConstr(quicksum(Y[c, p] for c in T[t] if (c, p) in Y) <= 1) for t in T for p in P}

    {(cu, p) : m.addConstr(quicksum(Y[c, p] for c in Curricula[cu] if (c, p) in Y)  <= 1) for cu in Curricula for p in P}


    {(c, d) : m.addConstr(quicksum(Y[c, p] for p in D[d] if (c, p) in Y) - DayInUse[c, d] >= 0) for c in Courses for d in range(len(D))}
        
        
    CalculationOfMNDViolation = {(c) : m.addConstr(quicksum(DayInUse[c, d] for d in range(len(D))) + DaysBelow[c] >= MND[c]) for c in Courses}

    # TODO Fix this constraint up
    # {(cu, p) : m.addConstr(CurrAssigned[cu, p] == )}
    # g1 = {(cu,p):  m.addConstr(-CurrAssigned[cu, p-1] + CurrAssigned[cu,p] - CurrAssigned[cu, p+1] - CurrAlone[cu, p] <= 0) for (cu, p) in CurrAssigned if (cu, p-1) in CurrAssigned and (cu, p+1) in CurrAssigned}
    # g2 = {(cu,p):  m.addConstr(CurrAssigned[cu,p] - CurrAssigned[cu, p+1] - CurrAssigned[cu, p] <= 0) for (cu, p) in CurrAssigned if (cu, p-1)  not in CurrAssigned}
    # g3 = {(cu,p):  m.addConstr(-CurrAssigned[cu, p-1] + CurrAssigned[cu,p]  - CurrAssigned[cu, p] <= 0) for (cu, p) in CurrAssigned if (cu, p+1) not in CurrAssigned}


    {(cu, d, p) : m.addConstr(CurrDayInUse[cu, d] >= CurrAssigned[cu, p]) for cu in Curricula for d in range(len(D)) for p in D[d]}

    # {(cu, d) : m.addConstr(quicksum(CurrAssigned[cu][p] for p in D[d]) == StudentLoad[cu, d]) for cu in Curricula for d in range(len(D)) }


    {(p, c) : m.addConstr(Y[p, c] <= TimeslotUsed[p]) for p in P for c in Courses if (p, c) in Y}
       
    # ConsecPeriods = {p: m.addConstr((TimeslotUsed[p] - TimeslotUsed[p - 1]) <= 0) for p in P if p != 1}

    m.addConstr(StudentMinMaxLoad == quicksum(StudentLoadViol[x] for x in StudentLoadViol))
    
    m.addConstr(MinimumWorkingDaysBelow == quicksum(DaysBelow[c] for c in Courses))

    m.addConstr(VarCurrCompactViol == quicksum(CurrAlone[x] for x in CurrAlone))

    m.setObjective(quicksum(s * R[s] for s in S))
    m.optimize()

if __name__ == "__main__":
    TeachingPeriodsvsQuality()