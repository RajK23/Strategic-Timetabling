import xml.etree.ElementTree as ET
import os
from gurobipy import * 
import operator
from collections import defaultdict

script_dir = os.path.dirname(__file__)
relpath = 'CCT/ITC2007/comp06.xml'
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
deltaMode = True

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
            # if not deltaMode:
            CAP[id] = int(attribs['size'])
            S.add(int(attribs['size']))
            # else:
                # CAP[id] = math.ceil(int(attribs['size']) / 25) * 25
                # S.add(CAP[id])

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
SGT = {s : [ss for ss in S if ss >= s] for s in S}
CGT = {s : [c for c in Courses if DEM[c] >= s] for s in S}
RGT = {s : [r for r in Rooms if CAP[r] >= s] for s in S}

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

def TimeslotConstraint(model, X):
    for course in PeriodConstraints:
        for day, period in PeriodConstraints[course]:
            con = day*Days + period + 1
            model.addConstr(X[course, con] == 0) 

def RoomPlanningProblem() :
    # need to implement the linear version of this algorithm maybe in the future
    RPP  = Model()
    R = {s : RPP.addVar(vtype= GRB.INTEGER) for s in S}
    RPP.setObjective(quicksum(s * R[s] for s in S))
    SetCovering = {s :  RPP.addConstr(len(P) * quicksum(R[ss]  for ss in SGT[s]) >= quicksum(L[c] for c in CGT[s])) for s in S}
    
    A ={s : RPP.addConstr(R[s] == sum(1 for r in Rooms if CAP[r] == s)) for s in S}

    RPP.optimize();

    for s in S:
        print(s, R[s].x)

#Can be solved if the Dem is pre sorted
def TeachingPeriodsProblem() :
    # Teaching Periods Problem
    TPP = Model()
    TP = {p : TPP.addVar(vtype = GRB.BINARY) for p in P}

    TPP.setObjective(quicksum(TP[p] for p in P), GRB.MINIMIZE)
    # TPP.addConstr(len(Rooms) * quicksum(TP[p] for p in TP) >= sum(L[c] for c in Courses))
    #SetCovering = {s: TPP.addConstr(len(RGT[s])*quicksum(TP[p] for p in P) >= quicksum(L[c] for c in CGT[s])) for s in S}
    
    # TPP.addConstr(len(Rooms) * quicksum(TP[p] for p in TP) >= sum(L[c] for c in Courses))
    for s in S:
        if s <= 30:
            TPP.addConstr((len(RGT[s]) - 1 )*quicksum(TP[p] for p in P) >= quicksum(L[c] for c in CGT[s]))        
        else:    
            TPP.addConstr((len(RGT[s])) * len(P)  >= sum(L[c] for c in Courses if DEM[c] >= s))

    # quicksum(TP[p] for p in P)
    ConsecPeriods = {p: TPP.addConstr((TP[p+1] - TP[p]) <= 0) for p in P if p+1 in P}
    TPP.optimize();

    for s in S:
        pass
        #print(s, (len(RGT[s])) ,(len(RGT[s])) * sum(TP[p].x for p in P) ,'>=', sum(L[c] for c in Courses if DEM[c] >= s))

    for p in P:
        if TP[p].x > 0.9:
            print("Allocated: " +  str(p))

def RoomsPlanningvsQuality():
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


    c4 = {(s, p): RPQ.addConstr(quicksum(X[c, p] for c in CGT[s])  <= quicksum(R[ss] for ss in SGT[s])) for s in S for p in P}
    #4e
    AllLecPlanned = {c : RPQ.addConstr(quicksum(X[c, p] for p in P) == L[c]) for c in Courses}
    #Must set z to one to be able to plannign a lecture in this class
    d = {(c, d) : RPQ.addConstr(quicksum(X[c, p] for p in D[d]) - Z[c, d] >= 0) for c in Courses for d in range(len(D))}
    e = {(c) : RPQ.addConstr(quicksum(Z[c, d] for d in range(len(D))) + W[c] >= MND[c]) for c in Courses}
    f = {(cu, p) : RPQ.addConstr(quicksum(X[c, p] for c in Curricula[cu]) - Q[cu, p]  == 0) for cu in Curricula for p in P}
    
    #Compactness
    g1 = {(cu,p):  RPQ.addConstr(-Q[cu, p-1] + Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1) in Q and (cu, p+1) in Q}
    g2 = {(cu,p):  RPQ.addConstr(Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1)  not in Q}
    g3= {(cu,p):  RPQ.addConstr(-Q[cu, p-1] + Q[cu,p]  - V[cu, p] <= 0) for (cu, p) in Q if (cu, p+1) not in Q}
    #Teacher can
    h = {(t, p): RPQ.addConstr(quicksum(X[c, p] for c in T[t]) <= 1) for t in T for p in P}

    #4f
    SetCovering = {s :  RPQ.addConstr(len(P) * quicksum(R[ss] * ss for ss in SGT[s]) >= quicksum(L[c] for c in CGT[s])) for s in S}
    A = {s : RPQ.addConstr(R[s] == sum(1 for r in Rooms if CAP[r] == s)) for s in S}

    BoundsSolver(fSeats, fQual, RPQ, "Room Planning", "Quality")
    # BiObjectiveSolver( fQual, fSeats, RPQ, 1)
    # RPQ.setObjective(fQual);
    # RPQ.optimize()

def RoomsPlanningvsTeachingPeriods():
    RPTP = Model ()
    fTime = quicksum(TP[p] for p in P)
    fSeats = quicksum(s * R[s] for s in S)
    X = {(c, p) : RPTP.addVar(vtype = GRB.BINARY) for c in Courses for p in P}
    Z = {(c, d) : RPTP.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    W = {(c) : RPTP.addVar(vtype = GRB.INTEGER) for c in Courses}
    Q = {(cu, p) : RPTP.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    V = {(cu, p) : RPTP.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    R = {s : RPTP.addVar(vtype= GRB.INTEGER) for s in S}
    TP = {p : RPTP.addVar(vtype = GRB.BINARY) for p in P}


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

    BiObjectiveSolver(fTime, fSeats, m, 1);

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
    # TPQ.setParam('OutputFlag', 0)
    # Look at adding in class constraints.
    # BiObjectiveSolver(fTime, fQual, TPQ, 1)
    BoundsSolver(fTime, fQual, TPQ, "TimePeriods", "Quality")

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

if __name__ == "__main__":
    RoomPlanningProblem()
