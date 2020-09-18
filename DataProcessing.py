import xml.etree.ElementTree as ET
import os
from gurobipy import * 
import operator
from collections import defaultdict

script_dir = os.path.dirname(__file__)
relpath = 'CCT/ITC2007/comp01.xml'
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


P = range(1, Days * PPD + 1)
D = [[(j+1) + i*(Days) for i in range(PPD)] for j in range(Days)]
SGT = {s : [ss for ss in S if ss >= s] for s in S}
CGT = {s : [c for c in Courses if DEM[c] >= s] for s in S}
RGT = {s : [c for c in Rooms if CAP[c] >= s] for s in S}

def QualityProblem() :
    Quality = Model()
    X = {(c, p) : Quality.addVar(vtype = GRB.BINARY) for c in Courses for p in P}
    #print(X)
    #Minimum working days
    Z = {(c, d) : Quality.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    W = {(c) : Quality.addVar(vtype = GRB.INTEGER) for c in Courses}
    Q = {(cu, p) : Quality.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    V = {(cu, p) : Quality.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}

    Quality.setObjective(quicksum(5 * W[c] for c in Courses) + quicksum(2 * V[cu,p] for cu in Curricula for p in P));
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
    g3= {(cu,p):  Quality.addConstr(-Q[cu, p-1] + Q[cu,p]  - V[cu, p] <= 0) for (cu, p) in Q if (cu, p+1) not in Q}
    #Teacher can
    h = {(t, p): Quality.addConstr(quicksum(X[c, p] for c in T[t]) <= 1) for t in T for p in P}

    #print(len(f));
    Quality.optimize();

    # for d in D:
    #     for p in d:
    #         print("Timeslot: ", str(p));
    #         for c in Courses:
    #             if (X[c, p].x > 0.9):
    #                 print("Course", c, "is planned in Timeslot" , p)

    #
    # for a in X:
    #     if (X[a].x > 0.9):
    #         print(a, X[a].x)
    #
    # print("----------")
    #for c in Courses:
        #if (W[c].x > 0.9):
        #print(c, W[c].x)
    # for c, d in Z:
    #     #if (Z[c, d].x > 0.9):
    #     print(sum(Z[c, d].x for d in range(len(D))), MND[c])

    # for a in Q:
    #     if (Q[a].x > 0.9):
    #         print(a, Q[a].x)

    # for a in V:
    #     if (V[a].x > 0.9):
    #         print(V[a].x)
    
def RoomPlanningProblem() :
    # need to implement the linear version of this algorithm maybe in the future
    RPP  = Model()
    R = {s : RPP.addVar(vtype= GRB.INTEGER) for s in S}
    RPP.setObjective(quicksum(s * R[s] for s in S))
    SetCovering = {s :  RPP.addConstr(len(P) * quicksum(R[ss] * ss for ss in SGT[s]) >= quicksum(L[c] for c in CGT[s])) for s in S}
    
    RPP.optimize();

    # for s in S:
    #     #if R[s].x != 0:
    #     print("Allocated  " +  str(int(R[s].x)) + " of size " + str(int(s)))

#Can be solved if the Dem is pre sorted
def TeachingPeriodsProblem() :
    # Teaching Periods Problem
    TPP = Model()
    TP = {p : TPP.addVar(vtype = GRB.BINARY) for p in P}
    TPP.setObjective(quicksum(TP[p] for p in P), GRB.MINIMIZE)
    SetCovering = {s: TPP.addConstr(len(RGT[s])*quicksum(TP[p] for p in P) >= quicksum(L[c] for c in CGT[s])) for s in S}
    ConsecPeriods = {p: TPP.addConstr((TP[p] - TP[p-1]) <= 0) for p in P if p != 1}
    return TPP;
    # for p in P:
    #     if T[p].x > 0.9:
    #         print("Allocated: " +  str(p))

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

    RPQ.setObjective(fQual);
    RPQ.optimize()

def TeachingPeriodsvsQuality():
    TPQ = Model()
    TP = {p : TPQ.addVar(vtype = GRB.BINARY) for p in P}
    X = {(c, p) : TPQ.addVar(vtype = GRB.BINARY) for c in Courses for p in P}
    #print(X)
    #Minimum working days
    Z = {(c, d) : TPQ.addVar(vtype = GRB.BINARY) for c in Courses for d in range(len(D))}
    W = {(c) : TPQ.addVar(vtype = GRB.INTEGER) for c in Courses}
    Q = {(cu, p) : TPQ.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}
    V = {(cu, p) : TPQ.addVar(vtype = GRB.BINARY) for cu in Curricula for p in P}


    fTime = quicksum(TP[p] for p in P)
    fQual = quicksum(5 * W[c] for c in Courses) + quicksum(2 * V[cu,p] for cu in Curricula for p in P)

    #5c
    c5 = {(c, p): TPQ.addConstr(X[c, p] - TP[p] <= 0) for c in Courses for p in P}
    
    RoomsAccomodate = {(s, p) : TPQ.addConstr(quicksum(X[c, p] for c in CGT[s]) <= len(RGT[s])) for s in S for p in P}
    AllLecPlanned = {c : TPQ.addConstr(quicksum(X[c, p] for p in P) == L[c]) for c in Courses}
    #Must set z to one to be able to plannign a lecture in this class
    d = {(c, d) : TPQ.addConstr(quicksum(X[c, p] for p in D[d]) - Z[c, d] >= 0) for c in Courses for d in range(len(D))}
    e = {(c) : TPQ.addConstr(quicksum(Z[c, d] for d in range(len(D))) + W[c] >= MND[c]) for c in Courses}
    f = {(cu, p) : TPQ.addConstr(quicksum(X[c, p] for c in Curricula[cu]) - Q[cu, p]  == 0) for cu in Curricula for p in P}
    
    #Compactness
    g1 = {(cu,p):  TPQ.addConstr(-Q[cu, p-1] + Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1) in Q and (cu, p+1) in Q}
    g2 = {(cu,p):  TPQ.addConstr(Q[cu,p] - Q[cu, p+1] - V[cu, p] <= 0) for (cu, p) in Q if (cu, p-1)  not in Q}
    g3= {(cu,p):  TPQ.addConstr(-Q[cu, p-1] + Q[cu,p]  - V[cu, p] <= 0) for (cu, p) in Q if (cu, p+1) not in Q}
    #Teacher can
    h = {(t, p): TPQ.addConstr(quicksum(X[c, p] for c in T[t]) <= 1) for t in T for p in P}

    SetCovering = {s: TPQ.addConstr(len(RGT[s])*quicksum(TP[p] for p in P) >= quicksum(L[c] for c in CGT[s])) for s in S}
    ConsecPeriods = {p: TPQ.addConstr((TP[p] - TP[p-1]) <= 0) for p in P if p != 1}
    
    TPQ.setObjective(fQual);

    TPQ.optimize();

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
    g3= {(cu,p):  RPTP.addConstr(-Q[cu, p-1] + Q[cu,p]  - V[cu, p] <= 0) for (cu, p) in Q if (cu, p+1) not in Q}
    #Teacher can
    h = {(t, p): RPTP.addConstr(quicksum(X[c, p] for c in T[t]) <= 1) for t in T for p in P}

    SetCoveringR = {s :  RPTP.addConstr(len(P) * quicksum(R[ss] * ss for ss in SGT[s]) >= quicksum(L[c] for c in CGT[s])) for s in S}

    SetCoveringTP = {s: RPTP.addConstr(len(RGT[s])*quicksum(TP[p] for p in P) >= quicksum(L[c] for c in CGT[s])) for s in S}
    ConsecPeriods = {p: RPTP.addConstr((TP[p] - TP[p-1]) <= 0) for p in P if p != 1}

    c4 = {(s, p): RPTP.addConstr(quicksum(X[c, p] for c in CGT[s])  <= quicksum(R[ss] for ss in SGT[s])) for s in S for p in P}
    c5 = {(c, p): RPTP.addConstr(X[c, p] - TP[p] <= 0) for c in Courses for p in P}

    RPTP.setObjective(fTime)

    RPTP.optimize();

#TODO : Final Multi-Objective Stuff
# Room planning vs Teaching Periods
# def RP_VS_TP():
#     RPTP = Model()
#     T = {p : RPTP.addVar(vtype = GRB.BINARY) for p in P}
#     R = {s : RPTP.addVar(vtype= GRB.INTEGER) for s in S}

#     # add in constraints
    

#     #Two objectives
#     fx = quicksum(T[p] for p in P)
#     fy = quicksum(s * R[s] for s in S)
    
    
#     RPTP.setObjective(fy);
#     RPTP.optimize() 
#     minY = RPTP.objVal

#     tbr = RPTP.addConstr(fy == minY) # tbr


#     RPTP.setObjective(fx)
#     RPTP.optimize()
#     maxX = RPTP.objVal;
#     RPTP.remove(tbr)
#     RPTP.optimize()
#     epsilon = RPTP.objVal
#     epsilonConstraint = RPTP.addConstr(fx <= epsilon)


#     # Do while loop 

#     while True:
#         RPTP.optimize()

#         fyHat = RPTP.objVal
#         epsilon = epsilon + delta
#         epsilonConstraint.RHS = epsilon + delta
        

#         if (epsilon > maxX or fyHat <= minY):
#             break;



    





if __name__ == "__main__":
    RoomsPlanningvsQuality()