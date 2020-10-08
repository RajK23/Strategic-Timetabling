

from gurobipy import *


m = Model() 

x1 = m.addVar() 
x2 = m.addVar()


u1 = x1
u2 = 3*x1 + 4*x2

m.addConstr(x1 <= 20)
m.addConstr(x2 <= 40)
m.addConstr(5 * x1 + 4* x2 <= 200)


m.setObjective(u2, GRB.MAXIMIZE)

epsilon = 8
delta = 1

m.addConstr(u1 == 9)

m.optimize()

print(x1.x)
print(x2.x)
print('--------')
print(u1.x)
print(u2.getValue())

# def BiObjectiveSolver(fx, fy, m, delta) :
#     m.setParam('OutputFlag', 0)
#     m.setObjective(fy, GRB.MAXIMIZE);
#     m.optimize();
#     minY = m.objVal

#     TEMP = m.addConstr(fy == minY)    
#     m.setObjective(fx, GRB.MAXIMIZE)
#     m.optimize()

#     #max x ← Minimize (f_x)
#     maxX = m.objVal;
#     m.remove(TEMP)

#     # Minimize f(x) again without the constraint
#     m.optimize()
#     epsilon = m.objVal
#     epsilonConstraint = m.addConstr(fx >= epsilon)
#     m.setObjective(fy, GRB.MAXIMIZE)


#     print("MIN Y", minY)
#     print("MAX X", minY)
#     print("EPSILON" , epsilon)
#     while True:
#         m.optimize()

#         #print('fyHat=', m.objVal)
#         fyHat = m.objVal
#         epsilon = epsilon + delta
#         #print('fyHat= ', x1.x)
#         print('x1=', x1.x)
#         print('x2=', x2.x)
#         epsilonConstraint.RHS = epsilon + delta
        
#         if (epsilon > maxX or fyHat <= minY):
#             print("BREAKING")
#             print(epsilon , "   ", maxX)
#             print(fyHat , "   ", minY)
#             break;



# BiObjectiveSolver(u1, u2, m, 0.01)