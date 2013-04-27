#!/usr/bin/env python
from KStandard import attributesFromDict
from numpy import *
from pulp import *


animals=["Cow","Pig","Sheep","Goose","Chicken"]
Money=100
TotalAnimals=100
costs={"Cow":10,
       "Pig":0.5,
       "Sheep":3,
       "Goose":0.33,
       "Chicken":0.69}

# Create the 'prob' variable to contain the problem data
prob = LpProblem("The Animal Purchasing Problem", LpMaximize)

animal_vars = LpVariable.dicts("Animals",animals,1,100,LpInteger)

# The objective function is added to 'prob' first
prob += lpSum([costs[i]*animal_vars[i] for i in animals]), "Total Cost of Animals"
# Then the constraints.
#for i in animals:
#    prob += lpSum([animal_vars[i]]) >= 1, "At least one of each animal"
prob += lpSum([animal_vars[i] for i in animals]) == TotalAnimals, "Total number of animals"
prob += lpSum([costs[i]*animal_vars[i] for i in animals]) <= Money, "Total Cost of Animals"

# The problem data is written to an .lp file
prob.writeLP("AnimalProblem.lp")

# The problem is solved using PuLP's choice of Solver
prob.solve()

# The status of the solution is printed to the screen
print "Status:", LpStatus[prob.status]

print "The animal costs are:"
for k,v in costs.iteritems():
    print "\t",k,":",v
print ""
print "To maximize the cost while buying 100 animals, the animals we should buy are:"
# Each of the variables is printed with it's resolved optimum value
for v in prob.variables():
    print "\t",v.varValue," ",v.name

# The optimised objective function value is printed to the screen    
print "\nTotal Cost is $", value(prob.objective)


