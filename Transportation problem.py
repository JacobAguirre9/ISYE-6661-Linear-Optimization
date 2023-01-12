# This example is taken from page 12 of the lecture notes found in this repository

import gurobipy as gp


def transportation_problem(I, J, s, d, c):
    """
    Solve the transportation problem using Gurobi.

    Parameters:
    - I: number of warehouses
    - J: number of customers
    - s: list of supply for each warehouse
    - d: list of demand for each customer
    - c: 2D list of costs for each warehouse-customer pair

    Returns:
    - A tuple containing the optimal objective value and a 2D list of flows on each warehouse-customer pair
    """
    model = gp.Model()
    # Create variables for the flow on each warehouse-customer pair
    x = model.addVars(I, J, lb=0, vtype=gp.GRB.CONTINUOUS, name="x")

    # Set the objective function
    model.setObjective(gp.quicksum(c[i][j] * x[i, j] for i in range(I) for j in range(J)), gp.GRB.MINIMIZE)

    # Add supply constraints
    for i in range(I):
        model.addConstr(gp.quicksum(x[i, j] for j in range(J)) <= s[i])

    # Add demand constraints
    for j in range(J):
        model.addConstr(gp.quicksum(x[i, j] for i in range(I)) >= d[j])

    # Add non-negativity constraint
    for i in range(I):
        for j in range(J):
            model.addConstr(x[i, j] >= 0)
    
    # Optimize the model
    model.optimize()

    # Extract the flows on each warehouse-customer pair
    flows = [[x[i, j].x for j in range(J)] for i in range(I)]

    return model.objVal, flows


  
  """
  
  Tractable example 
  
  # Define the number of warehouses and customers
I = 2
J = 2

# Define the supplies and demands
s = [100, 250]
d = [200, 50]

# Define the costs
c = [[4, 8], [6, 7]]

# Solve the transportation problem
obj_val, flows = transportation_problem(I, J, s, d, c)

# Print the results
print("Optimal objective value:", obj_val)
print("Flows on each warehouse-customer pair:", flows)

  
  """
