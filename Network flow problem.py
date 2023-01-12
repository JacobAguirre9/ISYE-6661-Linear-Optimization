import gurobipy as gp

def network_flow(nodes, arcs, capacities, source, sink, flows=None):
    """
    Compute the maximum flow in a network using Gurobi.

    Parameters:
    - nodes: list of nodes in the network
    - arcs: list of tuples representing the arcs in the network
    - capacities: dictionary with arc-capacity pairs
    - source: the source node
    - sink: the sink node
    - flows: dictionary with arc-flow pairs

    Returns:
    - A tuple containing the optimal objective value and a dictionary of flows on each arc
    """
    model = gp.Model()
    # Create variables for the flow on each arc
    if flows:
        flow = model.addVars(arcs, lb=0, ub=capacities, obj=flows, name="flow", vtype=gp.GRB.CONTINUOUS)
    else:
        flow = model.addVars(arcs, lb=0, ub=capacities, name="flow", vtype=gp.GRB.CONTINUOUS)
    # Add the flow conservation constraints
    for n in nodes:
        incoming = [(n, v) for (u, v) in arcs if v == n]
        outgoing = [(u, n) for (u, v) in arcs if u == n]
        if incoming and outgoing:
            model.addConstr(gp.quicksum(flow[u, n] for (u, n) in incoming) - 
                            gp.quicksum(flow[n, v] for (n, v) in outgoing) == 0)
        elif incoming:
            model.addConstr(gp.quicksum(flow[u, n] for (u, n) in incoming) == 0)
        elif outgoing:
            model.addConstr(gp.quicksum(flow[n, v] for (n, v) in outgoing) == 0)
    # Add the source and sink constraints
    model.addConstr(flow.sum(source,'*') == flow.sum('*',sink))
    # Optimize the model
    model.optimize()

    # Extract the flows on each arc
    arc_flows = {(i, j): flow[i, j].x for i, j in arcs}

    return model.objVal, arc_flows
  
  # You can then pass the nodes, arcs, capacities, source and sink of a network as well as any initial flows for the arcs to this 
  # function to get the solution for the network flow problem.

  # Tractable example below

  """
  # Define the nodes, arcs, capacities, source, and sink of a simple network
  nodes = [1, 2, 3, 4]  
  arcs = [(1, 2), (1, 3), (2, 3), (3, 4)]
  capacities = {(1, 2): 2, (1, 3): 4, (2, 3): 1, (3, 4): 2}
  source = 1
  sink = 4

  # Compute the maximum flow
  obj_val, flows = network_flow(nodes, arcs, capacities, source, sink)

  # Print the results
  print("Optimal objective value:", obj_val)
  print("Flows on each arc:", flows)

  """
