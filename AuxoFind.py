import cobra
from cobra import Reaction, Metabolite, Model, Gene
from collections import defaultdict
from warnings import warn
import math

def get_CEGs(cobra_model, growth_medium = []):
    """Identifies the set Conditionally Essential Genes (CEGs) required for growth in user-defined media by comparing gene essentiality in rich media. 
    
    Parameters
    ----------
    cobra_model : class:`~cobra.core.Model.Model` object
        The model to evaluate.
    growth_medium: list, optional
        List of exchange reactions modeling the growth medium. The default is the set of exchange reactions which are set to allow influx of extracellular components in the input model.


    Returns
    -------
    list
        The list of conditionally essential genes.
        
    References
    ----------
    .. [1] 
"""       
    check_BOF(cobra_model)
    
    growth_medium = [x.id for x in cobra_model.reactions if x.id.startswith('EX') if x.lower_bound < 0] if growth_medium == [] else growth_medium
        
    with cobra_model:
        for rx in cobra_model.reactions:
            if rx.id.startswith('EX_'):
                rx.lower_bound = 0
            if rx.id in growth_medium:
                rx.lower_bound = -10
    
        check_solution(cobra_model)
        
    
    MEGs = cobra.flux_analysis.find_essential_genes(cobra_model)
    with cobra_model:
        for rx in cobra_model.reactions:
            if rx.id.startswith('EX_'):
                rx.lower_bound = -10
        EGs = cobra.flux_analysis.find_essential_genes(cobra_model)
    CEGs = {x.id for x in MEGs - EGs}
    
    
    return CEGs

def check_BOF(cobra_model):
    """ Checks that the model has a single objective function."""
    BOF = [x.id for x in cobra_model.reactions if x.objective_coefficient == 1]
    if len(BOF) == 0:
        raise ValueError('The baseline model does not have any set objective function.')
        return
    if len(BOF) > 1:
        warn("The model has more than one objective function.", UserWarning)
    return

def check_solution(cobra_model):
    """ Checks that the constraint-based optimization has a non-negative non-zero non-null solution """
    sol = cobra_model.slim_optimize()
    if sol < 0.001 or math.isnan(sol):
        raise ValueError('The baseline model cannot simulate growth under the provided growth medium conditions.')
        return  
    
    
def get_auxotrophies(cobra_model, missing_genes, lower_bound = 0, growth_medium = [], biomass = '', force_to_zero = [], method = 'optimal_solution'):
    """Identifies set of exchange(s) to open in order to achieve non-negative optimal yield for single gene knockout mutant. 
    
    Parameters
    ----------
    cobra_model : class:`~cobra.core.Model.Model` object
        The model to evaluate.
    lower_bound : float, optional
        Minimum yield to be set as a lower bound the biomass objective function (defaults to half of the yield achieved by the baseline model).
    components_in_medium : list, optional
        List of exchange reactions modeling the baseline growth medium (without supplementation) (defaults to the lis of open exchanges).
    force_to_zero : list, optional
        List of exchange reactions to exclude from the solution space (defaults to empty list).
    biomass_reaction: str, optional
        The reaction ID for the biomass objective function (defaults to the first reaction with objective coefficient of 1).
    method: str, optional
        User can opt for single-search tree method to expose all possible epsilon-optimal alternative solutions using Gurobi’s PoolSearchMode (defaults to 'optimal_solution').
        
    Returns
    -------
    dict
        A dictionary mapping each CEG to a nutrient requirement.
        
    References
    ----------
    """  
    
    check_solution(cobra_model)
    check_BOF(cobra_model)
    
    lower_bound = cobra_model.slim_optimize()*0.5 if lower_bound == 0 else lower_bound  
    growth_medium = [x.id for x in cobra_model.reactions if x.id.startswith('EX') if x.lower_bound < 0] if growth_medium == [] else growth_medium
    biomass = [rx for rx in cobra_model.reactions if rx.objective_coefficient == 1][0].id if biomass == '' else biomass
    cobra_model.solver = 'gurobi'
    cobra_model.solver.configuration.tolerances.feasibility = 1e-9
    
    supplementations = {}

    for gene in missing_genes:
        with cobra_model:
            for rx in cobra.manipulation.find_gene_knockout_reactions(cobra_model, [gene]):
                rx.upper_bound = 0
                rx.lower_bound = 0
            sol = cobra_model.slim_optimize()
            if sol > 0.001:
                if not math.isnan(sol):
                    warn('The gene (%s) is not conditionally essential. Its knock-out has no effect on growth.'%gene)
                    continue
            
            supplementations[gene] = find_supplementations(cobra_model, lower_bound = lower_bound, growth_medium = growth_medium, biomass = biomass, force_to_zero = force_to_zero, method= method)

            
    return supplementations

def find_supplementations(cobra_model, lower_bound, growth_medium, biomass, force_to_zero, method):    
    
    """Finds required nutrient supplementations to support flux through the objective function once the model parameters have been set."""
    
    if method == 'single_search':
        cobra_model.solver.problem.setParam('PoolSearchMode', 2)
        cobra_model.solver.problem.setParam('PoolSolutions', 50)        
        
    indicators = []
    constraints = []
    for rx in cobra_model.reactions:
        if rx.id.startswith('EX_') and rx.id not in growth_medium+force_to_zero:
            rx.lower_bound = -10
            rx.upper_bound = 1000 
            
            indicator = cobra_model.problem.Variable('%s_iii'%rx.id , type = 'binary')
            indicators.append(indicator)

            new_cstr2 = cobra_model.problem.Constraint( -rx.flux_expression + rx.lower_bound*indicator ,ub = 0)
            constraints.append(new_cstr2)
            cobra_model.add_cons_vars([indicator, new_cstr2])
            
        elif rx.id in force_to_zero:
            rx.lower_bound = 0

    cobra_model.reactions.get_by_id(biomass).lower_bound = lower_bound

    cobra_model.objective = cobra_model.problem.Objective(-sum(indicators))
    solution = cobra_model.optimize(objective_sense="max")
    
    if method == 'single_search':
        d = {' AND '.join(x[1:]).replace('_iii',''):x[0] for x in get_active_decision_vars(cobra_model).values()}
        if '' in d.keys():
            del d['']
        indicator_results = [sols.split(' AND ') for sols, s in d.items() if s == max(d.values())]
    else:
        indicator_results = {ind.name[:-4]: ind.primal for ind in indicators if ind.primal != 0.0}.keys()
        
    cobra_model.remove_cons_vars(constraints+indicators)
    return indicator_results

def get_active_decision_vars(milp_problem):
    
    """Based on MILP solution find decision variables active"""
    out = defaultdict(list)

    num_sols = milp_problem.solver.problem.SolCount
    print('Number of solutions:', num_sols)
    for i in range(0, num_sols):
        milp_problem.solver.problem.setParam('solutionNumber', i)
        # Append the objective value
        out[i].append(milp_problem.solver.problem.PoolObjVal)
        for var in milp_problem.solver.problem.getVars():
            if var.Xn == 1 and var.VType == 'B':
                r_id = var.VarName
                out[i].append(r_id)

    return out
