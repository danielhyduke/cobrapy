from __future__ import with_statement
#cobra.flux_analysis.parsimonious.py
#Implement parsimonious flux balance analysis as described in {Lewis et al., 2010, Mol Syst Biol, 6, 390}
from ..manipulation.modify import convert_to_irreversible
from cobra import Metabolite, Reaction
def flux_balance_analysis(model, minimize_spontaneous_reactions=False, linear_objective=None):
    """Implements parsimonious flux balance analysis as described in {Lewis et al., 2010, Mol Syst Biol, 6, 390}

    Q: Would this be better as an option to Model.optimize() so we could allow for more control of what
    type of initial optimization is performed?

    model: :class:`~cobra.cobra.Model`

    minimize_spontaneous_reactions: Boolean.  If True then minimize spontaneous reactions in
    addition to catalyzed reactions.

    linear_objective: None or a dictionary of {:class:`~cobra.core.Reaction`: objective_coefficient}


    Returns the minimum flux that can support maximization of the linear objective.

    
    Q: Do we want to denote which reactions are active in the parsimonious flux space?  Possibly
    in a separate function?
    
    """
    irreversible = True

    _initial_model = model
    model = model.copy() #We'll be modifying model structure so we want an independent copy
    #1. Optimize for the initial linear objective
    if linear_objective is None:
        model.optimize()
    else:
        model.optimize(new_objective=linear_objective)


    #Constrain the 
    objective_reactions = [x for x in model.reactions if x.objective_coefficient != 0]
    linear_solution = model.solution.f
    for reaction in objective_reactions:
        reaction.lower_bound = reaction.upper_bound = reaction.objective_coefficient * linear_solution
        reaction.objective_coefficient = 0

    #2. Check whether the model is irreversible or not.
    #Preferable to work with model where all reactions are in one direction
    lower_bounds = [x.lower_bound for x in model.reactions]
    minimum_lower_bound = min(lower_bounds)
    maximum_lower_bound = max(lower_bounds)
    upper_bounds = [x.upper_bound for x in model.reactions]
    minimum_upper_bound = min(upper_bounds)
    maximum_upper_bound = max(upper_bounds)

    if minimum_lower_bound < 0 and maximum_upper_bound > 0:
        irreversible = False
    elif minimium_lower_bound < 0 and maximum_lower_bound > 0:
        irreversible = False
    elif minmium_upper_bound < 0 and maximum_upper_bound > 0:
        irreversible = False

    if not irreversible:
        convert_to_irreversible(model)

    #3. Calculate the minimum flux through the network subject to the initial optimization
    #Create a metabolite to measure the flux in the model
    flux_measure_metabolite = Metabolite('flux_measure')
    #Create a reaction to use to monitor the flux measurement metabolite
    flux_measure_reaction = Reaction('net_flux')
    flux_measure_reaction.add_metabolites({flux_measure_metabolite: -1})
    flux_measure_reaction.objective_coefficient = 1.
    model.add_reaction(flux_measure_reaction)
    
    model.add_metabolites([flux_measure_metabolite])
    if minimize_spontaneous_reactions:
        [x.add_metabolites({flux_measure_metabolite: 1}, add_to_container_model=False)
         for x in model.reactions]
    else:
        [x.add_metabolites({flux_measure_metabolite: 1}, add_to_container_model=False)
         for x in model.reactions if len(x.genes) > 0]
    
    
    model.optimize(objective_sense='minimize')
    return(model.solution.f, linear_solution)
    
    


