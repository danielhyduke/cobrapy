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


    #Constrain the variables in the linear objective to the 'optimal' value.
    objective_reactions = [x for x in model.reactions if x.objective_coefficient != 0]
    if len(objective_reactions) > 1:
        raise(Exception('parsimonious FBA only works with a single non-zero objective_coefficient'))
    #To run parsimonious FBA with a linear objective with multiple non-zero objective_coefficients
    #it would be necessary to perform sampling to 'account' for the possibility of multiple
    #equivalent optima.
    linear_solution = model.solution.f
    for reaction in objective_reactions:
        reaction_target_flux = linear_solution / reaction.objective_coefficient
        #Potential bug
        reaction.lower_bound = reaction_target_flux
        #reaction.upper_bound = reaction_target_flux
        
        reaction.objective_coefficient = 0

    #2. Convert model to irreversible form where all variable lower bounds are > 0.
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
    
    
    model.optimize(objective_sense='minimize', new_objective=flux_measure_reaction)
    from pdb import set_trace
    set_trace()
    return(model.solution.f, linear_solution)
    
    


