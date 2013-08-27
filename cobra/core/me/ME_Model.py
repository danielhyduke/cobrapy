#The basic question is whether we're going to eventually deprecate
#the Reaction object and GPRs in favor of CatalyzedReactions
from copy import deepcopy
from cobra import Model, DictList, Object
from cobra.core.me.CatalyzedReaction import CatalyzedReaction
from cobra.core.me.Enzyme import Modification
#from .CatalyzedReaction import CatalyzedReaction
class ME_Model(Model):
    """An ME_Model is a cobra.Model that provides support for
    explicit representation of gene products.

    Model: A string or a cobra.Model.  If a cobra.Model then the underlying structures
    are copied.  NOTE: It might be faster to use the guided copy which will then destroy
    the integrity of the Model

    """

        
    def __init__(self, id):
        if isinstance(id, Model):
            id = deepcopy(id)
            self._copy_parent_attributes(id)
            #Update the references for the Model.
            [[setattr(x, '_model', self) for x in y]
             for y in [self.reactions, self.genes, self.metabolites]]
        else:
            Model.__init__(self, id)
        self.subunits = DictList()
        self.complexes = DictList()
        self.catalysts = DictList()
        self.modifications = DictList()
        #If the id is a cobra.Model then the next line will convert all of
        #the Reactions to CatalyzedReactions
        [CatalyzedReaction(k)
         for k in self.reactions]
        
        self.genes = self.subunits

    def __setstate__(self, state):
        """Make sure all cobra.Objects in the model point to the model
        
        TODO: Make sure that the genes and metabolites referenced by
        the reactions.(genes|metabolites) point to the model's genes
        and metabolites.
        
        """
        self.__dict__.update(state)
        [[setattr(x, '_model', self)
          for x in self.__dict__[y]]
         for y in ['complexes', 'subunits', 'metabolites', 'reactions']]

        
    def copy(self):
        """Provides a partial 'deepcopy' of the ME_Model.  All of the Metabolite, Gene,
        Reaction, Subunit, Complex, and Catalyst objects are created anew but in a faster fashion than deepcopy

        #TODO: Define guided copy for Complexes, Catalysts, and Modifications

        #Calling a Model.copy before copying the complexes causes problems
        because of how Reactions are Catalyzed

        """
        from pdb import set_trace
        #Create a new ME_Model instance
        the_copy = Object.guided_copy(self)
        the_copy.compartments = deepcopy(self.compartments)
        #Create new cobra.Metabolites associated with the_copy ME_Model
        the_metabolites = DictList([x.guided_copy(the_copy)
                                    for x in self.metabolites])
        #Create new cobra.me.Enzyme.Subunits associated with the_copy ME_Model
        the_subunits = DictList([x.guided_copy(the_copy)
                                 for x in self.subunits])
        #Create new cobra.me.Enzyme.Modifications associated with the_copy ME_Model
        the_modifications = DictList([x.guided_copy(the_metabolites._object_dict)
                                 for x in self.modifications])
        #Add in an empty Modification that allows Complexes to become Catalysts without
        #any actual modifications
        _empty_modification = Modification()
        the_modifications.append(_empty_modification)
        #Create new cobra.me.Enzyme.Complexes associated with the_copy ME_Model
        the_complexes = DictList([x.guided_copy(the_copy, the_modifications._object_dict, the_subunits._object_dict)
                                 for x in self.complexes])

        #Create catalysts for all complexes.
        the_catalysts = DictList()
        for the_complex in the_complexes:
            if len(the_complex._modifications) == 0:
                #assume that the complex doesn't need to be modified to be activated
                #and set the stoichiometry for the empty modification to 0
                the_catalysts.append(the_complex.create_catalyst(_empty_modification, 0))
            else:
                #Create catalysts for each modification.
                #TODO: Deal with stoichiometry for modifications
                for the_modification in the_complex._modifications:
                    the_catalysts.append(the_complex.create_catalyst(the_modification, 1))

        #Create new cobra.me.CatalyzedReactions associated with the_copy ME_Model
        the_reactions = DictList([x.guided_copy(the_copy,
                                                the_metabolites._object_dict,
                                                the_subunits._object_dict,
                                                the_modifications._object_dict,
                                                the_complexes._object_dict)
                                  for x in self.reactions])

        the_copy.genes = the_copy.subunits = the_subunits
        the_copy.modifications = the_modifications
        the_copy.catalysts = the_catalysts
        the_copy.complexes = the_complexes
        the_copy.reactions = the_reactions
        return the_copy
