from cobra import Reaction, DictList
import re
from warnings import warn
from .Enzyme import Catalyst, Complex, Subunit
from .tools import split_logical_string
class CatalyzedReaction(Reaction):
    """Need to deal with complexes that have modifications.

    
    """
    def __init__(self, name=None, _create_catalyst_from_gene_association=True):
        """

        update_container_model_references: Boolean.  If name is a cobra.Reaction and True
        then update the name._model.reaction to refer to the new catalyzed reaction.


        _create_catalyst_from_gene_association:  Boolean.  If True then uses logic (aka gene reaction rule)
        to construct a catalyst.

        WARNING: If name is a cobra.Reaction and name._model is not None then name._model.reactions will
        be updated with CatalyzedReaction(name) and references to the cobra.Reaction will be removed by
        calling name.delete().
        
        """
        self._catalyst = set()
        #self._complex = set() removed complex because a catalyst catalyzes the reaction.
        self._genes = set() #Replacing the stoichiometry dict from superclass reaction with
        #a set because a single Gene may be in multiple Complexes, that catalyze the
        #CatalyzedReaction, with different Stoichiometries
        self.logic = None

        if not isinstance(name, Reaction):
            Reaction.__init__(self, name)
        else:
            self._copy_parent_attributes(name)
            #Clear the _metabolites attribute so that
            #they can be added and referential integrity can be
            #maintained.
            self._metabolites = {}
            self.add_metabolites(name._metabolites)
            self._species = self._metabolites

            #Deal with pre-0.3.0 Reactions where the _genes attribute was a dict
            if hasattr(name._genes, 'keys'):
                self._genes = set(name._genes.keys())
            [x._reaction.add(self) for x in self._genes]
            self.logic = self.gene_reaction_rule #Deprecated
            if self.logic in ['', ' ', None]:
                warn("%s does not have an associated gene logic string"%self.id)
            else:
                if self.logic.startswith('('):
                    #Deal with annoying gpr syntax
                    if self.logic.count('(')  == self.logic.count(')'):
                        if self.logic.endswith(')'):
                            self.logic = self.logic[1:-1].strip()

            if _create_catalyst_from_gene_association:
                self.__create_catalytic_reaction()
                #Assume that no modifications are necessary for the catalysts to be active

    @property
    def catalysts(self):
        return(list(self._catalyst))

    @property
    def complexes(self):
        complexes = set()
        [complexes.add(x.complex) for x in self.catalysts]
        return(list(complexes))

    def __setstate__(self, state):
        """
        """
        Reaction.__setstate__(self, state)
        #Possibly unnecessary as genes are currently the only form
        #of catalyst and their reaction set is updated in the call
        #to the parent Reaction.__setstate__
        for the_catalyst in state['_catalyst']:
            the_catalyst._reaction.add(self)
            for the_subunit in the_catalyst.complex.subunits:
                the_subunit._reaction.add(self)


        

    def __create_catalytic_reaction(self):
        """create_catalytic_reaction uses a cobra.Reaction to populate
        the CatalyzedReaction.  If the reaction has a model associated with
        it then the reaction updates the model.

        NOTE: We assume that the cobra.Reaction.gene_reaction_rule is in
        disjunctive normal form.

        NOTE: We have to deal with modifications required for the reaction.
        

        NOTE: If the reaction is associated with a model this function does
        not traverse the Genes awareness completely, so it might be best to just
        convert reaction._model.genes to subunits 
        
        """
        #Makes sure to account for changing the Gene to a Subunit in the
        #upstream portions. reaction._genes and reaction._model.genes
        #
        if self._model is not None:
            the_model = self._model
            #Clear references for the previous version of the reaction
            _old_reaction = the_model.reactions.get_by_id(self.id)
            _old_reaction.delete()
            #Update the Reaction in the Model to point to the CatalyzedReaction
            the_model.reactions._replace_on_id(self)
            if not hasattr(the_model, 'subunits'):
                the_model.subunits = DictList() #Potential for BUG
                #The subunits attribute is added to deal with cases where a gene
                #can encode multiple active subunits (i.e., RNA and protein)
                the_model.complexes = DictList()
                the_model.catalysts = DictList()
        else:
            the_model = None
        the_genes = self._genes
        if len(the_genes) == 0:
            warn("%s does not have associated genes so generating a CatalyzedReaction sans genes"%self.id)
            return()
        #_id_to_stoichiometry = dict([(k.id, v) for k, v in the_genes.iteritems()])
        #Replace gene objects with Subunits.  If a Model is provided
        self._genes = the_subunits = set()
        for the_gene in the_genes:
            if not isinstance(the_gene, Subunit):
                the_subunit = Subunit(the_gene)
            if the_model:
                try:
                    the_subunit = the_model.subunits.get_by_id(the_subunit.id)
                except:
                    the_model.subunits.append(the_subunit)
            #Associate the CatalyzedReaction with the Subunit
            the_subunit._reaction.add(self)
            the_subunits.add(the_subunit)
        #Dictionary to speed up translations
        subunit_id_to_object = dict([(k.id, k) for k in the_subunits])

        #Parse out logic.  Complexes, and Catalysts
        or_re = re.compile(' +or +', re.I)
        the_complexes = or_re.split(self.logic)
        for the_composition in the_complexes:
            the_composition.strip()
            if the_composition.startswith('('):
                #If the string starts with a parenthesis then assume that the
                #last character is a parenthesis that must me removed.
                if the_composition[-1] != ')':
                    the_composition = the_composition[1:].strip()
                else:
                    the_composition = the_composition[1:-1].strip()
            elif the_composition.count(')') > the_composition.count('(') and \
                     the_composition.endswith(')'):
                the_composition = the_composition[:-1].strip()
            
            
            the_composition = split_logical_string(the_composition)
            tmp_complex = Complex('')
            if the_model is None:
                [tmp_complex.add_subunit(subunit_id_to_object[k], v)
                 for k, v in the_composition.iteritems()]
            else:
                #We don't update the subunit awareness here because an equivalent
                #complex might already be represented in the_model
                [tmp_complex.add_subunit(subunit_id_to_object[k], v, False)
                 for k, v in the_composition.iteritems()]
                tmp_complex._model = the_model
                #TODO: Add an add_complex and add_catalyst function to ME_Model.
                ## from pdb import set_trace
                ## if tmp_complex.id == 'Complex__1_STM0999':
                ##     set_trace()
                try:
                    tmp_complex = the_model.complexes.get_by_id(tmp_complex.id)
                except:
                    the_model.complexes.append(tmp_complex)
                    [x._complex.add(tmp_complex) for x in tmp_complex._subunits.iterkeys()]
            tmp_catalyst = tmp_complex.create_catalyst()
            self.add_catalyst(tmp_catalyst)


    def add_catalyst(self, catalyst):
        """associates a catalyst with the reaction.
        
        catalyst: :class:`~cobra.core.me.Enzyme.Catalyst`

        TODO: Add in an option to check whether the subunits
        exist in the model if self._model is not None.
        
        """
        self._catalyst.add(catalyst)
        catalyst._reaction.add(self)
        map(self.add_gene, catalyst.complex.subunits)
        if self._model is not None:
            if self._model != catalyst._model:
                catalyst._model = self._model
                self._model.catalysts.append(catalyst)

    def remove_catalyst(self, catalyst):
        """
        TODO: Add in a parameter to remove the CatalyzedReaction from the model if no Catalysts are
        associated with the CatalyzedReaction after deleting the Catalyst
        
        """
        self._catalyst.remove(catalyst)
        #Remove the genes from the reaction that are not associated with any
        #remaining catalysts for the reaction.
        _active_subunits = set()
        [_active_subunits.update(k.complex.subunits) for k in self.catalysts]
        _inactive_subunits = set(catalyst.complex.subunits).difference(_active_subunits)
        map(self.remove_gene, _inactive_subunits)

        
        
    def guided_copy(self, model, metabolite_dict=None, subunit_dict=None,
                    modification_dict=None, catalyst_dict=None):
        """Should the complexes be guided_copy here?

        model: The new container cobra.Model

        metabolite_dict: A dictionary of the cobra.Metabolite objects that are in the model (Metabolite.id, Metabolite)

        subunit_dict: A dictionary of the cobra.me.Enzyme.Subunit objects that are in the model (Subunit.id, Subunit)

        catalyst_dict: A dictionary of the cobra.me.Enzyme.Catalyst objects that are in the model (Catalyst.id, Catalyst)

         modification_dict:  A dictionary of the cobra.me.Enzyme.Modification objects in the model (Modification.id, Modification).

         Q: How do we deal with a Complex of a specific composition mapping to multiple catalysts?
         
        """
        the_copy = Reaction.guided_copy(self, model, metabolite_dict, subunit_dict)
        warn("CatalyzedReaction.guided_copy has not yet been thoroughly auditted")
        #Associate ME_Model versions of the Objects with a copy of the CatalyzedReaction
        the_copy._catalyst = set()
        the_copy._genes = set()
        for the_catalyst in self._catalyst:
            model_catalyst = catalyst_dict[the_catalyst.id]
            model_catalyst._reaction.add(the_copy)
            the_copy._catalyst.add(model_catalyst)
            for the_gene in model_catalyst.complex.subunits:
                the_gene._reaction.add(the_copy)
                the_copy._genes.add(the_gene)
        return the_copy



    def add_to_model(self, model):
        """
        model: cobra.core.me.ME_Model
        
        """
        ## from pdb import set_trace
        ## set_trace()
        if self.model is not model:
            if self.model is not None:
                raise(Exception('CatalyzedReaction (%s) already associated with model (%s) cannot add to %s.'%(self.id, self.model.id,
                                                                                                  model.id)  +\
                                'Have not implemented code to change association to a different model'))
            self._model = model
            model.add_reaction(self)
            [x.add_to_model(model) for x in self._catalyst];


        
