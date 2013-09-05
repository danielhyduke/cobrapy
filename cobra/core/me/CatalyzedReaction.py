from cobra import Reaction, DictList
import re
from warnings import warn
from .Enzyme import Catalyst, Complex, Subunit
from .tools import split_logical_string
class CatalyzedReaction(Reaction):
    """Need to deal with complexes that have modifications.

    
    """
    def __init__(self, name=None):
        """

        update_container_model_references: Boolean.  If name is a cobra.Reaction and True
        then update the name._model.reaction to refer to the new catalyzed reaction.


        WARNING: If name is a cobra.Reaction and name._model is not None then name._model.reactions will
        be updated with CatalyzedReaction(name) and references to the cobra.Reaction will be removed by
        calling name.delete().
        
        """
        
        warn("CatalyzedReaction does not yet deal with modifications to " +\
             "Complexes b/c this information is not included in the gene " +\
             "reaction rules.")
        ## if isinstance(name, CatalyzedReaction):
        ##     warn("CatalyzedReactions cannot be converted into catalyzed reactions")
        ##     return name

        self._catalyst = set()
        self._complex = set()
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
            self.logic = self.gene_reaction_rule
            if self.logic in ['', ' ', None]:
                warn("%s does not have an associated gene logic string"%self.id)
            else:
                if self.logic.startswith('('):
                    #Deal with annoying gpr syntax
                    if self.logic.count('(')  == self.logic.count(')'):
                        if self.logic.endswith(')'):
                            self.logic = self.logic[1:-1].strip()
                        
            self.__create_catalytic_reaction()
                #Assume that no modifications are necessary for the catalysts to be active
                #NOTE: This introduces a potential bug where a catalyst of a fixed name can
                #catalyze multiple reactions
                #NOTE: This fails
                #self._catalyst = set([x.create_catalyst() for x in self._complex])

    def __setstate__(self, state):
        """
        """
        Reaction.__setstate__(self, state)
        #Possibly unnecessary as genes are currently the only form
        #of catalyst and their reaction set is updated in the call
        #to the parent Reaction.__setstate__
        for the_complex in state['_complex']:
            the_complex._reaction.add(self)
            for the_subunit in the_complex._subunits.keys():
                the_subunit._reaction.add(self)
        [k._reaction.add(self) for k in state['_catalyst']]


        

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
            self._complex.add(tmp_complex)
            tmp_complex._reaction.add(self)


    def guided_copy(self, model, metabolite_dict=None, subunit_dict=None,
                    modification_dict=None, complex_dict=None):
        """Should the complexes be guided_copy here?

        model: The new container cobra.Model

        metabolite_dict: A dictionary of the cobra.Metabolite objects that are in the model (Metabolite.id, Metabolite)

        subunit_dict: A dictionary of the cobra.me.Enzyme.Subunit objects that are in the model (Subunit.id, Subunit)

        complex_dict: A dictionary of the cobra.me.Enzyme.Complex objects that are in the model (Complex.id, Complex)

         modification_dict:  A dictionary of the cobra.me.Enzyme.Modification objects in the model (Modification.id, Modification).

         Q: How do we deal with a Complex of a specific composition mapping to multiple catalysts?
         
        """
        the_copy = Reaction.guided_copy(self, model, metabolite_dict, subunit_dict)

        #Associate ME_Model versions of the Objects with a copy of the CatalyzedReaction
        the_copy._complex = set()
        the_copy._genes = set()
        for the_complex in self._complex:
            model_complex = complex_dict[the_complex.id]
            model_complex._reaction.add(the_copy)
            the_copy._complex.add(model_complex)
            #Update awareness for the catalysts
            for the_catalyst in model_complex._catalysts:
                the_catalyst._reaction.add(the_copy)
                the_copy._catalyst.add(the_catalyst)
            ## for the_modification in the_complex._modifications:
            ##     try:
            ##         the_catalyst = model_complex._modification_to_catalyst[the_modification.id]
            ##     except:
            ##         the_catalyst = model_complex.create_catalyst(modification_dict[the_modification.id])
            ##     the_copy._catalyst.add(the_catalyst)
            for the_gene in model_complex._subunits:
                the_gene._reaction.add(the_copy)
                the_copy._genes.add(the_gene)
        return the_copy


