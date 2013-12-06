#cobra.core.Reaction.py
#######################
#BEGIN Class Reaction
#
#Is it better to restrict a Reaction to a single model or
#should we allow a Reaction to be associated with multiple models?
#
from collections import defaultdict
import re
from copy import deepcopy
from .Object import Object
from .Metabolite import Metabolite
from .Gene import Gene

from warnings import warn

class Reaction(Object):
    """Reaction is a class for holding information regarding
    a biochemical reaction in a cobra.Model object 

    """
    ## __slots__ = ['id', 'reversibility', '_metabolites', 'gene_reaction_rule',
    ##              'subsystem', '_genes', '_model',
    ##              'name', 'lower_bound', 'upper_bound', 'objective_coefficient',
    ##              'reaction', 'boundary']

    def __init__(self, name=None):
        """An object for housing reactions and associated information
        for cobra modeling..

        boundary: Boolean.  Indicates whether the reaction is at a boundary, i.e.,
        a source or a sink.
        
        """
        Object.__init__(self, name)
        self.gene_reaction_rule = '' #Deprecated
        self.subsystem = ''
        self._genes = set() #The cobra.Genes that are used to catalyze the reaction
        #reaction.  _ Indicates that it is not preferred to add a gene to a reaction
        #directly.
        #A dictionary of metabolites and their stoichiometric coefficients in
        #this reaction  _ Indicates that it is not preferred to add a metabolite to a reaction
        #directly.
        self._metabolites = {}
        self.name = name
        #self.model is None or refers to the cobra.Model that
        #contains self
        self.objective_coefficient = self.lower_bound = 0.
        self.upper_bound = 1000.
        self.reflection = None #Either None or if this reaction is irreversible then
        #a reaction in the model that is essentially self * -1
        self.variable_kind = 'continuous' #Used during optimization.  Indicates whether the
        #variable is modeled as continuous, integer, binary, semicontinous, or semiinteger.


    @property
    def reversibility(self):
        """This property removes the independence of the reversibility attribute and the reaction's
        current upper and lower bounds.

        reversibility is defined in the context of the current instantiation.
        
        """
        return self.lower_bound < 0 and self.upper_bound > 0
    
    @property
    def boundary(self):
        # single metabolite implies it must be a boundary
        if len(self._metabolites) == 1:
            return "system_boundary"
        # if there is more than one metabolite, if it ONLY produces or ONLY
        # consumes, it is also a boundary.
        all_stoichiometry = self._metabolites.values()
        if not min(all_stoichiometry) < 0 < max(all_stoichiometry):
            return "system_boundary"
        return None

    def _update_awareness(self):
        """Make sure all metabolites and genes that are associated with
        this reaction are aware of it.
        
        """
        [x._reaction.add(self) for x in self._metabolites]
        [x._reaction.add(self) for x in self._genes]

        

    def remove_from_model(self, model=None):
        """Removes the association

        model: cobra.Model object.  remove the reaction from this model.
        
        """
        # why is model being taken in as a parameter? This plays
        #back to the question of allowing a Metabolite to be associated
        #with multiple Models
        if model != self._model and model is not None:
            raise Exception('%s not in %s ergo it cannot be removed. (%s)'%(self,
                                                                  model,
                                                                  self._model))
                                                            
        new_metabolites = deepcopy(self._metabolites)
        new_genes = deepcopy(self._genes)
        self._model.reactions.remove(self)
        #Remove associations between the reaction and its container _model
        #and elements in the model
        self._model = None
        [x._reaction.remove(self)
         for x in self._metabolites.keys()]
        [x._reaction.remove(self)
         for x in self._genes]
        #Replace the model-linked metabolites with the new independent metabolites
        self._metabolites = {}
        self.add_metabolites(new_metabolites)
        #Replace the model-linked genes with new indepenent genes
        self._genes = set()
        [self.add_gene(k) for k in new_genes]

    def delete(self):
        """Removes all associations between a reaction and its container
        _model and metabolites and genes.
        
        TODO: Decide whether orphan metabolites should automatically be removed
        from the model.
        
        """
        self._model = None
        [x._reaction.remove(self)
         for x in self._metabolites.keys() if self in x._reaction]
        [x._reaction.remove(self)
         for x in self._genes if self in x._reaction]
        self._metabolites = {}
        self._genes = set()
        
        
    def __setstate__(self, state):
        """Probably not necessary to set _model as the cobra.Model that
        contains self sets the _model attribute for all metabolites and genes in the reaction.

        However, to increase performance speed we do want to let the metabolite and gene
        know that they are employed in this reaction

        """
        if "reaction" in state:
            state.pop("reaction")
        self.__dict__.update(state)
        for x in state['_metabolites']:
            setattr(x, '_model', self._model)
            x._reaction.add(self)
        for x in state['_genes']:
            setattr(x, '_model', self._model)
            x._reaction.add(self)

    
    def copy(self):
        """When copying a reaction, it is necessary to deepcopy the
        components so the list references aren't carried over.

        Additionally, a copy of a reaction is no longer in a cobra.Model.

        This should be fixed with self.__deecopy__ if possible
        """
        ## the_model = self._model
        ## self._model = None
        new_reaction = deepcopy(self)
        ## self._model = the_model
        return new_reaction

    def guided_copy(self, the_model, metabolite_dict, gene_dict=None):
        """Trying to make a faster copy procedure for cases where large
        numbers of metabolites might be copied.  Such as when copying reactions.

       the_model: The new container cobra.Model

       metabolite_dict: A dictionary of the cobra.Metabolite objects that are in
       the model (metabolite.id, metabolite)

        gene_dict: A dictionary of the cobra.Gene objects that are in the model (gene.id, gene)


        """
        the_copy = Object.guided_copy(self)
        #Replace the complex items in a faster fashion
        the_copy._model = the_model
        if gene_dict:
            the_copy._genes = set([gene_dict[k.id]
                                for k in self._genes])
        the_copy._metabolites = dict([(metabolite_dict[k.id], v)
                                      for k, v in self._metabolites.iteritems()])

        #make the metabolites and genes aware of the reaction
        [k._reaction.add(the_copy)
         for k in the_copy._genes]
        [k._reaction.add(the_copy)
         for k in the_copy._metabolites.keys()]

        return(the_copy)

    def pop(self, the_metabolite):
        """Remove a metabolite from the reaction and return the
        stoichiometric coefficient.

        the_metabolite: A cobra.Metabolite that is in the reaction
    
        
        """
        the_coefficient = self._metabolites.pop(the_metabolite)
        the_metabolite._reaction.remove(self)
        return the_coefficient
    
    def __add__(self, other_reaction):
        """Adds two reactions to each other.  Default behavior is
        to combine the metabolites but only use the remaining parameters
        from the first object.

        TODO: Deal with gene association logic from adding reactions.

        TODO: Simplify and add in an __iadd__
        
        """
        new_reaction = deepcopy(self)
        new_reaction.id = self.id + '_' + other_reaction.id
        new_reaction.add_metabolites(deepcopy(other_reaction._metabolites))
        new_reaction._genes.update(deepcopy(other_reaction._genes))
        #Make all the genes aware of this reaction
        [x._reaction.add(new_reaction) for x in new_reaction._genes]
        gpr_1 = new_reaction.gene_reaction_rule
        gpr_2 = other_reaction.gene_reaction_rule
        if gpr_1 != '' and gpr_2 != '':
            new_reaction.gene_reaction_rule = '%s and %s'%(gpr_1, gpr_2)
        elif gpr_2 != '':
            new_reaction.gene_reaction_rule = gpr_2
        return new_reaction
    
    def __sub__(self, other_reaction):
        """Subtracts two reactions.  Default behavior is
        to combine the metabolites but only use the remaining parameters
        from the first object.

        Note: This is equivalent to adding reactions after changing the sign
        of the metabolites in other_reaction
        
        """
        new_reaction = deepcopy(self)
        if self is other_reaction:
            other_reaction = deepcopy(other_reaction)
        new_reaction.id = self.id + '_' + other_reaction.id
        new_reaction.subtract_metabolites(deepcopy(other_reaction._metabolites))
        return new_reaction

    def __imul__(self, the_coefficient):
        """Allows the reaction coefficients to be rapidly scaled.
        
        """
        [self._metabolites.update({k: the_coefficient * v})
         for k, v in self._metabolites.items()]
        return self

    def __mul__(self, the_coefficient):
        """Allows a reaction to be multipled by a coeffient.
        
        Should this return a new reaction?
        
        """
        [self._metabolites.update({k: the_coefficient * v})
         for k, v in self._metabolites.items()]
        return self
        

    def parse_gene_association(self, the_type='gene'):
        """Extract all genes from the Boolean Gene_Association string.

        #Formerly, update_names
        """
        if the_type == 'gene':
            self._genes = set((re.compile(' {2,}').sub(' ', re.compile('\(| and| or|\+|\)').sub('', self.gene_reaction_rule))).split(' ' ))
            if '' in self._genes:
                self._genes.remove('')
            self._genes = set(map(Gene, self._genes))
            #Make the gene aware that it is involved in this reaction
            [x._reaction.add(self) for x in self._genes]


    def add_gene_reaction_rule(self, the_rule):
        """This adds a gene reaction rule.

        the_rule:  A boolean representation of the gene requirements for this reaction
        to be active as described in Schellenberger et al 2011 Nature Protocols 6(9):1290-307.

        Note that this method currently replaces any pre-existing rules
        
        """
        self.gene_reaction_rule = the_rule
        self.parse_gene_association()
        
    def get_reactants(self):
        from warnings import warn
        warn('Reaction.get_reactants is deprecated.  please use Reaction.reactants')
        return(self.reactants)

    def get_products(self):
        from warnings import warn
        warn('Reaction.get_products is deprecated.  please use Reaction.products')
        return(self.products)

    def get_gene(self):
        from warnings import warn
        warn('Reaction.get_gene is deprecated.  please use Reaction.genes')
        return(self.genes)
        

    @property
    def genes(self):
        """Return a list of genes for the reaction.

        """
        return list(self._genes)
        

    @property
    def products(self):
        """Return a list of products for the reaction
        
        """
        return([k for k, v in self._metabolites.items()
                if v > 0])

    @property
    def reactants(self):
        """Return a list of reactants for the reaction.

        """
        return([k for k, v in self._metabolites.items()
                if v < 0])
        
    @property
    def reaction(self):
        return self.build_reaction_string()

    @property
    def species(self):
        if hasattr(self, '_species'):
            return self._species.keys()
        else:
            return self._metabolites.keys()

    def get_coefficient(self, the_metabolite):
        """Return the stoichiometric coefficient for a metabolite in
        the reaction.

        the_metabolite: A metabolite Id.
        
        """
        _id_to_metabolites = dict([(x.id, x)
                                        for x in self._metabolites])

        if hasattr(the_metabolite, 'id'):
            the_metabolite = the_metabolite.id
        return self._metabolites[_id_to_metabolites[the_metabolite]]
    
    def get_coefficients(self, the_metabolites):
        """Return the stoichiometric coefficients for a list of
        metabolites in the reaction.

        the_metabolites:  A list of metabolite Ids.
        
        """
        return map(self.get_coefficient, the_metabolites)
    
    def add_metabolites(self, the_metabolites, combine=True, add_to_container_model=True):
        """Add metabolites and stoichiometric coefficients to the reaction.
        If the final coefficient for a metabolite is 0 then it is removed
        from the reaction.

        the_metabolites: A dict of cobra.Metabolites and their coefficient

        combine: Boolean. If True and a metabolite already exists in the
        reaction then the coefficients will be added.  If False the old
        metabolite will be discarded and the new one added.

        add_to_container_model: Boolean.  If True and this reaction is
        contained within a cobra.Model (i.e., self._model is a cobra.Model)
        then add the metabolite to the model.

        """
        _id_to_metabolites = dict([(x.id, x)
                                        for x in self._metabolites])
        new_metabolites = []
        for the_metabolite, the_coefficient in the_metabolites.items():
            #If a metabolite already exists in the reaction then
            #just add them.
            if the_metabolite.id in _id_to_metabolites:
                reaction_metabolite = _id_to_metabolites[the_metabolite.id]
                if combine:
                    self._metabolites[reaction_metabolite] += the_coefficient
                else:
                    self._metabolites[reaction_metabolite] = the_coefficient
            else:
                self._metabolites[the_metabolite] = the_coefficient
                #make the metabolite aware that it is involved in this reaction
                the_metabolite._reaction.add(self)
                new_metabolites.append(the_metabolite)
        for the_metabolite, the_coefficient in self._metabolites.items():
            if the_coefficient == 0:
                #make the metabolite aware that it no longer participates
                #in this reaction
                the_metabolite._reaction.remove(self)
                self._metabolites.pop(the_metabolite)
        _id_to_metabolites = dict([(x.id, x)
                                        for x in self._metabolites])
        if add_to_container_model and hasattr(self._model, 'add_metabolites'):
            self._model.add_metabolites(new_metabolites)
            

    def subtract_metabolites(self, metabolites):
        """This function will 'subtract' cobra.metabolites from a reaction, which
        means add the metabolites with -1*coefficient.  If the final coefficient
        for a metabolite is 0 then the metabolite is removed from the reaction.

        .. note:: That a final coefficient < 0 implies a reactant.

        metabolites: dict of {:class:`~cobra.core.Metabolite`: coefficient}
            These metabolites will be added to the reaction

        .. note:: This function uses deepcopy in case the reaction is being
        subtracted from itself.
        
        """
        metabolites = deepcopy(metabolites)
        metabolites = dict([(k, -v)
                                for k, v in metabolites.items()])
        self.add_metabolites(metabolites)



    def build_reaction_string(self, use_metabolite_names=False):
        """Generate a human readable reaction string.
        
        """
        _products = self.products
        _reactants = self.reactants
        _products.sort()
        _reactants.sort()
        reaction_string = ''
        for _species in _reactants:
            coefficient = repr(abs(self._metabolites[_species]))
            if use_metabolite_names:
                id = _species.name
            else:
                id = _species.id
            reaction_string += ' + %s %s'%(coefficient, id)

        if self.reversibility:
            reaction_string += ' <=> '
        else:
            if self.lower_bound < 0 and self.upper_bound <=0:
                reaction_string += ' <- '
            else:
                reaction_string += ' -> '                
        for _species in _products:
            coefficient = repr(self._metabolites[_species])
            if use_metabolite_names:
                id = _species.name
            else:
                id = _species.id
            reaction_string += ' + %s %s'%(coefficient, id)



        reaction_string = reaction_string.lstrip(' + ').rstrip(' + ')
        return(reaction_string)


    def check_mass_balance(self):
        """Makes sure that the reaction is elementally-balanced.

        """
        reaction_element_dict = defaultdict(list)
        for the_metabolite, the_coefficient in self._metabolites.items():
            if the_metabolite.formula is not None:
                [reaction_element_dict[k].append(the_coefficient*v)
                 for k, v in the_metabolite.formula.elements.items()]
        reaction_element_dict = dict([(k, sum(v))
                                      for k, v in reaction_element_dict.items()])
        if sum(map(abs, reaction_element_dict.values())) != 0:
            return [self.id, reaction_element_dict]
        else:
            return []
        
    def print_values(self):
        """Prints most of the contents of a reaction as a series of strings.
        
        """
        
        print "reaction:", self.id
        print "subsystem", self.subsystem
        print self.reaction
        print "metabolites:"
        print self._metabolites
        print "bounds: (%.2f, %.2f)" % (self.lower_bound, self.upper_bound)
        print "objective_coefficient", self.objective_coefficient
        print "gene reaction rule:", self.gene_reaction_rule
        print "all genes:", ", ".join(i.id for i in self._genes)


    def get_compartments(self):
        """
        """
        return(list(set([x.compartment for x in self._metabolites])))


    def remove_gene(self, cobra_gene):
        """Removes the association between a gene and a reaction

        cobra_gene: :class:`~cobra.core.Gene`. A gene that is associated with the reaction.
        
        """
        try:
            self._genes.remove(cobra_gene)
            cobra_gene._reaction.remove(self)
        except Exception, e:
            try:
                if hasattr(self._genes, 'keys'):
                    self._genes = set(self._genes.keys())
                    self.remove_gene(cobra_gene)
            except:
                raise Exception('Unable to remove gene %s from reaction %s: %s'%(cobra_gene.id, self.id, e))

    def add_gene(self, cobra_gene):
        """Associates a cobra.Gene object with a cobra.Reaction.

        cobra_gene: :class:`~cobra.core.Gene`. A gene to associate with the reaction.
        """
        try:
            self._genes.add(cobra_gene)
            cobra_gene._reaction.add(self)
            cobra_gene._model = self._model
        except Exception, e:
            try:
                if hasattr(self._genes, 'keys'):
                    self._genes = set(self._genes.keys())
                    self.add_gene(cobra_gene)
            except:
                raise Exception('Unable to add gene %s to reaction %s: %s'%(cobra_gene.id, self.id, e))
                            
