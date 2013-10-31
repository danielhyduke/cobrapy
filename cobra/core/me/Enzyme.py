#Match the enzymes from the S. LT2 and E. K12 models
from cobra import Gene, Species, Object
from warnings import warn
#TODO: Generalize the updating of ids after adding in elements.
class Modifier(Species):
    """A Modifier is a Species that modifies Complexes into Catalysts.
    
    The class is primarily used to track the relationship between the Modifier and
    associated modifications.  While not necessary for simple modifications, it will
    be useful for complex Modifications where n Modifications may use 1 species.
    
    """
    def __init__(self, id):
        raise Exception('not implemented yet')
        if isinstance(id, Species):
            self = id
            self._modifications = set()
            #update type
            #return(self)
        else:
            Species.__init__(self, id)
            self._modifications = set()
            
class Modification(Object):
    """This is a rough way for modifying enzymes.  Eventually, it will be better to build
    a reaction class for each type of modification and then use them to create Catalysts
    from Subunits and Modifications.
    
    """
    def __init__(self, id='Empty_Modification'):
        self.id = id
        self._modification_dict = {} #A dictionary that holds the modifications as the keys
        #and the stoichiometries as the values.
        self.logic = None
        self._targets = set() #Used to account for which Objects are modified by a modifier
        #This is helpful when running the remove_from model function.

    @property
    def modifiers(self):
        return(self._modification_dict.keys())

    @property
    def targets(self):
        return list(self._targets)

    def remove_target(self, target):
        """
        """
        self._targets.remove(target)
        target._modifications.pop(self)
        
    @property
    def composition(self):
        return(self._modification_dict)
    
    def add_species(self, species, stoichiometry):
        """Species must be a cobra model species

        stoichiometry: Integer.  The stoichiometry for the Subunit in the enzyme.
        
        """
        if isinstance(species, Species):
            self._modification_dict[species] = stoichiometry
        else:
            raise Exception("You must use cobra.Species when adding a species." +\
                            "You tried to add %s which is %s"%repr(species, type(species)))
        self.__update_id(self._modification_dict)
    def __update_id(self, stoichiometry_dict):
        """Ids are derived from their basic composition following this structure:
        modification__v_k.id where v is the stoichiometric coefficient and k an Object and
        each modification is separated by __.
        
        """
        if len(stoichiometry_dict) == 0:
            self.id = 'Empty_Modification'
        else:
            self.id = 'Modification'
        _tmp_list = stoichiometry_dict.items()
        _tmp_list.sort()
        for k, v in _tmp_list:
            self.id += '__%s_%s'%(repr(v), k.id)
        _update_logic(self, stoichiometry_dict)
        

    def guided_copy(self, species_dict):
        """Speed up copying.
        
        """
        warn("Modification guided copy does not play well with non-species modifications")
        the_copy = Object.guided_copy(self)
        the_copy._modification_dict = dict([(species_dict[k.id], v)
                                            for k, v in self._modification_dict.iteritems()])
        if len(the_copy._modification_dict) > 0:
            the_copy.__update_id(self._modification_dict)
        return the_copy



    def add_to_model(self, model):
        """
        model: cobra.core.me.ME_Model
        
        """
        if self not in model.modifications:
            model.add_metabolites(self.modifiers)
            [x.add_to_model(model) for x in self.targets];
            model.modifications.append(self)

#Blank modification that allows a Complex to become a Catalyst without a Modification
_empty_modification = Modification()

class Complex(Species):
    """
    TODO: Add in method for batch addition of subunits

    """

    def __init__(self, id=''):
        Species.__init__(self,id)
        del self._reaction #A complex isn't directly associated with reactions.
        self._subunits = {}
        #The next two attributes will speed up operations that are only interested in
        #the subset of subunits that falls into the RNA or Polypeptide class
        self._rna_subunits = set() #set which indicates which subunits are RNAs
        self._polypeptide_subunits = set() #set which indicates which subunits are polypeptides
        self._modifications = set() #The set of modifications which may be applied to this complex
        #to render it functional.  If None is in the set then the Complex can perform an operation
        #without any modifications.
        self._catalysts = set() #The set of Catalyst objects that are derived from the specific Complex
        self._subcomplexes = set() #Track complexes that were used to assemble a super-complex
        self._supercomplexes = set() #Track super-complexes that contain this complex
        self.pids = set() #set of protein ids associated with a complex
        self.logic = None

    @property
    def subunits(self):
        return self._subunits.keys()

    @property
    def catalysts(self):
        return list(self._catalysts)

    @property
    def modifications(self):
        return list(self._modifications)

    @property
    def reactions(self):
        _reactions = set()
        [_reactions.update(k.reactions) for k in self.catalysts]
        return(list(_reactions))
        
    def add_to_model(self, model):
        """
        model: cobra.core.me.ME_Model
        
        """
        if self.model is not model:
            if self.model is not None:
                raise(Exception('(%s) already associated with model (%s) cannot add to %s.'%(self.id, self.model.id,
                                                                                             model.id)  +\
                                'Have not implemented code to change association to a different model'))
            self._model = model
            model.complexes.append(self)
            [x.add_to_model(model) for x in self._subcomplexes];
            [x.add_to_model(model) for x in self._supercomplexes];
            [x.add_to_model(model) for x in self._subunits];
            [x.add_to_model(model) for x in self._catalysts];
            
            

    def add_subunit(self, subunit, stoichiometry=1, update_subunit_awareness=True):
        """
        subunit:  A Subunit object.  If a Complex then add the Subunits of the Complex instead
        of the Complex

        stoichiometry: Integer.  The stoichiometry for the Subunit in the enzyme.

        
        
        """
        #Deal with the case where Complexes are combined into a super complex
        if hasattr(subunit, '_subunits'):
            warn('You have added sub-Complex (%s) to Complex (%s).  We will extract the '%(subunit.id, self.name) +\
                 'Subunits from %s and add them to %s'%(subunit.id, self.name))
            [self.add_subunit(k, stoichiometry * v) for k, v in subunit._subunits.iteritems()]
            self._subcomplexes.add(subunit)
            subunit._supercomplexes.add(self)

        else:
            if subunit in self._subunits:
                warn('Updating %s stoichiometry from %s to %s in Complex %s'%(subunit.id,
                                                                              repr(stoichiometry),
                                                                              repr(stoichiometry + self._subunits[subunit]),
                                                                              self.id))
            else:
                self._subunits[subunit] = stoichiometry
                if update_subunit_awareness:
                    subunit._complex.add(self) #Make the subunit aware of the complexes in which it is used
                if subunit.is_polypeptide:
                    self._polypeptide_subunits.add(subunit)
                else:
                    self._rna_subunits.add(subunit)
            self.pids.add(subunit.pid)
            self._update_id(self._subunits)

    def update_subunit_stoichiometry(self, subunit_stoichiometry_dict, _update_model=True):
        """Function that lets one update the stoichiometry for a Subunit or
        Subunits in a list.

        subunit_stoichiometry_dict: A dictionary where the k is the Subunit and
        the value is the new stoichiometry.

        _update_model: Updates model.catalyst upon change of the catalyst id.
        
        """
        [self._subunits.update({k: v})
         for k, v in subunit_stoichiometry_dict.iteritems()]
        for the_catalyst in self._catalysts:
            _old_id = the_catalyst.id
            the_catalyst._update_id(subunit_stoichiometry_dict)
            the_catalyst._update_id(the_catalyst._modifications, prefix=the_catalyst.id)
            
            if self.model is not None and _old_id != the_catalyst.id and _update_model:
                _object_dict = self.model.catalysts._object_dict
                _object_dict.pop(_old_id)
                _object_dict[the_catalyst.id] = the_catalyst
                #self.model.catalysts._generate_index()
        self._update_id(self._subunits)
        
                                                                          
    def _update_id(self, stoichiometry_dict, prefix='Complex'):
        """Ids are derived from their basic composition following this structure:
        modification__v_k.name where v is the stoichiometric coefficient and k an Object and
        each modification is separated by __.

        """
        if self.model is not None:
            #Remove from the complex from the DictList because
            #the id might be changing due to the current convention of
            #constructing the id based on the complex composition.
            self.model.complexes._object_dict.pop(self.id)
            
        self.id = prefix
        _tmp_list = stoichiometry_dict.items()
        _tmp_list.sort()
        for k, v in _tmp_list:
            self.id += '__%s_%s'%(repr(v), k.id)
        _update_logic(self, stoichiometry_dict)
        if self.model is not None:
            #Add the complex back in with its new id
            self.model.complexes._object_dict[self.id] = self

        


    def create_catalyst(self, modification=_empty_modification, stoichiometry=1, data_source=None):
        """Modifies the complex according to the modification
        and creates a specific enzyme.

        modification: A Modification object or None.
        
        """
        if isinstance(modification, list) or isinstance(modification, tuple) or isinstance(modification, set):
            if len(modification) > 1:
                raise Exception("Complex.create_catalyst only supports the use of a single modification.  " +\
                                "A dirty way to do things at the moment is to consider merging modifications into a larger one. " +\
                                "Eventually, we'll add in support for multi modifications")
            else:
                modification = list(modification)[0]

            
        if modification in self.modifications:
            modification = [x for x in self.modifications if x.id == modification.id][0]
            from warnings import warn
            warn("%s has already been applied to complex %s"%(repr(modification),
                                                              self.id) +\
                 " returning the catalyst ")

            catalyst_id = list(set([x.id for x in modification.targets]).intersection([y.id for y in self.catalysts]))
          
            if len(catalyst_id) > 1:
                raise Exception("%s modifies complex %s into multiple catalysts %s"%(modification.id,
                                                                                     self.id,
                                                                                     repr(list(catalyst))))
            else:
                if len(catalyst_id) == 0:
                    from pdb import set_trace
                    set_trace()
                catalyst_id = catalyst_id[0]
                catalyst = [x for x in self.catalysts if x.id == catalyst_id]
                return(catalyst.pop())

                                                              
        catalyst = Catalyst(self)
        if modification is not None:
            if self.model is not None and hasattr(self.model, 'modifications'):
                try:
                    modification = self.model.modifications.get_by_id(modification.id)
                except KeyError:
                    self.model.modifications.append(modification)
            catalyst.modify(modification)
        if self.model is not None and hasattr(self.model, 'catalysts'):
            try:
                catalyst = self.model.catalysts.get_by_id(catalyst.id)
            except:
                self.model.catalysts.append(catalyst)
                catalyst._model = self.model
        if catalyst.complex is not None and catalyst.complex != self:
            raise Exception("Trying to add a catalyst (%s) for %s, but it's already associated with %s in %s"%(catalyst.id,
                                                                                                               self.id,
                                                                                                               catalyst.complex.id,
                                                                                                               self.model))
        else:       
            catalyst._complex = self  #Add in the reference to the complex
        if catalyst in self.catalysts:
            from warnings import warn
            warn("%s already activated for %s.  returning the already actived form"%(catalyst.id,
                                                                                     self.id))
            catalyst = [x for x in self.catalysts if x.id == catalyst.id][0]
            return(catalyst)

        self._catalysts.add(catalyst)
        self._modifications.add(modification)
        if data_source is not None:
            catalyst.source = data_source
        return(catalyst)


    ## def __getstate__(self):
    ##     """Is there anything that needs to be reset with doing a get state for
    ##     a complexe, such as catalyst or modification?
    ##     """
    ##     state = Species.__get__state(self)
    def __setstate__(self, state):
        """Probably not necessary to set _model as the cobra.Model that
        contains self sets the _model attribute for all associated subunits and catalysts


        """
        self.__dict__.update(state)
        for x in state['_subunits']:
            setattr(x, '_model', self._model)
            x._complex.add(self)
        for x in state['_catalysts']:
            setattr(x, '_model', self._model)
            x._complex = self

    def guided_copy(self, model, modification_dict, subunit_dict):
        """Note that this guided copy resets the reaction set to empty for the complex.
        model: The new cobra.Model container object

        modification_dict:  A dictionary of the cobra.me.Enzyme.Modification objects in the model (Modification.id, Modification).

subunit_dict: A dictionary of the cobra.me.Enzyme.Subunit objects that are in the model (Subunit.id, Subunit)        
        
        """
        ## from warnings import warn
        ## warn("Need to update Complex.guided_copy to account for modifications, catalysts, and subunits")
        the_copy = Species.guided_copy(self, model)
        #Associate the appropriate Subunits and stoichiometries with the_copy
        the_copy._subunits = {}
        [the_copy.add_subunit(subunit_dict[k.id], v)
            for k, v in self._subunits.iteritems()];

        the_copy._modifications = set([modification_dict[k.id]
                                        for k in self._modifications])
        
        #It's probably faster to regenerate the catalysts than copy them
        if len(the_copy._modifications) > 0:
            [the_copy.create_catalyst(k) for k in self._modifications]
        return the_copy
        
    
    def delete_catalyst(self, catalyst):
        """
        
        """
        catalyst.delete()
                    
    def remove_from_model(self):
        """

        """
        warn("%s.remove_from_model is not yet implemented"%(str(type(self))))

        
#class Catalyst(Complex): 
#pickling problems appear to arise from the Catalyst being a Complex
class Catalyst(Species):
    """An Catalyst is an active form of a Complex.  The Complex may or
    may not require Modifications to become a Catalyst.

    NOTE: Should probably have different kinds of Catalysts (i.e.,
    Enzyme, Transporter, ...)
    
    """
    def __init__(self, id):
        """

        id: String or :class:`cobra.me.Complex`

        """
        self._model = None
        self._complex = None
        self._reaction = set()
        if isinstance(id, Complex):
            #self._copy_parent_attributes(id)
            self.id = 'Catalyst' + id.id.lstrip('Complex')
            self._complex = id
            self._model = self._complex.model
        elif isinstance(id, str):
            Species.__init__(self, id)
            #Complex.__init__(self, id)
            #self._complex = None #Changed to deal with deepcopy / pickle issues
            #self._complex = self
        self._modifications = {} #The set of modifications that are required
        #to make the enzyme functional

    @property
    def modifications(self):
        return(self._modifications.keys())

    @property
    def reactions(self):
        return(list(self._reaction))

    @property
    def complex(self):
        return(self._complex)

    def add_reaction(self, reaction):
        """

        
        """
        if reaction.model is not None:
            if catalyst.model is not None and catalyst.model is not reaction.model:
                raise(Exception('Cannot add %s to %s because they are associated with '%(reaction.id,
                                                                                         self.id) + \
                                'different reactions %s and %s, respectively'%(reaction.model.id,
                                                                               self.model.id)))
            self.add_to_model(reaction.model)
            self._reaction.add(reaction)
            reaction.catalysts.add(self)

        
    def modify(self, modification, stoichiometry=1):
        """
        modification: A Modification object.

        stoichiometry: Integer.  The stoichiometry for the Modification.  Not sure
        if this will ever be greater than 1 but just in case.

        """
        if isinstance(modification, Modification):
            if modification in self._modifications:
                warn('Updating %s stoichiometry to %s from %s in Catalyst %s'%(modification.id,
                                                                             repr(stoichiometry),
                                                                             repr(self._modifications[modification]),
                                                                             self.id))
            else:
                self._modifications[modification] = stoichiometry
            self._update_id(self._modifications,prefix=self.id)
            modification._targets.add(self)

    def _update_id(self, stoichiometry_dict, prefix='Catalyst'):
        """Ids are derived from their basic composition following this structure:
        modification__v_k.name where v is the stoichiometric coefficient and k an Object and
        each modification is separated by __.

        """
        self.id = prefix
        _tmp_list = stoichiometry_dict.items()
        _tmp_list.sort()
        for k, v in _tmp_list:
            self.id += '__%s_%s'%(repr(v), k.id)
        _update_logic(self, stoichiometry_dict)




    
    def _remove_from_complex(self, complex):
        """
        """
        [x.remove_catalyst(self) for x in self.reactions]
        self._complex = None
        map(complex._modifications.remove, self.modifications)
        complex._catalysts.remove(self)


    def delete(self):
        model = self.model
        self._remove_from_complex(self.complex)
        [x.remove_target(self) for x in self.modifications]
        if model is not None:
            model.catalysts.remove(self)
            self._model = None
            for modification in self.modifications:
                if len(modification.targets) == 0:
                    model.modifications.remove(modification)
                    modification._model = None

        self._modifications = {}

    def add_to_model(self, model):
        if self.model is not model:
            if self.model is not None:
                raise(Exception('(%s) already associated with model (%s) cannot add to %s.'%(self.id, self.model.id,
                                                                                                  model.id)  +\
                                'Have not implemented code to change association to a different model'))

            model.catalysts.append(self)
            self._model = model
            self.complex.add_to_model(model)

            [x.add_to_model(model) for x in self.modifications]



        
class Subunit(Gene):
    """Subunits are Genes that have been transformed to a functional component state
    where they may be combined with other Subunits to form a Complex.

    A subunit represents an expressed sequence that can be added to a Complex.  It is either
    a polypeptide, aka protein, or a functional RNA (e.g., rRNA)
    
    """
    def __init__(self, id, name=None, pid=None, is_polypeptide=True):
        """

        pid: Integer.  Protein Identification number TODO: List source

        is_polypeptide: Boolean to distinguish between functional RNA units and polypeptides
        
        """
        if isinstance(id, Subunit):
            self._copy_parent_attributes(id)
        else:
            if isinstance(id, Gene):
                self._copy_parent_attributes(id)
                self._model = id._model
            else:
                Gene.__init__(self, id, name=name)
            self.pid = pid #Protein ID number
            self.is_polypeptide = is_polypeptide
        self._complex = set() #References the complexes that use this subunit
        #Whenever a subunit is created from another, reset the complexes
    
    def __getstate__(self):
        """Remove the references to container complexes when serializing to avoid
        problems associated with recursion.
        
        """
        state = Species.__getstate__(self)
        state['_complex'] = set()
        return state

    def guided_copy(self, the_model, the_complex=None):
        """Trying to make a faster copy procedure for cases where large
        numbers of speciess might be copied.  Such as when copying reactions.

        """
        the_copy = Gene.guided_copy(self, the_model)
        #Copy the more complex objects in a faster fashion

        the_copy._complex = set()
        return(the_copy)

    def remove_from_model(self):
        """

        """
        warn("%s.remove_from_model is not yet implemented"%(str(type(self))))

    def add_to_model(self, model):
        """
        model: cobra.core.me.ME_Model
        
        """

        if self.model is not model:
            if self.model is not None:
                raise(Exception('Gene (%s) already associated with model (%s) cannot add to %s.'%(self.id, self.model.id,
                                                                                                  model.id)  +\
                                'Have not implemented code to change association to a different model'))
            self._model = model
            model.subunits.append(self)
            model.genes = model.subunits
            [x.add_to_model(model) for x in self._complex]

    @property
    def complexes(self):
        return(self._complex)
#Module Functions
def _update_logic(the_object, stoichiometry_dict):
    """Construct a logical and string from a stoichiometry dictionary and
    attach to an object

    the_object: an object capable of having an attribute attached

    
    """
    tmp_list = stoichiometry_dict.items()
    tmp_list.sort()
    tmp_list = ['%s(%i)'%(k,v) for k, v in tmp_list]
    the_object.logic = reduce(lambda x, y: '%s AND %s'%(x, y), tmp_list)
