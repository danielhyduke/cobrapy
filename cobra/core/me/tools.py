#me.tools.py
import re
from warnings import warn
def split_logical_string(logical_string):
    """Splits AND-based logical strings.
    'xxafd(2) AND Fdas()' will be returned as a
    dict: {'xxafd': 2, 'Fdas': 1}.  The Default action
    is to set the value to 1 if there is nothing in the
    parantheses.

    """
    and_re = re.compile(' +AND +', re.I)
    the_atoms = and_re.split(logical_string)
    logical_dict = {}
    
    for the_atom in the_atoms:
        #Deal with annoying GPR syntax
        if the_atom.count('(') > the_atom.count(')'):
            the_atom = the_atom[1:].strip()
        if the_atom.count('(') < the_atom.count(')'):
            the_atom = the_atom[:-1].strip()
        the_atom.strip()
        if '(' in the_atom:
            the_key, the_value = the_atom.split('(')
            the_value = the_value.rstrip(')')
            if len(the_value) == 0:
                the_value = 1
            else:
                the_value = int(the_value)
        else:
            the_key = the_atom
            the_value = 1
        logical_dict[the_key] = the_value

    return logical_dict



def save_complexes(complex_dict, complex_file_handle, modification_file_handle,
                   homology_status, data_source='E_coli_ME_20121115'):
    """Saves information about Complees and Catalysts to text files in a format
    similar the JT ME format.
    
    """
    for the_complex in complex_dict.itervalues():
        try:
            the_reactions = reduce(lambda x, y: '%s; %s'%(x, y), the_complex._reaction)
        except:
            the_reactions = ''
        complex_file_handle.write('\n%s\t%s\t%s\t%s\t%s\t%s'%(the_complex.id,
                                                              the_complex.name,
                                                              the_complex.logic,
                                                              data_source,
                                                              homology_status,
                                                              the_reactions))
        
        for the_enzyme in the_complex._catalysts:
            #Note: Currently assume only one modification per enzyme
            #TODO: Deal with unmodified catalysts
            if len(the_enzyme._modifications) == 0:
                continue
            the_modification = the_enzyme._modifications.keys()[0]
            _unmapped_modifiers = [x.id for x in the_modification._modification_dict
                                   if hasattr(x, '_dummy_metabolite')]
            if len(_unmapped_modifiers) == 0:
                _unmapped_modifiers = ''
            else:
                _unmapped_modifiers = reduce(lambda x, y: '%s AND %s'%(x, y), _unmapped_modifiers)
            modification_file_handle.write('\n%s\t%s\t%s\t%s\t%s\t%s'%(the_enzyme.id,
                                                              the_complex.id,
                                                              the_modification.logic,
                                                              data_source,
                                                              the_modification.id,
                                                              _unmapped_modifiers))
