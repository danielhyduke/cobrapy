from __future__ import absolute_import

import json
from warnings import warn

from .. import Model, Metabolite, Reaction, Formula
from ..external.six import iteritems, string_types

# Detect numpy types to replace them.
try:
    from numpy import float_, bool_
except ImportError:
    class float_:
        pass

    class bool_:
        pass

_DEFAULT_REACTION_ATTRIBUTES = {
    'id', 'name', 'subsystem', 'lower_bound', 'upper_bound',
    'objective_coefficient', 'notes', 'gene_reaction_rule', 'variable_kind'}

_DEFAULT_METABOLITE_ATTRIBUTES = {
    'id', 'annotation', 'charge', 'compartment', 'formula', 'name', 'notes',
    '_bound', '_constraint_sense'}


def _fix_type(value):
    """convert possible types to str, float, and bool"""
    # Because numpy floats can not be pickled to json
    if isinstance(value, string_types):
        return str(value)
    if isinstance(value, float_):
        return float(value)
    if isinstance(value, bool_):
        return bool(value)
    return value


def _from_dict(obj):
    """build a model from a dict"""
    if 'reactions' not in obj:
        raise Exception('JSON object has no reactions attribute. Cannot load.')
    model = Model()
    # add metabolites
    new_metabolites = []
    for metabolite in obj['metabolites']:
        new_metabolite = Metabolite()
        for k, v in iteritems(metabolite):
            setattr(new_metabolite, k, _fix_type(v))
        new_metabolites.append(new_metabolite)
    model.add_metabolites(new_metabolites)
    # add reactions
    new_reactions = []
    for reaction in obj['reactions']:
        new_reaction = Reaction()
        for k, v in iteritems(reaction):
            if k == 'reversibility' or k == "reaction":
                continue
            elif k == 'metabolites':
                new_reaction.add_metabolites(
                    {model.metabolites.get_by_id(str(met)): coeff
                     for met, coeff in iteritems(v)})
            else:
                setattr(new_reaction, k, _fix_type(v))
        new_reactions.append(new_reaction)
    model.add_reactions(new_reactions)
    for k, v in iteritems(obj):
        if k in ['id', 'description', 'notes']:
            setattr(model, k, v)
    return model


def _to_dict(model):
    """convert the model to a dict"""
    reaction_attributes = _DEFAULT_REACTION_ATTRIBUTES
    metabolite_attributes = _DEFAULT_METABOLITE_ATTRIBUTES
    new_reactions = []
    new_metabolites = []
    for reaction in model.reactions:
        new_reaction = {key: _fix_type(getattr(reaction, key))
                        for key in reaction_attributes}
        # set metabolites
        mets = {str(met): coeff for met, coeff
                in iteritems(reaction._metabolites)}
        new_reaction['metabolites'] = mets
        new_reactions.append(new_reaction)
    for metabolite in model.metabolites:
        new_metabolite = {key: str(getattr(metabolite, key))
                          for key in metabolite_attributes}
        new_metabolites.append(new_metabolite)
    obj = {'reactions': new_reactions,
           'metabolites': new_metabolites,
           'id': model.id,
           'description': model.description,
           'notes': model.notes}
    return obj


def to_json(model):
    """Save the cobra model as a json string"""
    return json.dumps(_to_dict(model))


def from_json(jsons):
    """Load cobra model from a json string"""
    return _from_dict(json.loads(jsons))


def load_json_model(file_name):
    """Load a cobra model stored as a json file

    file_name : str or file-like object

    """
    # open the file
    should_close = False
    if isinstance(file_name, string_types):
        file_name = open(file_name, 'r')
        should_close = True

    model = _from_dict(json.load(file_name))

    if should_close:
        file_name.close()

    return model


def save_json_model(model, file_name):
    """Save the cobra model as a json file.

    model : :class:`~cobra.core.Model.Model` object

    file_name : str or file-like object

    """
    # open the file
    should_close = False
    if isinstance(file_name, string_types):
        file_name = open(file_name, 'w')
        should_close = True

    json.dump(_to_dict(model), file_name)

    if should_close:
        file_name.close()
