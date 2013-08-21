from __future__ import with_statement
import sys
from warnings import warn  # TODO - catch known warnings
from unittest import TestCase, TestLoader, TextTestRunner
from tempfile import gettempdir
from os import unlink
from copy import deepcopy
from os.path import join
try:  #skipIf is not in python 2.6 / 2.5
    from unittest import skipIf
except:
    try:  # should we make unittest2 an absolute requirement and use skipIf
        from unittest2 import skipIf
    except:
        skipIf = None

# deal with absolute imports by adding the appropriate directory to the path
if __name__ == "__main__":
    sys.path.insert(0, "../../../../")
    from cobra.core.me.test import data_directory, create_test_model
    from cobra import Object, Model, Metabolite, Reaction, io, DictList
    from cobra.core.me import ME_Model, CatalyzedReaction, Catalyst, Subunit, Enzyme, Complex, Modification
    sys.path.pop(0)
    #assert 0
else:
    from . import data_directory, create_test_model
    from .... import Object, Model, Metabolite, Reaction, io, DictList
    from .. import ME_Model, CatalyzedReaction, Catalyst, Subunit, Enzyme, Complex, Modification


class TestME_Model(TestCase):
    def setUp(self):
        self.m_model = create_test_model()
        #Note: self.me_model should not have any catalysts because there isn't any
        #information about modifications in the m_model
        self.me_model = ME_Model(self.m_model)



    def testIdentity(self):
        """When converting an M-Model to an ME-Model, the ids should be the same but the objects should be different.
        
        """
        self.assertEqual(self.m_model, self.me_model)
        self.assertIsNot(self.m_model, self.me_model)
        m_reaction = self.m_model.reactions[0]
        me_reaction = self.me_model.reactions.get_by_id(m_reaction.id)
        self.assertEqual(m_reaction, me_reaction)
        self.assertIsNot(m_reaction, me_reaction)
        m_gene = self.m_model.genes[0]
        me_gene = self.me_model.genes.get_by_id(m_gene.id)
        self.assertEqual(m_gene, me_gene)
        self.assertIsNot(m_gene, me_gene)

    def test_copy(self):
        """modifying copy should not modify the original"""
        # test that deleting reactions in the copy does not change the
        # number of reactions in the original model
        model_copy = self.me_model.copy()
        if len(self.me_model.catalysts) == 0:
            self.assertEqual(len(model_copy.catalysts), len(model_copy.complexes))
        old_reaction_count = len(self.me_model.reactions)
        self.assertEqual(len(self.me_model.reactions), len(model_copy.reactions))
        self.assertEqual(len(self.me_model.metabolites),
            len(model_copy.metabolites))
        model_copy.remove_reactions(model_copy.reactions[0:5])
        self.assertEqual(old_reaction_count, len(self.me_model.reactions))
        self.assertNotEqual(len(self.me_model.reactions),
            len(model_copy.reactions))


    def test_deepcopy(self):
        """Verify that reference structures are maintained when deepcopying.
        
        """
        model_copy = deepcopy(self.me_model)
        for gene, gene_copy in zip(self.me_model.genes, model_copy.genes):
            self.assertEqual(gene.id, gene_copy.id)
            reactions = list(gene.get_reaction())
            reactions.sort()
            reactions_copy = list(gene_copy.get_reaction())
            reactions_copy.sort()
            self.assertEqual(reactions, reactions_copy)
        for reaction, reaction_copy in zip(self.me_model.reactions, model_copy.reactions):
            self.assertEqual(reaction.id, reaction_copy.id)
            metabolites = reaction._metabolites.keys()
            metabolites.sort()
            metabolites_copy = reaction_copy._metabolites.keys()
            metabolites_copy.sort()
            self.assertEqual(metabolites, metabolites_copy)
        
        #TODO: Add in tests for Complexes, Catalysts, and Modifications

## class TestCobraIO(CobraTestCase):
##     pass
##     try:
##         from libsbml import SBMLReader as __SBMLReader
##         __test_sbml = True
##     except:
##         __test_sbml = False
##     if __test_sbml:
##         def test_sbml_read(self):
##             ## with catch_warnings(record=True) as w:
##             model = io.read_sbml_model(test_sbml_file)
##             self.assertEqual(len(model.reactions), len(self.model.reactions))
##             # make sure that an error is raised when given a nonexistent file
##             self.assertRaises(IOError, io.read_sbml_model,
##                               "fake_file_which_does_not_exist")
##         def test_sbml_write(self):
##             test_output_filename = join(gettempdir(), 'test_sbml_write.xml')
##             io.write_sbml_model(self.model, test_output_filename)
##             #cleanup the test file
##             unlink(test_output_filename)
##     try:
##         from cobra.io import save_matlab_model
##         __test_matlab = True
##     except:
##         __test_matlab = False
##     if __test_matlab:
##         def test_mat_read_write(self):
##             test_output_filename = join(gettempdir(), "test_mat_write.mat")
##             io.save_matlab_model(self.model, test_output_filename)
##             reread = io.load_matlab_model(test_output_filename)
##             self.assertEqual(len(self.model.reactions), len(reread.reactions))
##             self.assertEqual(len(self.model.metabolites), len(reread.metabolites))
##             for i in range(len(self.model.reactions)):
##                 self.assertEqual(len(self.model.reactions[i]._metabolites), \
##                     len(reread.reactions[i]._metabolites))
##                 self.assertEqual(self.model.reactions[i].id, reread.reactions[i].id)
##             unlink(test_output_filename)

# make a test suite to run all of the tests
loader = TestLoader()
#suite = loader.loadTestsFromModule(sys.modules[__name__])

def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()
