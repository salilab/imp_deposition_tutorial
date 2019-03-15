#!/usr/bin/env python

import unittest
import os
import sys
import subprocess
import ihm.reader
import ihm.cross_linkers

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..',
                                      '..', 'rnapolii'))

class Tests(unittest.TestCase):

    def test_mmcif(self):
        """Test generation of mmCIF file"""
        # Run modeling
        os.chdir(os.path.join(TOPDIR, 'modeling'))
#       p = subprocess.check_call(["python", 'modeling.py', "--mmcif",
#                                  "--dry-run"])
        self._check_mmcif_file("rnapolii.cif")

    def _check_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 1)
        self.assertEqual(s.citations[0].doi, '10.1074/mcp.M114.041673')
        self.assertEqual(len(s.software), 2)
        self.assertEqual(len(s.orphan_starting_models), 12)
        # Should be a single state, of one model
        self.assertEqual(len(s.state_groups), 1)
        self.assertEqual(len(s.state_groups[0]), 1)
        self.assertEqual(len(s.state_groups[0][0]), 1)
        model = s.state_groups[0][0][0][0]
        self.assertEqual(len(model._atoms), 0)
        self.assertEqual(len(model._spheres), 3949)
        # Should be 1 ensemble (cluster) of 100 models
        self.assertEqual([e.num_models for e in s.ensembles], [100])
        # Check localization densities
        self.assertEqual([len(e.densities) for e in s.ensembles], [1])
        # Check components
        self.assertEqual([len(e.sequence) for e in s.entities],
                         [1733, 1224, 318, 218, 215, 155, 171, 146, 122,
                          70, 120, 70])
        self.assertEqual([a.details for a in s.asym_units],
                         ['Rpb1.0', 'Rpb2.0', 'Rpb3.0', 'Rpb4.0', 'Rpb5.0',
                          'Rpb6.0', 'Rpb7.0', 'Rpb8.0', 'Rpb9.0', 'Rpb10.0',
                          'Rpb11.0', 'Rpb12.0'])
        # 3 restraints - 2 crosslinks, 1 EM3D map
        self.assertEqual(len(s.restraints), 3)
        xl1, xl2 = s.restraints[:2]
        self.assertEqual(xl1.linker.auth_name, ihm.cross_linkers.dss.auth_name)
        self.assertEqual(len(xl1.experimental_cross_links), 157)
        self.assertEqual(len(xl1.cross_links), 153)
        self.assertEqual(xl1.dataset.location.path, 'data/polii_xlinks.csv')
        self.assertEqual(xl2.linker.auth_name, ihm.cross_linkers.bs3.auth_name)
        em = s.restraints[2]
        self.assertEqual(em.dataset.location.path,
                         'data/emd_1883.map.mrc.gmm.50.txt')
        self.assertEqual(len(em.dataset.parents), 1)
        p = em.dataset.parents[0]
        self.assertEqual(p.location.db_name, 'EMDB')
        self.assertEqual(p.location.access_code, 'EMD-1883')

if __name__ == '__main__':
    unittest.main()
