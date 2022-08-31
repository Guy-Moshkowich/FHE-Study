import unittest
import sys
import os
sys.path.append(os.path.abspath("RLWE/"))
sys.path.append(os.path.abspath("Utils/"))
sys.path.append(os.path.abspath("BGV/"))

testmodules = [
    'ring_element_test',
    'bgv_test',
    'utils_test',
    ]


suite = unittest.TestSuite()
for t in testmodules:
    print(t)
    suite.addTest(unittest.defaultTestLoader.loadTestsFromName(t))
unittest.TextTestRunner().run(suite)