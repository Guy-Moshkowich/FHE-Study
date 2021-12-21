import unittest

testmodules = [
    'ring_test',
    'bgv_test',
    'utils_test'
    ]


suite = unittest.TestSuite()
for t in testmodules:
    print(t)
    suite.addTest(unittest.defaultTestLoader.loadTestsFromName(t))
unittest.TextTestRunner().run(suite)
