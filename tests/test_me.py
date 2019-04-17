from unittest import TestCase
import chimera

class MyTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testMe(self):
        self.assertEqual(42, 42)

    def testMeToo(self):
        self.assertEqual("Hello World", chimera.hello())

