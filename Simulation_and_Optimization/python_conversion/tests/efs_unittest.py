import unittest

from customize import efs_requirements

class TestEFS(unittest.TestCase):

    def test_nothing(self):
        self.assertTrue(True)


    def test_pipe_flow(self):
        print("")
        # create a material
        my_material = {'rho': 1, 'mu': 1}
        print(efs_requirements.pipe_flow(my_material, 2, 2))
        self.assertTrue(True)

    def test_head_loss(self):
        print("")
        print(efs_requirements.head_loss(2, 2, 2, 2))
        self.assertTrue(True)

if __name__ == "__main__":
    unittest.main()


