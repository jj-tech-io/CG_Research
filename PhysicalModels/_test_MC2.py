import unittest
import MC2


mc = MC2.Abs_Scat()
class _test_MC2(unittest.TestCase):
    def __init__(self) -> None:
        print(mc.get_alpha_base(400))
        print(mc.get_alpha_em(400))
        print(mc.get_alpha_ph(400))
    

if __name__ == "__main__":
    _test_MC2()