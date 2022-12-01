import unittest
import MC


mc = MC.Abs_Scat()
class testMC:
    def __init__(self) -> None:
        print(mc.get_alpha_base(400))
        print(mc.get_alpha_em(400))
        print(mc.get_alpha_ph(400))

if __name__ == "__main__":
    testMC()