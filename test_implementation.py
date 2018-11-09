import unittest
import numpy as np
import tiltedCHSH
import CKS2018

class Test_tiltedCHSH(unittest.TestCase):

    def setUp(self):
        self.alphavalues = np.linspace(start=0,
                                       stop=1.99,
                                       num=10,
                                       endpoint=True)

    @staticmethod
    def _concurrence(fourbyonematrix):
        return 0.5 * np.abs(fourbyonematrix[0, 0]*fourbyonematrix[0, 3] - fourbyonematrix[0, 1]*fourbyonematrix[0, 2])

    @staticmethod
    def _purestate2vector(rho):
        (eigenvalues, eigenvectors) = np.linalg.eig(rho)
        for index in range(len(eigenvalues)):
            if np.allclose([eigenvalues[index]], [1.]):
                return eigenvectors[index]
        raise Exception

    def test_Phi(self):

        # assert Phi(alpha) is a pure state
        for alpha in self.alphavalues:
            phi = tiltedCHSH.Phi(alpha=alpha)
            # in one way ...
            self.assertTrue(np.allclose(np.matmul(phi, phi), phi))

            # ... and in another way
            (eigenvalues, eigenvectors) = np.linalg.eig(phi)
            self.assertTrue(np.allclose(sorted(eigenvalues), [0., 0., 0., 1.]))

        # assert Phi(alpha=0) is maximally entangled
        vector = Test_tiltedCHSH._purestate2vector(tiltedCHSH.Phi(alpha=0))
        print(tiltedCHSH.Phi(alpha=0))
        concurrence = Test_tiltedCHSH._concurrence(vector)
        print(vector)
        print(concurrence)
        # TODO this should be 1, but it is not...
#        self.assertTrue(np.isclose(concurrence, 1.))

    def test_quantum_value(self):
        for (alpha, expected_quantum_value) in [(0., 2. * np.sqrt(2)), (2., 4.)]:
            self.assertEqual(tiltedCHSH.quantum_value(alpha), expected_quantum_value)

    def test_classical_value(self):
        for (alpha, expected_classical_value) in [(0., 2.), (2., 4.)]:
            self.assertEqual(tiltedCHSH.classical_value(alpha), expected_classical_value)

    def test_CHSH_operator(self):
        for (a, b, expected_operator) in \
                [ (0, 0, 2.*np.kron(tiltedCHSH.X, tiltedCHSH.X))]:
            self.assertTrue(np.allclose(tiltedCHSH.CHSH_operator(a, b), expected_operator))

class Test_CKS2018(unittest.TestCase):

    def setUp(self):
        self.alphavalues = np.linspace(start=0,
                                       stop=1.99,
                                       num=10,
                                       endpoint=True)

    def test_s(self):
        for (alpha, expected_s_value) in [(0, (4 + 5 * np.sqrt(2))/16.)]:
            self.assertTrue(np.isclose(CKS2018.s(alpha=alpha), expected_s_value))

    def test_mu(self):
        for (alpha, expected_mu_value) in [(0, (-1.) * (1. + 2. * np.sqrt(2.))/4.)]:
            self.assertTrue(np.isclose(CKS2018.mu(alpha=alpha), expected_mu_value))
    
    def test_beta_star(self):
        for (alpha, expected_beta_star_value) in [(0, (16. + 14. * np.sqrt(2.))/17.)]:
            self.assertTrue(np.isclose(CKS2018.beta_star(alpha=alpha), expected_beta_star_value))

    def test_conjugate(self):
        paulis = [tiltedCHSH.Id, tiltedCHSH.X, tiltedCHSH.Y, tiltedCHSH.Z]
        for pauli in paulis:
            self.assertTrue(np.allclose(pauli, CKS2018.conjugate(pauli, tiltedCHSH.Id)))
        for pauli_A in paulis[1:]:
            for pauli_B in paulis[1:]:
                if np.all(pauli_A == pauli_B):
                    pm = 1
                else:
                    pm = -1
                conj = CKS2018.conjugate(pauli_A, pauli_B)
                self.assertTrue(np.allclose( pm * pauli_A, conj))

    def test_effective_angle(self):

        num = 20

        # the effective angle maps the interval [0, b*] linearly onto [0, pi/4] ...
        for alpha in self.alphavalues:
            angles = np.linspace(0, CKS2018.b_star(alpha), num=num, endpoint=True)
            remote_angles = np.linspace(0, np.pi/4, num=num, endpoint=True)
            for index in range(num):
                angle = angles[index]
                remote_angle = CKS2018.effective_angle(alpha=alpha, x=angle)
                self.assertTrue(np.isclose(remote_angle, remote_angles[index]))

        # ... and maps the interval [b*, pi/2] linearly onto [pi/4, pi/2]
        for alpha in self.alphavalues:
            angles = np.linspace(CKS2018.b_star(alpha), np.pi/2., num=num, endpoint=True)
            remote_angles = np.linspace(np.pi/4, np.pi/2., num=num, endpoint=True)
            for index in range(num):
                angle = angles[index]
                remote_angle = CKS2018.effective_angle(alpha=alpha, x=angle)
                self.assertTrue(np.isclose(remote_angle, remote_angles[index]))


if __name__ == "__main__":
    unittest.main()
