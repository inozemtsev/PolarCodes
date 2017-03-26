import numpy as np
from scipy.optimize import brentq
from decoder import Decoder

__all__ = ['PolarBase', 'PolarCodeAWGN', 'PolarCodeBEC', 'PolarCodeBSC']


class CodeGenEmpirical:
    # simple empirical method

    def checknode_mean(self, t):
        return 2 * t - t**2

    def varnode_mean(self, t):
        return t**2


class CodeGenGaussian:
    # gaussian approximation

    def theta(self, x):
        return np.where(
            x >= 10,
            np.sqrt(np.pi / x) * np.exp(-x/4) * (1 - 10 / 7 / x),
            np.exp(-0.4527 * x**0.86 + 0.0218)
        )

    def theta1(self, x):
        # inverse function
        def _single_val(y):
            return brentq(lambda t: self.theta(t) - y, 1e-10, 1e10)
        return np.vectorize(_single_val)(x)

    def checknode_mean(self, t):
        return self.theta1(1 - (1 - self.theta(t))**2)

    def varnode_mean(self, t):
        return 2 * t


class PolarBase:
    def __init__(self, bhatt_param, N, K, kernel=np.array([[1, 0], [1, 1]]), mode=CodeGenGaussian()):
        if (np.log(N) / np.log(2)) % 1 != 0:
            raise ValueError('N is not a power of 2')

        self.N = N
        self.K = K
        self.kernel = kernel
        self.mode = mode

        n = int(np.log(N) / np.log(2))
        values = np.array([bhatt_param])

        for level in range(n):
            values = values[:, np.newaxis]
            values = np.hstack((mode.checknode_mean(values), mode.varnode_mean(values))).flatten()

        self.zvalues = values
        self.indices = np.argsort(self.zvalues)
        self.indices = np.array(list(map(lambda x: int(bin(x)[2:].zfill(n)[::-1], 2), self.indices)))
        self.frozen_indices = self.indices[K:]
        self.indices = self.indices[:K]

        for i in range(n-1):
            self.kernel = np.kron(self.kernel, kernel)

    def encode(self, message):
        if message.shape[0] != self.K:
            raise ValueError('Message length should be {}'.format(self.K))

        message_extended = np.zeros(self.N)
        message_extended[self.indices] = message
        encoded = (message_extended @ self.kernel) % 2
        return encoded
    
    def decode(self, message, L):
        u_frozen = np.zeros_like(message)
        u_frozen[self.indices] = -1
        codeword = Decoder(message, u_frozen, self.frozen_indices, L)
        return ((codeword @ self.kernel) % 2)[self.indices]


class PolarCodeAWGN(PolarBase):
    def __init__(self, Ec, N0, N, K, **kwargs):
        bhatt_param = np.exp(-Ec / N0)
        super(PolarCodeAWGN, self).__init__(bhatt_param, N, K, **kwargs)


class PolarCodeBEC(PolarBase):
    def __init__(self, epsilon, N, K, **kwargs):
        super(PolarCodeBEC, self).__init__(epsilon, N, K, **kwargs)
