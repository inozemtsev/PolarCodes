import numpy as np

__all__ = ['PolarBase', 'PolarCodeAWGN', 'PolarCodeBEC', 'PolarCodeBSC']

class PolarBase:
    def __init__(self, bhatt_param, N, K, kernel):
        if (np.log(N) / np.log(2)) % 1 != 0:
            raise ValueError('N is not a power of 2')
        
        self.N = N
        self.K = K
        self.kernel = kernel
        
        n = int(np.log(N) / np.log(2))
        values = np.array([bhatt_param])
        
        for level in range(n):
            values = values[:, np.newaxis]
            values = np.hstack((2 * values - values ** 2, values ** 2)).flatten()
            
        self.indices = np.argsort(values)[:K]
        self.indices = np.array(list(map(lambda x: int(bin(x).zfill(n)[:1:-1], 2), self.indices)))
        
        for i in range(n-1):
            self.kernel = np.kron(self.kernel, kernel)
            
    def encode(self, message):
        if message.shape[0] != self.K:
            raise ValueError('Message length should be {}'.format(self.K))
        
        message_extended = np.zeros(self.N)
        message_extended[self.indices] = message
        encoded = self.kernel * message_extended
        return encoded
    
class PolarCodeAWGN(PolarBase):
    def __init__(self, Ec, N0, N, K, kernel=np.array([[1, 1], [0, 1]])):
        bhatt_param = np.exp(-Ec / N0)
        super(PolarCodeAWGN, self).__init__(bhatt_param, N, K, kernel)
        
class PolarCodeBEC(PolarBase):
    def __init__(self, epsilon, N, K, kernel=np.array([[1, 1], [0, 1]])):
        super(PolarCodeAWGN, self).__init__(epsilon, N, K, kernel)
        
class PolarCodeBSC(PolarBase):
    def __init__(self, p, N, K, kernel=np.array([[1, 1], [0, 1]])):
        bhatt_param = (p * (1 - p)) ** 0.5
        super(PolarCodeAWGN, self).__init__(bhatt_param, N, K, kernel)