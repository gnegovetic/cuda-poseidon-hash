MODULUS = 8444461749428370424248824938781546531375899335154063827935233455917409239041

class Fq:
    def __init__(self, value):
        self.value = value % MODULUS
    def __add__(self, other): return Fq(self.value + int(other))
    def __sub__(self, other): return Fq(self.value - int(other))
    def __mul__(self, other): return Fq(self.value * int(other))
    def __truediv__(self, other): return self * Fq(pow(int(other), -1, MODULUS))
    def __pow__(self, power, modulo=None): return Fq(pow(self.value, int(power), MODULUS))
    def __eq__(self, other): return self.value == int(other)
    def __int__(self): return self.value
    def __repr__(self): return f"Fq({self.value})"

T = 3       # state width
RATE = 2    # elements per absorb/squeeze block
CAPACITY = 1
FULL_ROUNDS = 8
PARTIAL_ROUNDS = 57
ALPHA = 5
# Not real Poseidon parameters! Just for demonstration.
ARK = [[Fq(i + j + 1) for j in range(T)] for i in range(FULL_ROUNDS + PARTIAL_ROUNDS)]
MDS = [
    [Fq(2), Fq(1), Fq(1)],
    [Fq(1), Fq(2), Fq(1)],
    [Fq(1), Fq(1), Fq(2)],
]

def print_fqs(name: str, numbers: list[Fq]):
    ''' Splits and prints each state value into 12 little-endian 32-bit words'''
    print(f"{name}:")
    for val in numbers:
        words = [(int(val) >> (32 * i)) & 0xFFFFFFFF for i in range(12)]
        print(" ".join(f"{i}:0x{w:08x}" for i, w in enumerate(words)))

class PoseidonPermutation:
    def __init__(self):
        self.state = [Fq(0)] * T

    def permute(self):
        round_idx = 0
        # Full rounds (first half)
        for _ in range(FULL_ROUNDS // 2):
            self.apply_ark(round_idx)
            self.apply_sbox(full=True)
            self.apply_mds()
            round_idx += 1
        # Partial rounds
        for _ in range(PARTIAL_ROUNDS):
            self.apply_ark(round_idx)
            self.apply_sbox(full=False)
            self.apply_mds()
            round_idx += 1
        # Full rounds (second half)
        for _ in range(FULL_ROUNDS // 2):
            self.apply_ark(round_idx)
            self.apply_sbox(full=True)
            self.apply_mds()
            round_idx += 1

    def apply_ark(self, round_idx):
        for i in range(T):
            self.state[i] += ARK[round_idx][i]

    def apply_sbox(self, full=True):
        if full:
            for i in range(T):
                self.state[i] = self.state[i] ** ALPHA
        else:
            self.state[0] = self.state[0] ** ALPHA

    def apply_mds(self):
        new_state = []
        for i in range(T):
            acc = Fq(0)
            for j in range(T):
                acc += MDS[i][j] * self.state[j]
            new_state.append(acc)
        self.state = new_state

class PoseidonSponge:
    def __init__(self):
        self.perm = PoseidonPermutation()
        self.rate = RATE
        self.capacity = CAPACITY
        self.state = [Fq(0)] * T
        self.pos = 0
        self.mode = "absorbing"

    def absorb(self, inputs):
        for x in inputs:
            if self.pos == self.rate:
                self.permute()
                self.pos = 0
            self.state[self.pos] += Fq(x)
            self.pos += 1

    def permute(self):
        self.perm.state = self.state.copy()
        self.perm.permute()
        self.state = self.perm.state.copy()

    def squeeze(self, num_elements):
        outputs = []
        if self.mode != "squeezing":
            self.permute()
            self.pos = 0
            self.mode = "squeezing"
        while len(outputs) < num_elements:
            if self.pos == self.rate:
                self.permute()
                self.pos = 0
            outputs.append(int(self.state[self.pos]))
            self.pos += 1
        return outputs

    def hash(self, inputs, num_outputs=1):
        """
        Convenience function: absorbs input list and squeezes num_outputs field elements.
        Typical hash usage: num_outputs=1.
        """
        self.state = [Fq(0)] * T
        self.pos = 0
        self.mode = "absorbing"
        self.absorb(inputs)
        # print_fqs("After absorbing", self.state)
        outputs = self.squeeze(num_outputs)
        print_fqs("After squeeze", self.state)
        return outputs if num_outputs > 1 else outputs[0]

    def hash_bytes(self, data: bytes, num_outputs=1):
        """
        Hash bytes: splits the input into chunks < modulus and feeds as field elements.
        Default: each field element is 47 bytes (376 bits, < 377-bit modulus).
        """
        CHUNK_SIZE = 47
        elements = []
        for i in range(0, len(data), CHUNK_SIZE):
            chunk = data[i:i+CHUNK_SIZE]
            elem = int.from_bytes(chunk, byteorder='little')
            elements.append(elem)
        return self.hash(elements, num_outputs)


if __name__ == "__main__":

    # Hash 256 bytes (0..255)
    # sponge = PoseidonSponge()
    # digest = sponge.hash_bytes(bytes(range(256)), num_outputs=2)

    # print("Poseidon hash_bytes(0..255) = ")
    # print(f"[{digest[0]:064x}],")
    # print(f"[{digest[1]:064x}]")

    #print_fqs("modulus", [MODULUS])

    # Test with a single field element
    sponge = PoseidonSponge()
    single_element_hash = sponge.hash([77], num_outputs=1)
    print_fqs("Hash", [single_element_hash])
