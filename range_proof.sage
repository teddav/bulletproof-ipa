import random

load("ipa.sage")

p = 929
Fp = GF(p)
E = EllipticCurve(Fp, [5, 15])
Fr = GF(E.gens()[0].order())
R.<X> = Fr[]


def inner_product(a, b):
    assert len(a) == len(b)
    return sum([_a * _b for (_a, _b) in zip(a, b)])


def fold(vec, val):
    res = []
    for i in range(0, len(vec), 2):
        res.append(vec[i] * val + vec[i+1] * (1 / val))
    return res


n = 8  # number of bits
print(f"We will be proving that v is between 0 and {pow(2, n)}\n")

G = E.random_point()
H = E.random_point()

Gs = [E.random_point() for _ in range(n)]
Hs = [E.random_point() for _ in range(n)]

v = Fr(random.randint(0, pow(2, n)))
print("v =", v)

v_bin = bin(v)[2:].zfill(n)[::-1][:n]
print("v_bin = ", v_bin)

print("\nWe can commit to v from the start")
blinding_gamma = Fr.random_element()
V = v * G + blinding_gamma * H
print(f"v commitment (V): {V}\n")

aL = vector([Fr(int(bit)) for bit in v_bin])
assert v == sum([aL[i] * 2 ^ i for i in range(n)])

# (2^0, 2^1, 2^2, ..., 2^(n-1))
vec_2n = vector([Fr(2 ^ i) for i in range(n)])
assert v == inner_product(aL, vec_2n)

# (1, 1, 1, ..., 1)
vec_1n = vector([Fr(1)] * n)

aR = aL - vec_1n
assert inner_product(aL, aR) == 0

print("aL = ", aL)
print("aR = ", aR)

blinding_alpha = Fr.random_element()

A = inner_product(aL, Gs) + inner_product(aR, Hs) + blinding_alpha * H
print("A = ", A)

# linear terms for left and right polys
sL = vector([Fr.random_element() for i in range(n)])
sR = vector([Fr.random_element() for i in range(n)])

blinding_beta = Fr.random_element()

S = inner_product(sL, Gs) + inner_product(sR, Hs) + blinding_beta * H
print("S = ", S)

print("\nProver sends A, S, V to Verifier")
print("Verifier sends random challenges y and z\n")
y = Fr.random_element()
vec_y_n = vector([y ^ i for i in range(n)])

z = Fr.random_element()
vec_z_1n = vector([z] * n)

main_inner_product = inner_product(aL - vec_z_1n, aR.pairwise_product(vec_y_n) + vec_y_n * z + z ^ 2 * vec_2n)

delta_y_z = (z - z ^ 2) * inner_product(vec_1n, vec_y_n) - z ^ 3 * inner_product(vec_1n, vec_2n)

assert main_inner_product == z ^ 2 * v + delta_y_z
print("Combined inner product = z ^ 2 * v + delta_y_z.\nWe can continue...\n")

lX = aL - vec_z_1n + sL * X
print("lX = ", lX)
rX = vec_y_n.pairwise_product(aR + vec_z_1n) + z ^ 2 * vec_2n + vec_y_n.pairwise_product(sR * X)
print("rX = ", rX)

tX = inner_product(lX, rX)
print(f"tX = {tX}\n")

[t0, t1, t2] = tX.coefficients(sparse=False)

print("Notice that t0 is the inner product we're trying to prove\n")
assert t0 == main_inner_product

blinding_tau1 = Fr(123)
blinding_tau2 = Fr(456)
T1 = t1 * G + blinding_tau1 * H
T2 = t2 * G + blinding_tau2 * H
print("T1 = ", T1)
print("T2 = ", T2)

print("\nVerifier sends challenge x\n")
x = Fr.random_element()
print("x = ", x)

# evaluate left and right polys at u
lx = lX(x)
rx = rX(x)
tx = tX(x)
print("lx = ", lx)
print("rx = ", rx)
print("tx = ", tx)
assert tx == lx * rx

print("\nProver sends proof_blindings_mu and proof_blindings_tau to Verifier")
proof_blindings_mu = blinding_alpha + blinding_beta * x
proof_blindings_tau = z ^ 2 * blinding_gamma + blinding_tau1 * x + blinding_tau2 * x ^ 2

print("\nVerifier computes a new basis vector")
vec_y_n_inv = vec_y_n.apply_map(lambda x: 1/x)
H_y_minus1 = [vec_y_n_inv[i] * Hs[i] for i in range(n)]

print("\nFinal verification:")

check1_lhs = tx * G + proof_blindings_tau * H
check1_rhs = V * z ^ 2 + delta_y_z * G + T1 * x + T2 * x ^ 2
assert check1_lhs == check1_rhs
print("Check 1 ✅")

P = - proof_blindings_mu * H + A + S*x + inner_product(-vec_z_1n, Gs) + inner_product(z * vec_y_n + z ^ 2 * vec_2n, H_y_minus1)
assert P == inner_product(lx, Gs) + inner_product(rx, H_y_minus1) 
print("Check 2 ✅")

print("Verification successful")
