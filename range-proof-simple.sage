import random

load("ipa.sage")

p = 929
Fp = GF(p)
E = EllipticCurve(Fp, [5, 15])
Fr = GF(E.gens()[0].order())


n = 8  # number of bits
print(f"We will be proving that v is between 0 and {pow(2, n)}\n")

# (2^0, 2^1, 2^2, ..., 2^(n-1))
vec_2n = vector([Fr(2 ^ i) for i in range(n)])
# (1, 1, 1, ..., 1)
vec_1n = vector([Fr(1)] * n)

v = Fr(random.randint(0, pow(2, n)))
print("v =", v)

v_bin = bin(v)[2:].zfill(n)[::-1][:n]
print("v_bin = ", v_bin)

aL = vector([Fr(int(bit)) for bit in v_bin])
assert v == sum([aL[i] * 2 ^ i for i in range(n)])

assert v == inner_product(aL, vec_2n)

# Define aR
aR = aL - vec_1n
assert inner_product(aL, aR) == 0

print("aL = ", aL)
print("aR = ", aR)

# Define generators
G = E.random_point()
H = E.random_point()

Gs = [E.random_point() for _ in range(n)]
Hs = [E.random_point() for _ in range(n)]

# Commit to v
print("\nWe can commit to v from the start")
blinding_gamma = Fr.random_element()
V = v * G + blinding_gamma * H
print(f"v commitment (V): {V}\n")


blinding_alpha = Fr.random_element()
A = inner_product(aL, Gs) + inner_product(aR, Hs) + blinding_alpha * H
print("A = ", A)
print("\nProver sends A, V to Verifier")

print("Verifier sends random challenges y and z\n")
y = Fr.random_element()
vec_y_n = vector([y ^ i for i in range(n)])

z = Fr.random_element()
vec_z_1n = vector([z] * n)

l = aL - vec_z_1n
r = aR.pairwise_product(vec_y_n) + vec_y_n * z + z ^ 2 * vec_2n
main_inner_product = inner_product(l, r)

delta_y_z = (z - z ^ 2) * inner_product(vec_1n, vec_y_n) - z ^ 3 * inner_product(vec_1n, vec_2n)
t = z ^ 2 * v + delta_y_z

assert main_inner_product == t
print("Combined inner product = z ^ 2 * v + delta_y_z.\nWe can continue...\n")

vec_y_n_inv = vec_y_n.apply_map(lambda x: 1/x)
H_y_minus1 = [vec_y_n_inv[i] * Hs[i] for i in range(n)]

P = A - inner_product(vec_z_1n, Gs) + inner_product(z * vec_y_n + z ^ 2 * vec_2n, H_y_minus1) - blinding_alpha * H
assert P == inner_product(l, Gs) + inner_product(r, H_y_minus1)

print("\nFinally we run the IPA with: P + t * Q")
Q = E.random_point()

ipa_proof = ipa(l, r, Gs, H_y_minus1, t, Q, Fr)
P_full = P + t * Q
verify(Gs, H_y_minus1, P_full, ipa_proof[0], ipa_proof[1], ipa_proof[2], ipa_proof[3], ipa_proof[4], Q, n, Fr)
print("IPA proof ✅")
