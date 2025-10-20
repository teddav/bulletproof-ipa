load("ipa.sage")

p = 929
Fp = GF(p)
E = EllipticCurve(Fp, [5, 15])
G = E.gens()[0]
Fr = GF(G.order())
R.<x> = Fr[]


def inner_product(a, b):
    assert len(a) == len(b)
    return sum([_a * _b for (_a, _b) in zip(a, b)])


def fold(vec, val):
    res = []
    for i in range(0, len(vec), 2):
        res.append(vec[i] * val + vec[i+1] * (1 / val))
    return res


n = 8  # number of bits

H = G * Fr(123)
Q = G * Fr(456)
B = G * Fr(789)

Gs = [i * 5 * G for i in range(1, n + 1)]
Hs = [i * 5 * H for i in range(1, n + 1)]

v = 159
v_bin = bin(v)[2:].zfill(n)[::-1][:n]
print("v_bin = ", v_bin)
aL = vector([Fr(int(bit)) for bit in v_bin])
assert v == sum([aL[i] * 2 ^ i for i in range(n)])

# (2^0, 2^1, 2^2, ..., 2^(n-1))
vec_2 = vector([Fr(2 ^ i) for i in range(n)])
assert v == inner_product(aL, vec_2)

# (1, 1, 1, ..., 1)
vec_1 = vector([Fr(1) for _ in range(n)])

aR = aL - vec_1
assert inner_product(aL, aR) == 0

print("aL = ", aL)
print("aR = ", aR)

# blinding_alpha = Fr.random_element()
blinding_alpha = Fr(580)

A = inner_product(aL, Gs) + inner_product(aR, Hs) + blinding_alpha * B
print("A = ", A)

# linear terms for left and right polys
# sL = vector([Fr.random_element() for i in range(n)])
# sR = vector([Fr.random_element() for i in range(n)])
sL = vector([Fr(12), Fr(838), Fr(190), Fr(510),
            Fr(232), Fr(512), Fr(1), Fr(813)])
sR = vector([Fr(566), Fr(789), Fr(732), Fr(902),
            Fr(910), Fr(387), Fr(819), Fr(868)])
# blinding_beta = Fr.random_element()
blinding_beta = Fr(543)

S = inner_product(sL, Gs) + inner_product(sR, Hs) + blinding_beta * B
print("S = ", S)

# blinding_gamma = Fr.random_element()
blinding_gamma = Fr(218)
V = v * G + blinding_gamma * B

print("\nProver sends A, S, V to Verifier")
print("Verifier sends random challenges y and z\n")
# y = Fr.random_element()
# z = Fr.random_element()
y = Fr(800)
z = Fr(900)
y_n = vector([y ^ i for i in range(n)])
z_1 = vector([z] * n)

lhs = inner_product(aL - z_1, aR.pairwise_product(y_n) + y_n * z + z ^ 2 * vec_2)

delta_y_z = (z - z ^ 2) * inner_product(vec_1, y_n) - z ^ 3 * inner_product(vec_1, vec_2)
rhs = z ^ 2 * v + delta_y_z

print("lhs = ", lhs)
print("rhs = ", rhs)

assert lhs == rhs
print("Verification successful\n")

lx = aL - z_1 + sL * x
print("lx = ", lx)
rx = y_n.pairwise_product(aR + z_1) + z ^ 2 * vec_2 + y_n.pairwise_product(sR * x)
print("rx = ", rx)

tx = inner_product(lx, rx)
print(f"tx = {tx}\n")

[t0, t1, t2] = tx.coefficients()

print("Notice that t0 is the inner product we're trying to prove\n")
assert t0 == lhs

blinding_t1 = Fr(123)
blinding_t2 = Fr(456)
T1 = t1 * G + blinding_t1 * B
T2 = t2 * G + blinding_t2 * B
print("T1 = ", T1)
print("T2 = ", T2)

print("\nVerifier sends challenge u\n")
# u = Fr.random_element()
u = Fr(100)
print("u = ", u)

# evaluate left and right polys at u
lu = lx(u)
ru = rx(u)
tu = tx(u)
print("lu = ", lu)
print("ru = ", ru)
print("tu = ", tu)

assert tu == lu * ru

print("\nProver sends proof_blindings1 and proof_blindings2 to Verifier")
proof_blindings1 = blinding_alpha + blinding_beta * u
proof_blindings2 = z ^ 2 * blinding_gamma + blinding_t1 * u + blinding_t2 * u ^ 2

print("\nVerifier computes a new basis vector")
y_n_inv = y_n.apply_map(lambda x: 1/x)
H_y_minus1 = [y_n_inv[i] * Hs[i] for i in range(n)]

print("\nFinal verification:")
print("Step 1:")
final_lhs = A + S*u + inner_product(-z_1, Gs) + inner_product(z * y_n + z ^ 2 * vec_2, H_y_minus1)
final_rhs = inner_product(lu, Gs) + inner_product(ru, H_y_minus1) + proof_blindings1 * B

assert final_lhs == final_rhs
print("Step 1 ✅")

final2_lhs = tu * G + proof_blindings2 * B
final2_rhs = V * z ^ 2 + delta_y_z * G + T1 * u + T2 * u ^ 2

assert final2_lhs == final2_rhs
print("Step 2 ✅")

print("Verification successful")
