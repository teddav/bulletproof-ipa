import math

p = 929
Fp = GF(p)
E = EllipticCurve(Fp, [5, 15])
G = E.gens()[0]
Fr = GF(G.order())
assert is_prime(G.order())

print("We define Q as a random point on the curve")
Q = E.random_point()


def inner_product(a, b):
    assert len(a) == len(b)
    return sum([_a * _b for (_a, _b) in zip(a, b)])


def fold(vec, val):
    assert len(vec) % 2 == 0
    half = len(vec) // 2
    left = vec[:half]
    right = vec[half:]
    return [left[i] * val + right[i] * (1 / val) for i in range(half)]


def get_LR(a, b, Gs, Hs):
    assert len(a) == len(b)
    assert len(a) % 2 == 0
    half = len(a) // 2
    aL = a[:half]
    aR = a[half:]
    bL = b[:half]
    bR = b[half:]
    GL = Gs[:half]
    GR = Gs[half:]
    HL = Hs[:half]
    HR = Hs[half:]
    L = inner_product(aL, GR) + inner_product(bR, HL) + \
        inner_product(aL, bR) * Q
    R = inner_product(aR, GL) + inner_product(bL, HR) + \
        inner_product(aR, bL) * Q
    return L, R


# Length of the vectors (this needs to be even length)
N = 8

a = [Fr.random_element() for _ in range(N)]
b = [Fr.random_element() for _ in range(N)]

# random elements
Gs = [E.random_point() for _ in range(N)]
Hs = [E.random_point() for _ in range(N)]

c = inner_product(a, b)
print("We want to prove that c = a * b")
print("c = ", c)

P = inner_product(a, Gs) + inner_product(b, Hs) + c * Q
print("P = ", P)

print("\nRound 1:")

L1, R1 = get_LR(a, b, Gs, Hs)

x1 = Fr.random_element()
print(f"\nRandom challenge x1: {x1}\n")

a_prime1 = fold(a, x1)
b_prime1 = fold(b, 1 / x1)
G_prime1 = fold(Gs, 1 / x1)
H_prime1 = fold(Hs, x1)

P_prime1 = inner_product(a_prime1, G_prime1) + inner_product(b_prime1,
                                                             H_prime1) + inner_product(a_prime1, b_prime1) * Q

assert P_prime1 == P + x1 ^ 2 * L1 + x1 ^ (-2) * R1
print("Round 1 successful")
print("Verifier receives L1, R1\n")

print("\nRound 2:")
L2, R2 = get_LR(a_prime1, b_prime1, G_prime1, H_prime1)

x2 = Fr.random_element()
print(f"\nRandom challenge x2: {x2}\n")

a_prime2 = fold(a_prime1, x2)
b_prime2 = fold(b_prime1, 1 / x2)
G_prime2 = fold(G_prime1, 1 / x2)
H_prime2 = fold(H_prime1, x2)

P_prime2 = inner_product(a_prime2, G_prime2) + inner_product(b_prime2,
                                                             H_prime2) + inner_product(a_prime2, b_prime2) * Q

assert P_prime2 == P_prime1 + x2 ^ 2 * L2 + x2 ^ (-2) * R2
print("Round 2 successful")
print("Verifier receives L2, R2\n")

print("\nRound 3:")
L3, R3 = get_LR(a_prime2, b_prime2, G_prime2, H_prime2)

x3 = Fr.random_element()
print(f"\nRandom challenge x3: {x3}\n")

a_prime3 = fold(a_prime2, x3)
b_prime3 = fold(b_prime2, 1 / x3)
G_prime3 = fold(G_prime2, 1 / x3)
H_prime3 = fold(H_prime2, x3)

P_prime3 = inner_product(a_prime3, G_prime3) + inner_product(b_prime3,
                                                             H_prime3) + inner_product(a_prime3, b_prime3) * Q

assert P_prime3 == P_prime2 + x3 ^ 2 * L3 + x3 ^ (-2) * R3
print("Round 3 successful")
print("Verifier receives a_prime3, b_prime3, L3, R3\n")

print("\nVerification:")
print("Verifier now has: a_prime3, b_prime3, (L1, R1), (L2, R2), (L3, R3)\n")

print("He recomputes P_prime3:")
print("First needs to recompute G_prime3 and H_prime3, from the initial Gs and Hs")
verifier_G_prime1 = fold(Gs, 1 / x1)
verifier_G_prime2 = fold(verifier_G_prime1, 1 / x2)
verifier_G_prime3 = fold(verifier_G_prime2, 1 / x3)
assert verifier_G_prime3 == G_prime3
verifier_H_prime1 = fold(Hs, x1)
verifier_H_prime2 = fold(verifier_H_prime1, x2)
verifier_H_prime3 = fold(verifier_H_prime2, x3)
assert verifier_H_prime3 == H_prime3
print("verifier_G_prime3 = ", verifier_G_prime3)
print("verifier_H_prime3 = ", verifier_H_prime3)

verifier_P_prime3 = inner_product(a_prime3, verifier_G_prime3) + inner_product(
    b_prime3, verifier_H_prime3) + inner_product(a_prime3, b_prime3) * Q
print("verifier_P_prime3 = ", verifier_P_prime3)

final_check = verifier_P_prime3 \
    - (x1 ^ 2 * L1 + x1 ^ (-2) * R1) \
    - (x2 ^ 2 * L2 + x2 ^ (-2) * R2) \
    - (x3 ^ 2 * L3 + x3 ^ (-2) * R3)
assert final_check == P

print("Verification successful!")

print("\nBut recomputing P_prime3 is expensive for the verifier")
print("Let's improve it")
xs = [x1, x2, x3]

s_G = [Fr(1) for _ in range(N)]
s_H = [Fr(1) for _ in range(N)]

for i in range(N):
    i_bits = bin(i)[2:].zfill(int(math.log(N, 2)))
    for j in range(int(math.log(N, 2))):
        bit_j = int(i_bits[j])
        if bit_j == 1:
            s_G[i] *= xs[j]
            s_H[i] *= 1 / xs[j]
        else:
            s_G[i] *= 1 / xs[j]
            s_H[i] *= xs[j]

print("s_G = ", s_G)
print("s_H = ", s_H)

G_final = inner_product(s_G, Gs)
H_final = inner_product(s_H, Hs)
assert G_final == G_prime3[0]
assert H_final == H_prime3[0]

print("Verification (with optimization) successful!")
