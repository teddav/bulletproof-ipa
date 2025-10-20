p = 929
Fp = GF(p)


def fold(vec, val):
    assert len(vec) % 2 == 0
    half = len(vec) // 2
    left = vector(vec[:half])
    right = vector(vec[half:])
    return left * val + right * (1 / val)


def get_LR(a, b):
    assert len(a) == len(b)
    assert len(a) % 2 == 0
    half = len(a) // 2
    aL = vector(a[:half])
    aR = vector(a[half:])
    bL = vector(b[:half])
    bR = vector(b[half:])
    L = aL * bR
    R = aR * bL
    return L, R

# We have 2 vectors `a_real` and `b_real` we want to prove the inner product of:
a_real = vector(Fp, [458, 593, 443, 699, 142, 338, 714, 91])
b_real = vector(Fp, [721, 164, 121, 550, 545, 789, 261, 510])

c = a_real * b_real
print("We want to prove that c = a * b")
print("c = ", c)

# we could easily find 2 fake vectors whose inner product equals `c` too:
a = vector(Fp, [0, 0, 0, 0, 0, 0, 0, 1])
b = vector(Fp, [0, 0, 0, 0, 0, 0, 0, c])

assert a * b == c

print("\nRound 1:")
c1 = a * b
print("inner product round 1 = ", c1)
L1, R1 = get_LR(a, b)

x1 = Fp.random_element()
print(f"\nRandom challenge x1: {x1}\n")

a_prime1 = fold(a, x1)
b_prime1 = fold(b, 1 / x1)
print("a_prime1 = ", a_prime1)
print("b_prime1 = ", b_prime1)

c_prime1 = c1 + x1 ^ 2 * L1 + x1 ^ (-2) * R1
assert c_prime1 == a_prime1 * b_prime1
print("Round 1 successful")
print("Verifier receives L1, R1\n")

print("\nRound 2:")
c2 = a_prime1 * b_prime1
print("inner product round 2 = ", c2)
L2, R2 = get_LR(a_prime1, b_prime1)

x2 = Fp.random_element()
print(f"\nRandom challenge x2: {x2}\n")

a_prime2 = fold(a_prime1, x2)
b_prime2 = fold(b_prime1, 1 / x2)
print("a_prime2 = ", a_prime2)
print("b_prime2 = ", b_prime2)

c_prime2 = c2 + x2 ^ 2 * L2 + x2 ^ (-2) * R2
assert c_prime2 == a_prime2 * b_prime2
print("Round 2 successful")
print("Verifier receives L2, R2\n")

print("\nRound 3:")
c3 = a_prime2 * b_prime2
print("inner product round 3 = ", c3)
L3, R3 = get_LR(a_prime2, b_prime2)

x3 = Fp.random_element()
print(f"\nRandom challenge x3: {x3}\n")

a_prime3 = fold(a_prime2, x3)
b_prime3 = fold(b_prime2, 1 / x3)
print("a_prime3 = ", a_prime3)
print("b_prime3 = ", b_prime3)

c_prime3 = c3 + x3 ^ 2 * L3 + x3 ^ (-2) * R3
assert c_prime3 == a_prime3 * b_prime3
print("Round 3 successful", c_prime3)
print("Verifier receives L3, R3\n")

print("a and b now have length 1")
print("We can send a_prime3 and b_prime3 to the verifier")

print("\nVerifier now has: a_prime3, b_prime3, L3, R3, L2, R2, L1, R1\n")

assert c == a_prime3 * b_prime3 \
    - (x1 ^ 2 * L1 + x1 ^ (-2) * R1) \
    - (x2 ^ 2 * L2 + x2 ^ (-2) * R2) \
    - (x3 ^ 2 * L3 + x3 ^ (-2) * R3)

print("Verification successful!")
