import math

def inner_product(a, b):
    assert len(a) == len(b)
    return sum([_a * _b for (_a, _b) in zip(a, b)])


def fold(vec, val):
    assert len(vec) % 2 == 0
    half = len(vec) // 2
    left = vec[:half]
    right = vec[half:]
    return [left[i] * val + right[i] * (1 / val) for i in range(half)]


def get_LR(a, b, Gs, Hs, Q):
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

def round(a, b, Gs, Hs, P, Q, Fr):
    L, R = get_LR(a, b, Gs, Hs, Q)
    x = Fr.random_element()
    a_prime = fold(a, x)
    b_prime = fold(b, 1 / x)
    G_prime = fold(Gs, 1 / x)
    H_prime = fold(Hs, x)
    
    P_prime = inner_product(a_prime, G_prime) + inner_product(b_prime ,H_prime) + inner_product(a_prime, b_prime) * Q
    assert P_prime == P + x ^ 2 * L + x ^ (-2) * R
    return x, L, R, a_prime, b_prime, G_prime, H_prime, P_prime

def ipa(a, b, Gs, Hs, c, Q, Fr):
    _c = inner_product(a, b)
    assert c == _c
    P = inner_product(a, Gs) + inner_product(b, Hs) + c * Q

    xs = []
    Ls = []
    Rs = []
    while len(a) > 1:
        x,L,R,a,b,Gs,Hs,P = round(a, b, Gs, Hs, P, Q, Fr)
        xs.append(x)
        Ls.append(L)
        Rs.append(R)
    # print("a = ", a)
    # print("b = ", b)
    # print("Gs = ", Gs)
    # print("Hs = ", Hs)
    # print("P = ", P)
    return xs, Ls, Rs, a, b

def folding_weights(xs, N, Fr):
    s_G = [Fr(1) for _ in range(N)]
    s_H = [Fr(1) for _ in range(N)]

    for i in range(N):
        i_bits = bin(i)[2:].zfill(int(math.log(N, 2)))
        for k in range(int(math.log(N, 2))):
            bit_j = int(i_bits[k])
            if bit_j == 1:
                s_G[i] *= xs[k]
                s_H[i] *= 1 / xs[k]
            else:
                s_G[i] *= 1 / xs[k]
                s_H[i] *= xs[k]

    return s_G, s_H

def verify(Gs, Hs, P, xs, Ls, Rs, a_final, b_final, Q, N, Fr):
    s_G, s_H = folding_weights(xs, N, Fr)
    # print("s_G = ", s_G)
    # print("s_H = ", s_H)

    G_final = inner_product(s_G, Gs)
    H_final = inner_product(s_H, Hs)
    # print("G_final = ", G_final)
    # print("H_final = ", H_final)

    P_final = inner_product(a_final, [G_final]) + inner_product(b_final, [H_final]) + inner_product(a_final, b_final) * Q
    # print("P_final = ", P_final)
    assert P == P_final - sum([xs[i] ^ 2 * Ls[i] + xs[i] ^ (-2) * Rs[i] for i in range(len(xs))])

def run_ipa():
    p = 929
    Fp = GF(p)
    E = EllipticCurve(Fp, [5, 15])
    G = E.gens()[0]
    Fr = GF(G.order())
    assert is_prime(G.order())

    print("We define Q as a random point on the curve")
    Q = E.random_point()

    # Length of the vectors (this needs to be even length)
    N = 8

    a = [Fr.random_element() for _ in range(N)]
    b = [Fr.random_element() for _ in range(N)]

    # random elements
    Gs = [E.random_point() for _ in range(N)]
    Hs = [E.random_point() for _ in range(N)]

    xs, Ls, Rs, a_final, b_final = ipa(a, b, Gs, Hs, inner_product(a, b), Q, Fr)

    # Verification
    P = inner_product(a, Gs) + inner_product(b, Hs) + inner_product(a, b) * Q
    verify(Gs, Hs, P, xs, Ls, Rs, a_final, b_final, Q, N, Fr)
    print("Verification successful!")
# run_ipa()