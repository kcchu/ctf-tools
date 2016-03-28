# RSA Utilities
#
# Other factorization algorithms included with Sage:
# qsieve - The best algorithm for factoring numbers of the form \(pq\) up to
#          around 100 digits.
# ecm.factor - The best algorithm for factoring numbers of the form \(n=pm\),
#              where \(p\) is not “too big”.

def rsa_make_d(p, q, e):
    """Computes private exponent d given p, q and e.
    """
    return inverse_mod(e, (p - 1) * (q - 1))

def factor_rsa_modulus_n(N, e, d):
    """Factorize the RSA modulus N given the public and private exponents, e
    and d.

    Source: https://crypto.stanford.edu/~dabo/papers/RSA-survey.pdf
    CTF: BCTF 2016 crypto 500 Hyper RSA
    """
    k = d * e - 1
    while True:
        g = randint(2, N-1)
        t = k
        while t % 2 == 0:
            t = t // 2
            x = power_mod(g, t, N)
            y = gcd(x - 1, N)
            if x > 1 and y > 1:
                p = y
                q = N // y
                if p > q:
                    return (p, q)
                else:
                    return (q, p)

def factor_fermat(N):
    """Factorize N into p and q for p and q share half of their leading bits.
    i.e., if the gap between p and q is below the square root of p

    Source: http://facthacks.cr.yp.to/fermat.html
    CTF: BKP CTF 2016 Bob's Hat
    """
    if N <= 0: return [N]
    if is_even(N): return [2,N/2]
    a = ceil(sqrt(N))
    while not is_square(a^2 - N):
        a = a + 1
    b = sqrt(a^2-N)
    return [a - b,a + b]

def factor_lattice(N, nearp, howclose, t, k):
    """Finds p very quickly if p has about half as many bits as N and half of
    the leading bits of p are known.

    Source: http://facthacks.cr.yp.to/lattice.html
    """
    R.<x> = PolynomialRing(ZZ)
    f = howclose * x + nearp
    M = matrix(t)
    for i in range(t):
        M[i] = (f ^ i * N ^ max(k - i, 0)).coefficients(sparse=False) + [0] * (t - 1 - i)
    M = M.LLL()
    Q = sum(z * (x / howclose) ^ i for i, z in enumerate(M[0]))
    for r, multiplicty in Q.roots():
        if nearp + r > 0:
            g = gcd(N, nearp + r)
            if g > 1: return [g, N / g]
    return [1, N]

def factor_rsa_wiener(N, e):
    """Wiener's attack: Factorize the RSA modulus N given the public exponents
    e when d is small.

    Source: https://crypto.stanford.edu/~dabo/papers/RSA-survey.pdf
    CTF: BKP CTF 2016 Bob's Hat
    """
    N = Integer(N)
    e = Integer(e)
    cf = (e / N).continued_fraction().convergents()
    for f in cf:
        k = f.numer()
        d = f.denom()
        if k == 0:
            continue
        phi_N = ((e * d) - 1) / k
        b = -(N - phi_N + 1)
        dis = b ^ 2 - 4 * N
        if dis.sign() == 1:
            dis_sqrt = sqrt(dis)
            p = (-b + dis_sqrt) / 2
            q = (-b - dis_sqrt) / 2
            if p.is_integer() and q.is_integer() and (p * q) % N == 0:
                p = p % N
                q = q % N
                if p > q:
                    return (p, q)
                else:
                    return (q, p)

def modular_sqrt(a, p):
    """Returns modular square root of an integer number a modulo a prime number
    p.

    Source: http://www.mersennewiki.org/index.php/Modular_Square_Root
    """

    if a == 0:
        return (0, )

    if legendre_symbol(a, p) != 1:
        raise ValueError("a is not a quadratic residue modulo p")

    if p == 2:
        return (a % p,)

    if p % 4 == 3:
        r = power_mod(a, (p + 1) // 4, p)
        return (r, p - r)

    if p % 8 == 5:
        v = power_mod(2 * a, (p - 5) // 8, p)
        i = (2 * a * power_mod(v, 2, p)) % p
        r = a * v * (i - 1) % p
        return (r, p - r)

    if p % 8 == 1: # Shanks' method
        q = p - 1
        e = 0
        while q % 2 == 0:
            q //= 2
            e += 1
        while True:
            x = Integer(randint(2, p - 1))
            z = power_mod(x, q, p)
            if power_mod(z ^ 2, e - 1, p) != 1:
                break
        y = z
        r = e
        x = power_mod(a, (q - 1) // 2, p)
        v = (a * x) % p
        w = (v * x) % p
        while w != 1:
            k = 1
            ws = power_mod(w, 2, p)
            while power_mod(ws, k, p) != 1:
                k += 1
            d = power_mod(y ^ 2, r - k - 1, p)
            y = power_mod(d, 2, p)
            r = k
            v = (d * v) % p
            w = (w * y) % p
        return (v, p - v)

    raise ValueError("Cannot find modular square root")

def modular_sqrt_rabin(a, p, q):
    """Returns modular square root of an integer number a modulo the product of
    two primes p and q. (i.e. decryption of Rabin)

    CTF: HITCON 2015 Quals crypto 314 Rsabin
    """
    rp = modular_sqrt(a, p)
    rq = modular_sqrt(a, q)
    r = []
    for a in rp:
        for b in rq:
            r.append(crt(a, b, p, q))
    return r
