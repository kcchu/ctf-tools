# RSA Utilities
#
# Other factorization algorithms included with Sage:
# qsieve - The best algorithm for factoring numbers of the form \(pq\) up to
#          around 100 digits.
# ecm.factor - The best algorithm for factoring numbers of the form \(n=pm\),
#              where \(p\) is not “too big”.

def factor_rsa_modulus_n(N, e, d):
    """Factorize the RSA modulus N given the public and private exponents, e
    and d.

    Source: https://crypto.stanford.edu/~dabo/papers/RSA-survey.pdf
    Reference: BCTF 2016 crypto 500 Hyper RSA
    """
    k = d * e - 1
    while True:
        g = randint(2, N-1)
        t = k
        while t % 2 == 0:
            t = Integer(t / 2)
            x = power_mod(g, t, N)
            y = gcd(x - 1, N)
            if x > 1 and y > 1:
                p = y
                q = Integer(N / y)
                if p > q:
                    return (p, q)
                else:
                    return (q, p)

def factor_fermat(N):
    """Factorize N into p and q for p and q share half of their leading bits.
    i.e., if the gap between p and q is below the square root of p

    Source: http://facthacks.cr.yp.to/fermat.html
    Reference: BKP CTF 2016 Bob's Hat
    """
    if N <= 0: return [N]
    if is_even(N): return [2,N/2]
    a = ceil(sqrt(N))
    while not is_square(a^2 - N):
        a = a + 1
    b = sqrt(a^2-N)
    return [a - b,a + b]
