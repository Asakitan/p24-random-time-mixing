# Exact uniformly mixing random readout on the 24-site path

**Byline:** Chain.Ren · Independent research

**Version:** [proof-v1.0.1](https://github.com/Asakitan/p24-random-time-mixing/releases/tag/proof-v1.0.1) · 2026-11-27 (Australia/Sydney). OpenAI Codex contributed to the derivation, exposition and computation.

**Theorem.** Let $A$ be the adjacency Hamiltonian of the path on 24 sites with unit coupling, and let $M(t)_{ij}=|(\exp(itA))_{ij}|^2$. There exist at most 145 nonnegative real readout times and probabilities on them such that $\mathbb E M(t)=J/24$. A single law works for all initial sites. The times are obtained by existence, not tabulated; the explicit discrete certificate below lives on a phase torus rather than on the real time axis.

**Problem.** [Kang, arXiv:2609.14463v1, Theorem and Outlook](https://arxiv.org/html/2609.14463v1) gives existence for $q=N+1$ a prime, twice a prime, or a power of two; its negative families are multiples of 15 or 21, or multiples of 6 with $N\ge11$. The paper leaves the other moduli unresolved there. The present result handles $q=25$ only. See the [search record](SOURCES.md).

## 1. Sufficient cosine moments

Write $\theta_r=2\cos(r\pi/25)$, $u_r(k)=\sqrt{2/25}\sin(rk\pi/25)$, $E_r=u_ru_r^\top$ for $1\le k,r\le24$, and let $R$ reverse the site order. These are the path adjacency eigenvalues and projectors, by the sine addition formula and orthogonality. Let $\circ$ denote entrywise multiplication and $D_0=\sum_r E_r\circ E_r$.

For $i,j$ in the site range, use $\sin^2 v=(1-\cos 2v)/2$, the product-to-sum identity and $\sum_{r=1}^{24}\cos(2\pi mr/25)=25[25\mid m]-1$. Then
$$
(D_0)_{ij}=\frac{2+[i=j]+[i+j=25]}{50}.
$$
Since $Ru_r=(-1)^{r+1}u_r$, the sums of the odd and even projectors are $(I+R)/2$ and $(I-R)/2$. Their entrywise squares sum to $(I+R)/2$; call this last matrix $S$. Equivalently, $S=\sum_{r\equiv s\bmod2}E_r\circ E_s$, including the diagonal pairs.

Set $d_{rs}=-1/24$ for distinct modes of the same index parity and $d_{rs}=0$ for modes of opposite parity. If the expected cosines of all mode differences equal $d_{rs}$, expanding $\exp(itA)=\sum_r e^{it\theta_r}E_r$ and combining conjugate terms gives
$$
\mathbb E M(t)=D_0-(S-D_0)/24=J/24.
$$
Only sufficiency is needed. This uses the parity moments of Kang's path analysis.

## 2. Exact frequency coordinates

Put $z=e^{i\pi/25}$, $x=z+z^{-1}$ and $\beta=(1,x,\ldots,x^9)$. Since $z$ is a primitive 50th root and $[\mathbb Q(z):\mathbb Q]=\varphi(50)=20$, while $z$ has degree two over the real field $\mathbb Q(x)$ by $z^2-xz+1=0$ and nonreality, $x$ has degree 10. Its minimal polynomial is
$$
F(X)=X^{10}-10X^8+35X^6-X^5-50X^4+5X^3+25X^2-5X-1.
$$
Indeed $z^{10}F(z+z^{-1})=\Phi_{50}(z)=z^{20}-z^{15}+z^{10}-z^5+1$ as a Laurent polynomial identity before evaluation. The ten coordinates of $\beta$ are rationally independent. Starting with $H_0=2,H_1=X$, compute $H_{r+1}=XH_r-H_{r-1}$ modulo $F$, and let $v_r$ be the ten coefficients in ascending order. Then $v_r\cdot\beta=H_r(x)=z^r+z^{-r}=\theta_r$.

For $r>s$, reduce $v_r-v_s$ modulo overall sign by taking the lexicographically smaller of the vector and its negative. The verifier enumerates the 276 pairs, giving 144 distinct representatives $W$, with consistent target $d_w$ on every duplicate. Using signs in this way does not change a cosine.

## 3. Certificate of an interior point

For $g\in\{0,1,\ldots,5\}^{10}$ define
$$
B_w(g)=2\cos(2\pi w\cdot g/6).
$$
Evaluation needs only the integer table $(2,1,-1,-2,-1,1)$. [certificate.json](certificate.json) lists 145 such points $g_j$, positive integer numerators $m_j$, and denominator $D=\sum m_j$. Form $C$ with columns $(1,(B_w(g_j))_{w\in W})$, sorting $W$ lexicographically. [verify.py](verify.py) checks by direct integer arithmetic
$$
12\sum_j m_j B_w(g_j)=
\begin{cases}-D&\text{equal index parity},\\0&\text{opposite index parity}.\end{cases}
$$
It also computes $\det C\bmod 65537=33681$ by modular elimination, so $C$ is nonsingular over the reals. The displayed equations divide by $24D$ to give the desired expectations of the cosines. All 145 coefficients $m_j/D$ are strictly positive. The unscaled moment columns are thus the vertices of a full-dimensional simplex in $\mathbb R^{144}$ containing $d$ in its interior. A floating-point optimizer helped select points, but neither rounding its weights nor trusting its output is a premise of the proof.

## 4. From the phase simplex to actual readout times

For any $\phi\in(\mathbb R/2\pi\mathbb Z)^{10}$ let $b(\phi)=(\cos(w\cdot\phi))_{w\in W}$. By the continuous-time Kronecker theorem and rational independence of $\beta$, $t\beta\bmod2\pi$ for $t\ge0$ is dense in the full ten-torus. In particular each vertex $b(2\pi g_j/6)$ can be approximated arbitrarily closely by $b(t_j\beta)$, with its own nonnegative $t_j$.

Adjoin a top row of ones to the unscaled vertex matrix. Its inverse exists, and its inverse applied to $(1,d)$ has only strictly positive coordinates. Inversion and matrix-vector multiplication are continuous on nonsingular matrices. Therefore all sufficiently small simultaneous vertex perturbations keep the matrix invertible and these coordinates positive. Pick the 145 approximating times inside such a neighborhood. With the *actual*, perturbed matrix, those inverse coordinates are real probabilities summing to one and give exactly $d$, not just an approximation. Section 1 now implies the theorem. No upper bound on the largest time or minimum-support assertion is involved. ∎

**Verification.** Run `python -B verify.py` with ordinary Python 3.9+ (standard library, assertions enabled). The script reconstructs the algebraic frequencies, all targets and the complete integer certificate.
