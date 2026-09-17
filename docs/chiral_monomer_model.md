# Chiral monomers: what changes in the model, term by term

Design note, 2026-09-17. Nothing here is implemented. Equations of section 2 are the code as it is
(`src/pair_potential.cpp`, `src/chain.cpp`, `src/input.h`).

## 1. Variables

    block i = 1..N (one block = 7 monomers)
    r_i            position                         b_i = r_{i+1} - r_i   bond vector
    s_i            state:  0 = coil,  +1 = R helix,  -1 = L helix
    theta_i        valence angle at block i (between -b_{i-1} and b_i)
    u_i            axis = unit vector along  b_{i-1}/|b_{i-1}| + b_i/|b_i|

    NEW, fixed for the whole run (quenched):
    x_k = +1 / -1  chirality of monomer k
    X_i = sum of the 7 x_k of block i     ->  X_i in {-7,-5,-3,-1,+1,+3,+5,+7}

## 2. The Hamiltonian now (achiral monomers)

    H0 =  sum_i  1/2 k_c (|b_i| - l_c)^2                 c = pair class of (s_i, s_i+1):   CC, CH, HH, RL
        + sum_i  1/2 kappa_t (theta_i - theta0_t)^2      t = triple class of (s_i-1, s_i, s_i+1)
        + sum_i  eps(s_i, s_i+1)                         HH: -J0   CH: +J1   RL: +J2   CC: 0
        + sum_i  E_helix |s_i|
        + sum_{j > i+1}  U_nb(i,j)

    U_nb(i,j) = U_iso(r)                                  if i or j is coil      (isotropic core)
              = hf_scale [ A I_a(r) + eps_s chi C I_c(r) ]   if both are helices

    r = |r_j - r_i|,  rh = (r_j - r_i)/r,  e1 = rh.u_i,  e2 = rh.u_j
    psi = atan2( (u_i x u_j).rh ,  u_i.u_j - e1 e2 )
    (labels 1,2 are swapped so that S = e1 + e2 >= 0; psi does not change under the swap)
    S = e1 + e2,  D = e1 - e2
    A = a0(S,D) + a1(S,D) cos psi + a2(S,D) cos 2psi          achiral sector  (even in psi)
    C = b1(S,D) sin psi + b2(S,D) sin 2psi                     chiral sector   (odd in psi)
    every coefficient is  q(S,D) = c0 + c1 S + c2 D + c3 S^2 + c4 D^2 + c5 S D
    two coefficient tables: HOMO (R.R and L.L) and RL
    chi = +1 for R.R,  -1 for L.L;   for R.L pairs chi = +1 and psi -> -psi when the first helix is L
    I(r) = sqrt(pi/2) sigma exp(r0^2 / 2 sigma^2) erfc( r / sqrt(2) sigma ),  r0 = 1.7 a, cutoff 3 a
    sigma_a = 0.93 a (t100), 0.87 a (t45);  sigma_c = 0.70 a / 0.51 a (medians)
    chiral sector only inside the fitted domain (|e| <= 0.5, r in [1.5, 2.5] a, shifted to 0 at 2.5 a)

Symmetries of H0:  rotations, translations, exchange 1<->2, and
mirror of all coordinates together with R<->L. The last one gives  P(M) = P(-M)  exactly,
U_RR(psi) = U_LL(-psi)  and  U_RL(e1,e2,psi) = U_RL(-e2,-e1,-psi) [exchange + mirror].

## 3. With chiral monomers -- the model AS BUILT

"As built" = the state of a block selects its fine-grained Hamiltonian (R block: dihedral well at +phi0,
L block: well at -phi0, coil block: K* = 0). This is how every table above was measured.
The monomer chirality then only changes the statistical weight of the two helical states:

    H = H0  -  h  sum_i  s_i X_i

    h = half of the free-energy difference, per monomer, between the monomer sitting in the R helix
        and in the L helix.  One new parameter.  Coil blocks (s = 0) feel nothing.

Every other term is EXACTLY unchanged: k_c, l_c, kappa_t, theta0_t, J0, J1, J2, E_helix, U_iso,
A, C, I_a, I_c, chi, both tables. Reason: an R block has the same shape whatever monomers built it,
and the pair potential is made by the shapes only.

Symmetry: the term -h s X is odd under R<->L and contains no coordinates. It breaks
"mirror + R<->L". What remains is "mirror + R<->L + (X -> -X)", which maps the chain onto its
mirror-image polymer, a different system. P(M) = P(-M) is lost. h X is the field conjugate to
m = (n_R - n_L)/N. The two pair-table relations stay exactly true.

## 4. Additional changes if the fine-grained model is ONE double-well dihedral potential

Fine-grained dihedral energy  V(phi) = V_sym(phi) - delta x_k g(phi),  g odd (e.g. sin phi),
well stiffness K*, small parameter  eps_d = delta / K*.

Expanding around the wells, phi = s phi0 + eta:

    constant :  - delta x s g(phi0)                  -> the on-site term of section 3, h = delta g(phi0)
    linear   :  well moves, |phi| = phi0 + Delta     Delta_i = eps_d g'(phi0) s_i X_i / 7
    quadratic:  stiffness  K* (1 - eps_d g''(phi0) s_i x_k)

Every table parameter p that depends on the helix shape then becomes (first order in eps_d):

    bonds      l_HH      -> l_HH      + (dl_HH/dphi0) (Delta_i + Delta_i+1)/2         same for k_HH, l_CH, l_RL ...
    bends      theta0_HHH-> theta0_HHH+ (dtheta0/dphi0) (mean Delta of the helical blocks of the triple)
    achiral    A         -> A + A'_1 Delta_i + A'_2 Delta_j                           A'_n = dA/dphi0 of helix n
    chiral     chi C     -> chi ( C + C'_1 Delta_i + C'_2 Delta_j )
               for a same-handed pair  chi Delta = eps_d g' X / 7 :  the new piece is ODD in X and
               multiplies the odd-in-psi sector. This is the only new CHIRAL piece of the pair potential.
    stiffness  all of the above also get a term  (dp/dK*) ( - eps_d g'' s X / 7 ) K*
    J0 J1 J2   unchanged at block resolution; at monomer resolution a domain wall prefers to sit
               where x changes sign (cannot be represented with blocks of 7)
    E_helix    -> E_helix + 7 (delta^2 / 2 kT) [ <g^2>_coil - var_well(g) ]             uniform, second order
    coil core  isotropic -> no axis -> no chiral sector is possible; size change is even in X, order delta^2
    x1 x2 P2   appears only at SECOND order and carries s1 s2:
               eps_d^2 g'^2 (s_i X_i)(s_j X_j)/49 * d2A/dphi0_1 dphi0_2

Symmetry relations in this picture:

    U_RR(X1,X2; e1,e2,psi) = U_LL(-X1,-X2; e1,e2,-psi)
    U_RL(X1,X2; e1,e2,psi) = U_RL(-X2,-X1; -e2,-e1,-psi)
    exchange 1<->2 with the labels carried along: exact as before

Size: eps_d = 0.01 for delta = 1 kT at K* = 100 (the stiffness of the umbrella campaign).
NOT KNOWN: none of the derivatives dp/dphi0, dp/dK* has been measured. The two campaigns (t100, t45)
differ in theta0, not in the dihedral. In this picture the corrections cannot be written as numbers
without a new umbrella measurement at shifted phi0.

### Limits of section 4 (what is rigorous and what is not)

Rigorous: (i) the symmetry relations; (ii) that a first-order term exists, is linear in each x_k and carries
the factor s (perturbation theory in delta, smooth dependence of the pair free energy on the block shape);
(iii) the helix geometry Omega, d, rho from theta0, phi0 (checked against an explicitly built chain).

Approximations made to reach the compact form with one A' and one C':
  - the SUM X_i is used. Exactly, the first-order change is  eps_d s_i sum_k x_k A'_k  with a different weight
    A'_k for each of the 7 monomer positions in the block (a dihedral in the middle of the block and one at its
    end deform the block differently). A' = mean of A'_k is exact only if all positions weigh the same.
    (For the on-site term the sum IS exact: every monomer contributes -h s x_k independently.)
  - dihedrals span 4 beads and straddle block boundaries; assigning each to one block is a choice.
  - the radial widths sigma_a, sigma_c also depend on the shape, so I_a and I_c change too (omitted above).
  - harmonic, stiff wells (K* = 100) are assumed for the size of the shift.
  - spherical beads: the monomer chirality acts through the dihedral only.

Not tested: no number in section 4 has been checked against a simulation. Test = the pair free energy of two
blocks with shifted phi0 (umbrella, or rigid 7-bead fragments as a first estimate).

## 5. Not covered

Monomers that are themselves chiral-shaped objects (not spheres): every bead-bead contact becomes chiral,
also for coil pairs. That needs a new fine-grained model and a new measurement; it cannot be written by hand.
