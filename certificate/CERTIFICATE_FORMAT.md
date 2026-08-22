# Human-readable Bernstein certificate format

## Purpose

`grunbaum_certificate.txt` is a language-independent representation of the
complete binary subdivision tree and the separating covector attached to every
terminal box. It contains proof data only:

- the coordinate bisected at each internal node; and
- five exact integer entries at each leaf.

The tensor coefficients, chart polynomials, Bernstein conversion, interval
control data, de Casteljau subdivisions, and inequalities are reconstructed by
the checker.

## Header

The first line is

```text
G4BCERT-TEXT-1 nodes=77713 leaves=38857 witness_scale=4503599627370496
```

Here `witness_scale=2^52` records how candidate directions were quantized. The
checker treats each displayed witness as an exact integer vector. A common
positive scaling would not change its validity.

## Preorder node records

Every subsequent nonempty line is one node in preorder, left child before
right child. There are two record types.

### Internal node

```text
split x1
```

The allowed coordinates are

```text
x1 x2 x3 y1 y2 y3
```

The current box is bisected at the midpoint of that coordinate. The left child
is the lower half and the right child is the upper half.

The root is `[0,1]^6` in Bernstein coordinates. It corresponds to the
quaternion chart `[-1,1]^6` by

```text
chart_coordinate = 2 * bernstein_coordinate - 1.
```

### Leaf

```text
leaf w0 w1 w2 w3 w4
```

The five decimal integers form an exact covector `w`. On that box the checker
proves

```text
w . C_alpha > 0
```

for all `5^6 = 15625` vector Bernstein control points `C_alpha`. Since every
value of the five-component polynomial map is a convex combination of these
control points, the map cannot equal zero on the box.

## Exact root control data

The checker performs the following operations independently for each of the
five components.

1. It derives the six-variable power-basis polynomial from `5A`, `B`, and the
   quaternionic frame formulas.
2. It converts the polynomial on `[-1,1]^6` to the common tensor-product
   Bernstein basis of degree four in every coordinate. The exact one-variable
   conversion matrix, multiplied by six, is

   ```text
    6 -6  6 -6  6
    6 -3  0  3 -6
    6  0 -2  0  6
    6  3  0 -3 -6
    6  6  6  6  6
   ```

   After six axes the common omitted denominator is `6^6`.
3. It divides every component by the maximum absolute value of that
   component's exact root Bernstein coefficients. This positive componentwise
   normalization does not change the common zero set.
4. It stores an outward interval enclosure on the fixed scale `2^52`:

   ```text
   L[c,alpha] <= 2^52 * C[c,alpha] <= U[c,alpha].
   ```

## Midpoint subdivision

At `split <axis>`, the checker applies degree-four midpoint de Casteljau
subdivision along that coordinate. The exact child coefficients are repeated
averages of the parent coefficients. Lower endpoints are rounded downward and
upper endpoints upward, so the exact child control points remain enclosed.

## Exact leaf inequality

For a leaf witness `w`, the rigorous lower bound for control index `alpha` is

```text
sum_c w[c] * (L[c,alpha] if w[c] >= 0 else U[c,alpha]).
```

The checker evaluates this with signed `__int128` arithmetic and accepts the
leaf only when the integer is strictly positive for every control index.

## Tree statistics

```text
nodes:          77713
internal nodes: 38856
leaves:         38857
maximum depth:  60
```

For a full binary tree, `nodes = 2*leaves - 1`.
