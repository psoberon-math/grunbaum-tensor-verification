# Format of `grunbaum_leaf_boxes.csv`

This CSV is a redundant, explicit rendering of the terminal boxes in
`grunbaum_certificate.txt`. It is intended to make independent implementations
possible without parsing a preorder tree.

There are 38,857 data rows, one per leaf. The header is

```text
leaf_id,node_id,depth,
x1_k,x1_d,x2_k,x2_d,x3_k,x3_d,
y1_k,y1_d,y2_k,y2_d,y3_k,y3_d,
w0,w1,w2,w3,w4
```

For each coordinate, the pair `(k,d)` denotes the exact dyadic interval

```text
[k / 2^d, (k+1) / 2^d]
```

in the Bernstein cube `[0,1]^6`. To recover the corresponding interval in the
quaternion chart `[-1,1]^6`, apply `x -> 2x-1` to both endpoints.

The last five columns are the exact integer separating covector. Their
semantics and the exact root normalization are specified in
`CERTIFICATE_FORMAT.md`.

The boxes have disjoint interiors and their exact dyadic volumes sum to one.
The CSV was regenerated from the text certificate and compared byte for byte
with the copy at the repository root.
