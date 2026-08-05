# Rank-30 dual-CAS continuation

This branch is an isolated computational continuation of the rank-30 elliptic-curve research pipeline. It does not change the profile README and does not claim a new record.

`payload.tar.gz.b64` decodes to a deterministic archive containing the exact curve/point data and separately written Python, SageMath and Magma verifiers. The decoded archive SHA-256 is:

```text
e373381129472924604689df1443d68dc2cb51c05390dc52193761e60523a28d
```

The workflow executes:

- the dependency-free finite-reduction independence certificate;
- SageMath 10.9 global-minimal-model, torsion, conductor, root-number and 29-by-29 height-matrix computations;
- SageMath saturation at all primes below `2^12`;
- exact and height-based checks submitted as separate programs to the official free Magma calculator.

Timeouts and unavailable proprietary computations are recorded as unresolved, never promoted to a pass.
