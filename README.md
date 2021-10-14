# VMEC_9011

This is one of the first versions of VMEC, denoted `VMEC9011`.

## Notes

From documentation in `vmec.f`:

> THE POLOIDAL ANGLE IS DETERMINED BY MINIMIZING \<M\> = m\*\*2 S(m) , WHERE S(m) = Rm\*\*2 + Zm\*\*2.

This is different from what is published in Hirshman & Meier (1985)
and could explain why `xmpq(m,1) = m(m-1)` instead of `m^p(m^q-<M>)`.
