#!/usr/bin/env python3
"""
Python script to compute Montgomery constants for the scalar field.
Run with: sage -python compute_constants.py
Or just copy-paste into a SageMath notebook/terminal.
"""

from sage.all import *

# Field modulus from the curve specification
p = Integer(0xf06e44682c2aa440f5f26a5ae1748ff85ccc2efc3068faf2154ff8a2e94d81)

print("="*80)
print("Montgomery Constants Computation for Scalar Field")
print("="*80)
print()

# Verify p is prime
print(f"Field modulus p = 0x{int(p):x}")
print(f"p in decimal = {p}")
print(f"Bit length of p: {p.nbits()}")
print(f"Is p prime? {is_prime(p)}")
print()

# Montgomery parameter R = 2^256
R_exp = 256
R = Integer(2)**Integer(R_exp)
print(f"R = 2^{R_exp}")

# Compute R mod p
R_mod_p = R % p
print(f"R mod p = 0x{int(R_mod_p):x}")
print()

# Compute R^2 mod p (for converting to Montgomery form)
R2_mod_p = (R * R) % p
print(f"R^2 mod p = 0x{int(R2_mod_p):x}")
print()

# Compute R^3 mod p (sometimes useful for optimizations)
R3_mod_p = (R * R * R) % p
print(f"R^3 mod p = 0x{int(R3_mod_p):x}")
print()

# Compute MU = -p^{-1} mod 2^64
# This is used in Montgomery reduction
mod_2_64 = Integer(2)**Integer(64)
p_inv_mod_2_64 = inverse_mod(p, mod_2_64)
MU = (-p_inv_mod_2_64) % mod_2_64
print(f"MU = -p^{{-1}} mod 2^64 = 0x{int(MU):016x}")
print()

# Helper function to convert 256-bit number to [u64; 4] array in little-endian
def to_u64_array(n):
    """Convert a 256-bit number to [u64; 4] array in little-endian format."""
    n = int(n)  # Convert Sage Integer to Python int
    limbs = []
    for i in range(4):
        limb = (n >> (64 * i)) & ((1 << 64) - 1)
        limbs.append(limb)
    return limbs

# Convert to Rust array format
print("="*80)
print("Rust Constants (copy these into your scalarfield.rs file)")
print("="*80)
print()

print("// Field modulus: p = 0x{:x}".format(int(p)))
print("const MODULUS: [u64; 4] = [")
p_limbs = to_u64_array(p)
for i, limb in enumerate(p_limbs):
    print(f"    0x{limb:016x},")
print("];")
print()

print("// R = 2^256 mod p (Montgomery parameter)")
print("const R: [u64; 4] = [")
r_limbs = to_u64_array(R_mod_p)
for i, limb in enumerate(r_limbs):
    print(f"    0x{limb:016x},")
print("];")
print()

print("// R^2 = 2^512 mod p (for Montgomery conversion)")
print("const R2: [u64; 4] = [")
r2_limbs = to_u64_array(R2_mod_p)
for i, limb in enumerate(r2_limbs):
    print(f"    0x{limb:016x},")
print("];")
print()

print("// R^3 = 2^768 mod p (for efficient conversion)")
print("#[allow(dead_code)]")
print("const R3: [u64; 4] = [")
r3_limbs = to_u64_array(R3_mod_p)
for i, limb in enumerate(r3_limbs):
    print(f"    0x{limb:016x},")
print("];")
print()

print("// -p^{{-1}} mod 2^64 (Montgomery parameter mu)")
print(f"const MU: u64 = 0x{MU:016x};")
print()

# Verify the constants
print("="*80)
print("Verification")
print("="*80)
print()

# Verify R * R^{-1} ≡ 1 (mod p)
R_inv = inverse_mod(R, p)
print(f"R * R^{{-1}} mod p = {(R * R_inv) % p} (should be 1)")

# Verify R^2 calculation
print(f"R^2 mod p correct: {R2_mod_p == pow(R, 2, p)}")

# Verify MU calculation
verification = (MU * p) % mod_2_64
print(f"MU * p mod 2^64 = 0x{int(verification):016x}")
print(f"Should be 2^64 - 1 = 0x{int(mod_2_64-1):016x}")
print(f"MU verification: {verification == (mod_2_64 - 1)}")
print()

# Additional useful information
print("="*80)
print("Additional Information")
print("="*80)
print()

# Check if p is a safe prime or has special structure
print(f"p - 1 = 2^k * q where:")
p_minus_1 = p - 1
k = 0
q = p_minus_1
while q % 2 == 0:
    k += 1
    q = q // 2
print(f"  k = {k} (2-adicity)")
print(f"  q = 0x{q:x}")
print()

# Find a generator (useful for FFTs and other operations)
print("Finding a small generator of the multiplicative group...")
for g in range(2, 20):
    g_sage = Integer(g)
    if power_mod(g_sage, p-1, p) == 1 and power_mod(g_sage, (p-1)//2, p) != 1:
        print(f"Generator g = {g}")
        g_mont = (g_sage * R_mod_p) % p
        print(f"Generator in Montgomery form: 0x{int(g_mont):x}")
        g_mont_limbs = to_u64_array(g_mont)
        print(f"const GENERATOR: Self = ScalarField {{ limbs: [0x{g_mont_limbs[0]:016x}, 0x{g_mont_limbs[1]:016x}, 0x{g_mont_limbs[2]:016x}, 0x{g_mont_limbs[3]:016x}] }};")
        break
print()

print("="*80)
print("Script completed successfully!")
print("="*80)
