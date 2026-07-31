"""
beckmann-core: shared library for the Beckmann rearrangement prediction
products (beckmann-nbo, beckmann-pyscf). Every function here is
method-agnostic (works the same whether the downstream product uses
Gaussian/NBO7 or PySCF) and takes explicit paths/arguments -- no shared
filesystem conventions or DATA_INPUT/DATA_OUTPUT-style constants live here;
each product defines its own.
"""
