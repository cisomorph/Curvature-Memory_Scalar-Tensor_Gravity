"""Make the shared, audited modules importable in every test."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'Ordered_Simulations', 'common'))
