"""Compatibility shim for notebook - imports from the package."""

from novoboard.config import vocab, vocab_reverse, vocab_size  # noqa: F401

# Print statements for notebook compatibility (original behavior)
print("Training vocab_reverse ", vocab_reverse)
print("Training vocab ", vocab)
print("Training vocab_size ", vocab_size)
