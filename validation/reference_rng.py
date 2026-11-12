"""Independent scalar ChaCha8, Bernstein layout (64-bit counter/stream).

Algorithm: https://cr.yp.to/chacha/chacha-20080128.pdf. No third-party code
copied. This intentionally uses scalar integer operations, independently of
rand_chacha/ppv-lite86. Not a production or cryptographic API.

Since 0.7.0 one engine draws every replicate from the stream of
`crn_craft::rules::stream`: ChaCha8 keyed by SHA-256 of the domain
`crn-craft.rules/0.1` (no terminating NUL) followed by the root seed as eight
little-endian bytes, with the replicate index as the 64-bit stream id, read
as little-endian 64-bit words from word position zero. `words(root,
replicate)` is that stream and `rng-vectors-rules.json` holds its known
answers (nine streams, 40 words each); the Rust test
`crates/crn-craft/tests/boundaries_v1.rs` checks the engine's stream against
the same file. The network engine's stream (`crn-craft.chacha8/0.1`) and its
`rng-vectors.json` were removed with that engine in task M8b.
"""

import hashlib
import json
from pathlib import Path
import struct

MASK = (1 << 32) - 1


def rotate(value, bits):
    return ((value << bits) | (value >> (32 - bits))) & MASK


def block(key, counter, stream):
    initial = list(struct.unpack("<4I", b"expand 32-byte k"))
    initial += list(struct.unpack("<8I", key))
    initial += [counter & MASK, counter >> 32, stream & MASK, stream >> 32]
    x = initial.copy()

    def quarter(a, b, c, d):
        x[a] = (x[a] + x[b]) & MASK
        x[d] = rotate(x[d] ^ x[a], 16)
        x[c] = (x[c] + x[d]) & MASK
        x[b] = rotate(x[b] ^ x[c], 12)
        x[a] = (x[a] + x[b]) & MASK
        x[d] = rotate(x[d] ^ x[a], 8)
        x[c] = (x[c] + x[d]) & MASK
        x[b] = rotate(x[b] ^ x[c], 7)

    for _ in range(4):
        for indices in [(0, 4, 8, 12), (1, 5, 9, 13), (2, 6, 10, 14), (3, 7, 11, 15)]:
            quarter(*indices)
        for indices in [(0, 5, 10, 15), (1, 6, 11, 12), (2, 7, 8, 13), (3, 4, 9, 14)]:
            quarter(*indices)
    words = [(a + b) & MASK for a, b in zip(x, initial, strict=True)]
    return [words[i] | (words[i + 1] << 32) for i in range(0, 16, 2)]


DOMAIN = b"crn-craft.rules/0.1"


def words(root, replicate, domain=DOMAIN):
    key = hashlib.sha256(domain + root.to_bytes(8, "little")).digest()
    counter = 0
    while True:
        yield from block(key, counter, replicate)
        counter += 1


def vectors(domain=DOMAIN):
    rows = []
    for root in [0, 42, 2026110401]:
        for replicate in [0, 1, (1 << 64) - 1]:
            stream = words(root, replicate, domain)
            rows.append(
                {
                    "root": root,
                    "replicate": replicate,
                    "words": [next(stream) for _ in range(40)],
                }
            )
    return {"schema": "crn-craft.rng-vectors/0.1", "streams": rows}


if __name__ == "__main__":
    here = Path(__file__).parent
    expected = json.loads((here / "rng-vectors-rules.json").read_text())
    if expected != vectors():
        raise RuntimeError("independent ChaCha8 known-answer fixture differs")
    print(
        json.dumps(
            {
                "schema": "crn-craft.rng-check/0.1",
                "streams": 9,
                "words": 360,
                "passed": True,
            }
        )
    )
