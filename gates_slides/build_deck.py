#!/usr/bin/env python3
"""Assemble the self-contained two-slide Gates deck.

Inlines every illustration as base64 so index.html can be emailed or opened
anywhere. Art comes from prepare_art.py (illustrations) and make_results.py
(result panels computed from the trained model).

Run from gates_slides/:  python3 build_deck.py  ->  index.html
"""
import base64
from pathlib import Path

HERE = Path(__file__).resolve().parent
ART = HERE / "art"

IMAGES = {
    "{{QUANTUM_BRANCH}}": "quantum_branch.png",
    "{{AI_BRANCH}}": "ai_arm_v2.png",
    "{{MOLECULE}}": "pic1.png",
    "{{PARITY}}": "parity.png",
}


def b64(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()


def main():
    html = (HERE / "template.html").read_text()
    for token, name in IMAGES.items():
        if token in html:
            html = html.replace(token, b64(ART / name))
    out = HERE / "index.html"
    out.write_text(html)
    print(f"wrote {out} ({out.stat().st_size/1e6:.2f} MB, fully self-contained)")


if __name__ == "__main__":
    main()
