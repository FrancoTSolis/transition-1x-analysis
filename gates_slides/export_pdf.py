#!/usr/bin/env python3
"""Export the HTML deck to a vector PDF (13.333 x 7.5 in, 16:9).

Chromium prints the page rather than screenshotting it, so all text stays
real text (selectable, searchable, sharp at any zoom) and the SVG diagrams
stay vector. Only the four embedded illustrations remain raster, since they
are bitmaps to begin with.

Run from gates_slides/:  python3 export_pdf.py  ->  quantum_ai_codesign.pdf
"""
from pathlib import Path

from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
OUT = HERE / "quantum_ai_codesign.pdf"


def main():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto((HERE / "index.html").as_uri())
        page.emulate_media(media="print")
        page.wait_for_timeout(400)
        page.pdf(path=str(OUT), width="13.333in", height="7.5in",
                 print_background=True, margin={"top": "0", "bottom": "0",
                                                "left": "0", "right": "0"})
        browser.close()
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
