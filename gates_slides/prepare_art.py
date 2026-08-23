#!/usr/bin/env python3
"""Prepare the illustration panels used by the deck.

The supplied illustrations are light artwork on a solid black background,
which cannot sit on a light slide. Two conversions are applied:

  line art (molecule)  -> invert lightness only, keeping hue/saturation, so
                          light-grey bonds become dark grey while the red O
                          and blue H stay red and blue.
  filled art (landscape) -> flood-fill the OUTER black region only, so the
                          illustration's own dark base plate survives.

Crops of the optical-lattice and quantum-circuit figure are taken as-is
(already dark-on-white).

Run from gates_slides/:  python3 prepare_art.py
"""
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

SRC = Path("/home/fts/.cursor/projects/"
           "xuanwu-tank-east-fts-projects-transition-1x-analysis/assets")
OUT = Path(__file__).resolve().parent / "art"


def rgb_to_hls(rgb):
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    mx, mn = np.max(rgb, -1), np.min(rgb, -1)
    l = (mx + mn) / 2
    d = mx - mn
    s = np.zeros_like(l)
    nz = d > 1e-9
    denom = np.where(l > 0.5, 2.0 - mx - mn, mx + mn)
    s[nz] = d[nz] / np.maximum(denom[nz], 1e-9)
    h = np.zeros_like(l)
    ir, ig, ib = (mx == r) & nz, (mx == g) & nz, (mx == b) & nz
    h[ir] = ((g - b)[ir] / d[ir]) % 6
    h[ig] = ((b - r)[ig] / d[ig]) + 2
    h[ib] = ((r - g)[ib] / d[ib]) + 4
    return h / 6.0, l, s


def _hue(p, q, t):
    t = t % 1.0
    out = np.copy(p)
    m1 = t < 1 / 6
    m2 = (t >= 1 / 6) & (t < 1 / 2)
    m3 = (t >= 1 / 2) & (t < 2 / 3)
    out[m1] = p[m1] + (q[m1] - p[m1]) * 6 * t[m1]
    out[m2] = q[m2]
    out[m3] = p[m3] + (q[m3] - p[m3]) * (2 / 3 - t[m3]) * 6
    return out


def hls_to_rgb(h, l, s):
    q = np.where(l < 0.5, l * (1 + s), l + s - l * s)
    p = 2 * l - q
    out = np.stack([_hue(p, q, h + 1 / 3), _hue(p, q, h), _hue(p, q, h - 1 / 3)], -1)
    grey = s < 1e-9
    out[grey] = l[grey][..., None]
    return np.clip(out, 0, 1)


def invert_lightness(img: Image.Image) -> Image.Image:
    rgb = np.asarray(img.convert("RGB"), dtype=np.float64) / 255.0
    h, l, s = rgb_to_hls(rgb)
    return Image.fromarray((hls_to_rgb(h, 1.0 - l, s) * 255).astype(np.uint8))


def drop_outer_black(img: Image.Image, bg=(255, 255, 255), thresh=42):
    """Replace only the background region connected to the image border."""
    rgb = np.asarray(img.convert("RGB")).astype(np.int16)
    dark = rgb.max(axis=2) < thresh
    lab, n = ndimage.label(dark)
    border = set(lab[0].tolist()) | set(lab[-1].tolist()) \
        | set(lab[:, 0].tolist()) | set(lab[:, -1].tolist())
    border.discard(0)
    outer = np.isin(lab, list(border))
    out = np.asarray(img.convert("RGB")).copy()
    out[outer] = bg
    return Image.fromarray(out)


def white_to_alpha(img: Image.Image, thresh=246) -> Image.Image:
    """Make the outer white background transparent, so the artwork can sit
    directly on a tinted panel instead of inside a white box."""
    rgb = np.asarray(img.convert("RGB"))
    light = rgb.min(axis=2) >= thresh
    lab, _ = ndimage.label(light)
    border = set(lab[0].tolist()) | set(lab[-1].tolist()) \
        | set(lab[:, 0].tolist()) | set(lab[:, -1].tolist())
    border.discard(0)
    outer = np.isin(lab, list(border))
    out = np.dstack([rgb, np.where(outer, 0, 255).astype(np.uint8)])
    return Image.fromarray(out, "RGBA")


def trim(img: Image.Image, bg=255, pad=14) -> Image.Image:
    a = np.asarray(img.convert("RGB")).astype(np.int16)
    content = np.abs(a - bg).max(axis=2) > 12
    ys, xs = np.where(content)
    if len(xs) == 0:
        return img
    box = (max(int(xs.min()) - pad, 0), max(int(ys.min()) - pad, 0),
           min(int(xs.max()) + pad, img.width), min(int(ys.max()) + pad, img.height))
    return img.crop(box)


def pad_to_ratio(img: Image.Image, ratio=4 / 3, bg=(255, 255, 255)):
    w, h = img.size
    tw, th = (w, int(round(w / ratio))) if w / h > ratio else (int(round(h * ratio)), h)
    canvas = Image.new("RGB", (max(tw, w), max(th, h)), bg)
    canvas.paste(img, ((canvas.width - w) // 2, (canvas.height - h) // 2))
    return canvas


def stack(top: Image.Image, bottom: Image.Image, gap=52, bg=(255, 255, 255)):
    """Stack two trimmed elements centred in one panel image."""
    w = max(top.width, bottom.width)
    canvas = Image.new("RGB", (w, top.height + gap + bottom.height), bg)
    canvas.paste(top, ((w - top.width) // 2, 0))
    canvas.paste(bottom, ((w - bottom.width) // 2, top.height + gap))
    return canvas


def trim_alpha(img: Image.Image, pad=24) -> Image.Image:
    """Crop a transparent-background PNG to its visible content."""
    img = img.convert("RGBA")
    box = img.getchannel("A").getbbox()
    if box is None:
        return img
    x0, y0, x1, y1 = box
    return img.crop((max(x0 - pad, 0), max(y0 - pad, 0),
                     min(x1 + pad, img.width), min(y1 + pad, img.height)))


def main():
    OUT.mkdir(exist_ok=True)

    # supplied transparent artwork: use as-is, just trimmed
    supplied = Path(__file__).resolve().parent / "image_sources_gates_presentation"
    for name in ("quantum_branch.png", "ai_branch.png", "pic1.png"):
        src = supplied / name
        if src.exists():
            out = trim_alpha(Image.open(src))
            out.save(OUT / name)
            print(name, out.size)

    yizhou = Image.open(SRC / "image-3066c19f-2451-49d9-b5c0-0542a98220be.png")
    lattice = pad_to_ratio(trim(yizhou.crop((0, 10, 600, 370))))
    white_to_alpha(lattice).save(OUT / "quantum_lattice.png")

    circuit = trim(yizhou.crop((40, 380, 640, 540)))
    white_to_alpha(circuit).save(OUT / "quantum_circuit.png")

    molecule = invert_lightness(Image.open(SRC / "pic1-a01292c3-f2ac-4a5a-9973-020742ed191b.png"))
    white_to_alpha(pad_to_ratio(trim(molecule))).save(OUT / "molecule.png")

    landscape = drop_outer_black(
        Image.open(SRC / "pic2-f99fdbdb-fb78-434f-8708-f034bcc7eee3.png"))
    white_to_alpha(pad_to_ratio(trim(landscape))).save(OUT / "landscape.png")

    for f in sorted(OUT.glob("*.png")):
        print(f.name, Image.open(f).size)


if __name__ == "__main__":
    main()
