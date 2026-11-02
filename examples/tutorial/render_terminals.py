#!/usr/bin/env python3
"""Render the tutorial's terminal transcripts (*.terminal.html) to PNG.

Needs Playwright with Chromium (`pip install playwright` and
`playwright install chromium`, or a preinstalled browser):

    python examples/tutorial/render_terminals.py docs/images/tutorial
"""

import argparse
from pathlib import Path

from playwright.sync_api import sync_playwright


def shrink(path):
    """Store a terminal screenshot as a 128-color palette PNG, if Pillow is present."""
    try:
        from PIL import Image
    except ImportError:
        return
    with Image.open(path) as image:
        quantized = image.convert("RGB").quantize(colors=128, dither=Image.Dither.NONE)
    quantized.save(path, optimize=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--keep-html", action="store_true")
    args = parser.parse_args()
    pages = sorted(args.folder.glob("*.terminal.html"))
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(
            viewport={"width": 980, "height": 400}, device_scale_factor=1.5
        )
        for source in pages:
            page.goto(source.resolve().as_uri())
            target = source.with_name(source.name.replace(".terminal.html", ".png"))
            page.locator(".window").screenshot(path=str(target))
            shrink(target)
            if not args.keep_html:
                source.unlink()
            print(target)
        browser.close()


if __name__ == "__main__":
    main()
