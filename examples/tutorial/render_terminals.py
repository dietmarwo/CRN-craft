#!/usr/bin/env python3
"""Render the tutorial's terminal transcripts (*.terminal.html) to PNG.

Needs Playwright with Chromium (`pip install playwright` and
`playwright install chromium`):

    python examples/tutorial/render_terminals.py docs/images/tutorial

To use a browser that is already installed, such as Chrome, name it with
`--browser PATH` or the environment variable `CRN_TUTORIAL_BROWSER`; Playwright
then launches that executable instead of its own download:

    python examples/tutorial/render_terminals.py docs/images/tutorial \\
        --browser /usr/bin/google-chrome
"""

import argparse
import os
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
    parser.add_argument(
        "--browser",
        default=os.environ.get("CRN_TUTORIAL_BROWSER"),
        help="path of an installed Chromium-based browser (default: the one "
        "Playwright downloaded; env CRN_TUTORIAL_BROWSER)",
    )
    args = parser.parse_args()
    pages = sorted(args.folder.glob("*.terminal.html"))
    with sync_playwright() as playwright:
        options = {} if args.browser is None else {"executable_path": args.browser}
        browser = playwright.chromium.launch(**options)
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
