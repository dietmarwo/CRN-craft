"""Extract the eLife source tables and HB observations without spreadsheet software."""

import argparse
import csv
import io
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4486817/fullTextXML"
HB = "https://cdn.elifesciences.org/articles/07464/elife-07464-fig4-data1-v3.zip"
TABLE = "{urn:oasis:names:tc:opendocument:xmlns:table:1.0}"
TEXT = "{urn:oasis:names:tc:opendocument:xmlns:text:1.0}"
OFFICE = "{urn:oasis:names:tc:opendocument:xmlns:office:1.0}"


def read(url, local):
    if local:
        return local.read_bytes()
    with urlopen(url, timeout=35) as response:
        return response.read()


def table1(raw):
    root = ET.fromstring(raw)
    table = root.find(".//table-wrap[@id='tbl1']")
    kind = None
    for row in table.iter("tr"):
        cells = ["".join(td.itertext()).strip() for td in row]
        if not cells:
            continue
        if cells[0] == "Time/dpc":
            kind = "copies" if cells[1] == "E(m)" else "normalized_variance"
        elif len(cells) == 4 and cells[1] and kind:
            source = re.match(r"[A-Za-z]+", cells[3]).group()
            age = 100 if "*" in cells[0] else float(cells[0])
            # Table 1 labels these early points by division, then rounds
            # their age to two decimals. Use the end of that same division.
            if "†" in cells[3]:
                nearest = round(age * 24 / 7) * 7 / 24
                age = nearest if abs(age - nearest) <= 0.0051 else age
            n = re.match(r"[0-9]+", cells[2])
            yield kind, age, float(cells[1]), int(n.group()) if n else "", source


def hb_rows(raw):
    with ZipFile(io.BytesIO(raw)) as bundle:
        byte_sheet = bundle.read("newheteroplasmymeasurements.ods")
    with ZipFile(io.BytesIO(byte_sheet)) as ods:
        doc = ET.fromstring(ods.read("content.xml"))
    rows = []
    for row in doc.find(f".//{TABLE}table").findall(f"{TABLE}table-row"):
        cells = []
        for c in row.findall(f"{TABLE}table-cell"):
            val = c.get(f"{OFFICE}value") or " ".join(
                "".join(t.itertext()) for t in c.findall(f"{TEXT}p")
            )
            cells.extend(
                [val] * min(16, int(c.get(f"{TABLE}number-columns-repeated", 1)))
            )
        rows.append(cells[:16] + [""] * max(0, 16 - len(cells)))
    headers = {r[0]: r for r in rows if r[0]}
    start = next(i for i, r in enumerate(rows) if r[0].startswith("Individual"))
    groups, observations = [], []
    for i in range(1, 16):
        label = headers["Reference"][i]
        litter = re.search(r"\b\d{4}\b", label).group()
        age = int(headers["Age / dpc"][i])
        n = int(headers["n"][i])
        groups.append(
            (
                i,
                litter,
                age,
                n,
                headers["E(h)"][i],
                headers["V(h)"][i],
                headers["V'(h)"][i],
            )
        )
        observations.extend(
            (i, litter, age, j + 1, float(row[i]) / 100)
            for j, row in enumerate(rows[start : start + n])
        )
    return groups, observations


def convert(xml=None, hb=None):
    groups, cells = hb_rows(read(HB, hb))
    for name, header, rows in [
        (
            "reference-table1.csv",
            "observable,day,value,n,study",
            table1(read(XML, xml)),
        ),
        (
            "reference-table2.csv",
            "group,litter,day,n,mean,variance,normalized_variance",
            groups,
        ),
        ("reference-hb.csv", "group,litter,day,cell,heteroplasmy", cells),
    ]:
        with (ROOT / name).open("w", newline="") as file:
            out = csv.writer(file, lineterminator="\n")
            out.writerow(header.split(","))
            out.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--xml", type=Path)
    p.add_argument("--hb", type=Path)
    a = p.parse_args()
    convert(a.xml, a.hb)
