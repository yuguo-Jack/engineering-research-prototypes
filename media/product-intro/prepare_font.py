"""Make the redistributable, renamed CJK font subset used by the film."""
from pathlib import Path
import argparse

from fontTools import subset
from fontTools.ttLib import TTFont

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Official OFL Noto Sans SC variable TTF")
    parser.add_argument("--license", type=Path, required=True, help="Its original SIL OFL 1.1 notice")
    args = parser.parse_args()
    text = (HERE / "index.html").read_text(encoding="utf-8")
    text += (HERE / "storyboard.json").read_text(encoding="utf-8")
    text += "".join(chr(value) for value in range(32, 127))
    options = subset.Options()
    options.flavor = "woff2"
    options.name_IDs = [0, 1, 2, 3, 4, 5, 6, 13, 14, 16, 17]
    font = TTFont(args.source)
    subsetter = subset.Subsetter(options=options)
    subsetter.populate(text=text)
    subsetter.subset(font)
    names = {1: "ERP Sans", 4: "ERP Sans", 6: "ERPSans", 16: "ERP Sans"}
    for record in font["name"].names:
        if record.nameID in names:
            record.string = names[record.nameID].encode(record.getEncoding())
    target = HERE / "fonts"
    target.mkdir(exist_ok=True)
    font.flavor = "woff2"
    output = target / "erp-sans.woff2"
    font.save(output)
    notice = "\n".join(line.rstrip() for line in args.license.read_text(encoding="utf-8").splitlines()) + "\n"
    (target / "OFL.txt").write_text(notice, encoding="utf-8", newline="\n")
    print(f"ERP Sans subset: {output.stat().st_size} bytes; source copyright/license retained")


if __name__ == "__main__":
    main()
