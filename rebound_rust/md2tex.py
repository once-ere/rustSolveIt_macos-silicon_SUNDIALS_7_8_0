#!/usr/bin/env python3
"""md2tex.py - convert rebound_rust.md to rebound_rust.tex.

This handles exactly the Markdown subset used by the port documentation:
ATX headings, fenced code blocks, GitHub pipe tables, bullet and numbered
lists, block quotes, horizontal rules, and the inline spans **bold**,
*italic*, `code` and [text](url).

It is deliberately a small, readable script rather than a dependency on
pandoc, so that regenerating the PDF needs nothing beyond Python and a
LaTeX installation.

Usage:  python md2tex.py rebound_rust.md rebound_rust.tex

Part of the rebound_rs documentation toolchain, GPL-3.0-or-later.
"""

import io
import re
import sys
import datetime

# ---------------------------------------------------------------- preamble

PREAMBLE = r"""\documentclass[11pt]{article}
\usepackage[margin=2.4cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{textcomp}
\usepackage{microtype}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{ragged2e}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{enumitem}
\usepackage[hidelinks]{hyperref}
\usepackage{parskip}

\lstset{
  basicstyle=\ttfamily\footnotesize,
  breaklines=true,
  breakatwhitespace=false,
  columns=fullflexible,
  keepspaces=true,
  frame=single,
  framerule=0.3pt,
  rulecolor=\color{black!30},
  backgroundcolor=\color{black!4},
  aboveskip=6pt, belowskip=6pt,
  literate={-}{{-}}1 {~}{{\textasciitilde}}1,
}

\newcommand{\code}[1]{\texttt{\small #1}}

% Table columns that wrap instead of overflowing the page.
\newcolumntype{L}[1]{>{\RaggedRight\arraybackslash}p{#1}}

\setlist[itemize]{leftmargin=1.4em, itemsep=1pt, topsep=3pt}
\setlist[enumerate]{leftmargin=1.7em, itemsep=1pt, topsep=3pt}

% The section titles already carry their own numbers ("15.9 REBOUNDx binary
% files..."), and the prose cross-references those numbers. Letting LaTeX add
% a second, different set would produce "16 15.9 REBOUNDx binary files" and
% make every cross-reference wrong. So headings are numbered by the document,
% not by LaTeX; they still appear in the table of contents.
\setcounter{secnumdepth}{-1}
\setcounter{tocdepth}{2}

\title{\textbf{rebound\_rust} \\[2pt] \large Provenance of the Pure-Rust Ports of
REBOUND 5.1.1 and REBOUNDx 5.1.0}
\author{rustSolveIt / Windows~11 port verification record}
\date{DATE_PLACEHOLDER}

\begin{document}
\maketitle
\tableofcontents
\clearpage
"""

POSTAMBLE = "\n\\end{document}\n"

# ------------------------------------------------------------- unicode map

UNICODE = {
    "—": "---",            # em dash
    "–": "--",             # en dash
    "‘": "`", "’": "'",
    "“": "``", "”": "''",
    "…": "\\ldots{}",
    "§": "\\S{}",
    "×": "\\(\\times\\)",
    "≤": "\\(\\leq\\)",
    "≥": "\\(\\geq\\)",
    "→": "\\(\\rightarrow\\)",
    "←": "\\(\\leftarrow\\)",
    "≈": "\\(\\approx\\)",
    "±": "\\(\\pm\\)",
    "µ": "\\(\\mu\\)",
    "φ": "\\(\\varphi\\)",
    "λ": "\\(\\lambda\\)",
    "Ω": "\\(\\Omega\\)",
    "ω": "\\(\\omega\\)",
    "°": "\\(^\\circ\\)",
    " ": "~",
    "−": "\\(-\\)",
    "✓": "yes",
    "✗": "no",
    "é": "\\'e",
    "ü": '\\"u',
    "ö": '\\"o',
    "ä": '\\"a',
    "©": "\\textcopyright{}",
    "ý": "\\'y",
    "ø": "\\o{}",
    "Š": "\\v{S}",
    "·": "\\(\\cdot\\)",
    "⚠": "\\textbf{!}",
    "⁻": "\\(^{-}\\)",
    "¹": "\\(^{1}\\)",
    "²": "\\(^{2}\\)",
    "³": "\\(^{3}\\)",
    "½": "\\(\\tfrac{1}{2}\\)",
    "α": "\\(\\alpha\\)",
    "β": "\\(\\beta\\)",
    "γ": "\\(\\gamma\\)",
    "Δ": "\\(\\Delta\\)",
    "δ": "\\(\\delta\\)",
    "ε": "\\(\\varepsilon\\)",
    "π": "\\(\\pi\\)",
    "σ": "\\(\\sigma\\)",
    "τ": "\\(\\tau\\)",
    "θ": "\\(\\theta\\)",
    "∞": "\\(\\infty\\)",
    "≠": "\\(\\neq\\)",
    "∼": "\\(\\sim\\)",
    "•": "\\textbullet{}",
}

SPECIALS = {
    "\\": r"\textbackslash{}",
    "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
    "_": r"\_", "{": r"\{", "}": r"\}",
    "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
}


# Plain-ASCII stand-ins for use inside verbatim listings, where a LaTeX
# command would be printed literally instead of being interpreted.
ASCII_FOLD = {
    "—": "--", "–": "-", "−": "-", "·": ".", "×": "x",
    "‘": "'", "’": "'", "“": '"', "”": '"',
    "…": "...", "§": "section ", "≤": "<=", "≥": ">=",
    "→": "->", "←": "<-", "≈": "~", "±": "+/-", "≠": "!=",
    "∞": "inf", "°": " deg", "•": "*", "✓": "yes", "✗": "no",
    "⚠": "!", "¹": "1", "²": "2", "³": "3", "⁻": "-",
}


def ascii_fold(text):
    for u, r in ASCII_FOLD.items():
        text = text.replace(u, r)
    return "".join(ch if ord(ch) < 128 else "?" for ch in text)


def esc(text):
    """Escape LaTeX specials in ordinary prose.

    Order matters: the specials are escaped FIRST and the Unicode
    substitutions applied afterwards. Doing it the other way round would
    feed the backslashes and braces of the substitutions (\\(\\times\\),
    \\S{}, ...) straight back into the special-character escaper and
    print them literally.
    """
    out = []
    for ch in text:
        out.append(SPECIALS.get(ch, ch))
    text = "".join(out)
    for u, r in UNICODE.items():
        text = text.replace(u, r)
    return text


def esc_code(text):
    """Escape for use inside \\texttt{} - keeps characters visible."""
    out = []
    for ch in text:
        out.append(SPECIALS.get(ch, ch))
    text = "".join(out)
    for u, r in UNICODE.items():
        # inside code, a dash is just a dash
        text = text.replace(u, "-" if u in ("—", "–") else r)
    return text


# One pass over the four inline constructs. Bold and italic are matched
# with a lazy body so that a span containing a `code` fragment - which is
# very common in this document, e.g. **Zero `unsafe`, zero warnings.** -
# is still recognised as one span and processed recursively.
SPAN = re.compile(
    r"(?P<code>`[^`]+`)"
    r"|(?P<link>\[[^\]]*\]\([^)]*\))"
    r"|(?P<bold>\*\*.+?\*\*)"
    r"|(?P<ital>(?<![\*\w])\*(?!\s)[^*]+?\*(?!\w))"
)


def inline(text):
    """Convert inline Markdown spans to LaTeX, escaping everything else."""
    pieces = []
    pos = 0
    for m in SPAN.finditer(text):
        pieces.append(esc(text[pos:m.start()]))
        pos = m.end()
        if m.group("code"):
            pieces.append("\\texttt{" + esc_code(m.group("code")[1:-1]) + "}")
        elif m.group("link"):
            lm = re.match(r"\[([^\]]*)\]\(([^)]*)\)", m.group("link"))
            label, url = lm.group(1), lm.group(2)
            if url.startswith(("http://", "https://")):
                safe = url.replace("%", "\\%").replace("#", "\\#")
                pieces.append("\\href{" + safe + "}{" + inline(label) + "}")
            else:
                pieces.append(inline(label))
        elif m.group("bold"):
            pieces.append("\\textbf{" + inline(m.group("bold")[2:-2]) + "}")
        else:
            pieces.append("\\emph{" + inline(m.group("ital")[1:-1]) + "}")
    pieces.append(esc(text[pos:]))
    return "".join(pieces)


HEADING = {1: "section", 2: "section", 3: "subsection", 4: "subsubsection", 5: "paragraph"}

BULLET = re.compile(r"^\s*[-*+]\s+")
NUMBER = re.compile(r"^\s*\d+\.\s+")


def gather_list(lines, i, n, marker):
    """Collect one Markdown list into a list of item strings.

    Items in this document are frequently several lines long and separated
    from each other by a blank line. A blank line therefore does NOT end the
    list - only a blank line followed by something that is not another item
    of the same kind does. Getting this wrong closes and reopens the list at
    every item, which restarts an enumerate at "1." each time.
    """
    items = []
    current = []
    while i < n:
        line = lines[i]
        stripped = line.strip()

        if marker.match(line):
            if current:
                items.append(" ".join(current))
            current = [marker.sub("", line).strip()]
            i += 1
            continue

        if stripped == "":
            # Look ahead: does the list continue after the blank line?
            j = i
            while j < n and lines[j].strip() == "":
                j += 1
            if j < n and (marker.match(lines[j]) or lines[j].startswith(("    ", "\t"))):
                i = j
                continue
            break

        # An indented continuation line belongs to the current item.
        if current and line.startswith((" ", "\t")):
            current.append(stripped)
            i += 1
            continue

        break

    if current:
        items.append(" ".join(current))
    return items, i


def split_row(line):
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [c.strip() for c in line.split("|")]


def convert(md):
    lines = md.split("\n")
    out = []
    i = 0
    n = len(lines)
    in_toc = False

    while i < n:
        line = lines[i]
        stripped = line.strip()

        # ---- fenced code block
        if stripped.startswith("```"):
            i += 1
            body = []
            while i < n and not lines[i].strip().startswith("```"):
                body.append(lines[i])
                i += 1
            i += 1
            # A listing body is passed through verbatim, so LaTeX escaping
            # must NOT be applied - but a stray non-ASCII character would
            # still stop the compile, so those are folded to ASCII.
            out.append("\\begin{lstlisting}")
            out.extend(ascii_fold(b) for b in body)
            out.append("\\end{lstlisting}")
            continue

        # ---- headings
        m = re.match(r"^(#{1,5})\s+(.*)$", stripped)
        if m:
            level = len(m.group(1))
            title = m.group(2).strip()
            # Drop the hand-written table of contents; \tableofcontents does it.
            if title.lower().startswith("table of contents"):
                in_toc = True
                i += 1
                continue
            in_toc = False
            if level == 1 and title.lower().startswith("part "):
                out.append("\\clearpage")
                out.append("\\section*{" + inline(title) + "}")
                out.append("\\addcontentsline{toc}{section}{" + inline(title) + "}")
            else:
                cmd = HEADING[level]
                out.append("\\" + cmd + "{" + inline(title) + "}")
            i += 1
            continue

        if in_toc:
            # skip the manual TOC body
            if stripped == "" or stripped.startswith(("*", "-", "1", "2", "3", "4",
                                                      "5", "6", "7", "8", "9")):
                i += 1
                continue
            in_toc = False

        # ---- horizontal rule
        if stripped in ("---", "***", "___"):
            out.append("\\bigskip\\hrule\\bigskip")
            i += 1
            continue

        # ---- table
        if stripped.startswith("|") and i + 1 < n and re.match(
                r"^\|[\s:|-]+\|?$", lines[i + 1].strip()):
            header = split_row(lines[i])
            i += 2
            rows = []
            while i < n and lines[i].strip().startswith("|"):
                rows.append(split_row(lines[i]))
                i += 1
            ncol = len(header)
            # Size the columns from the widest cell each one actually holds,
            # so a table of short labels does not get the same layout as a
            # table of full sentences. The weights are then normalised to the
            # text width, with a floor so no column collapses to nothing.
            if ncol == 1:
                spec = "L{0.95\\linewidth}"
            else:
                widths = []
                for c in range(ncol):
                    longest = len(header[c])
                    for r in rows:
                        if c < len(r):
                            longest = max(longest, len(r[c]))
                    widths.append(max(longest, 6))
                total = float(sum(widths))
                fracs = [max(0.10, 0.95 * w / total) for w in widths]
                scale = 0.95 / sum(fracs)
                spec = "".join("L{%.3f\\linewidth}" % (f * scale) for f in fracs)
            out.append("\\begin{longtable}{" + spec + "}")
            out.append("\\toprule")
            out.append(" & ".join(inline(h) for h in header) + " \\\\")
            out.append("\\midrule")
            out.append("\\endfirsthead")
            out.append("\\toprule")
            out.append(" & ".join(inline(h) for h in header) + " \\\\")
            out.append("\\midrule")
            out.append("\\endhead")
            for r in rows:
                r = (r + [""] * ncol)[:ncol]
                out.append(" & ".join(inline(c) for c in r) + " \\\\")
            out.append("\\bottomrule")
            out.append("\\end{longtable}")
            continue

        # ---- block quote
        if stripped.startswith(">"):
            body = []
            while i < n and lines[i].strip().startswith(">"):
                body.append(lines[i].strip()[1:].strip())
                i += 1
            out.append("\\begin{quote}")
            out.append(inline(" ".join(body)))
            out.append("\\end{quote}")
            continue

        # ---- bullet list
        if BULLET.match(line):
            items, i = gather_list(lines, i, n, BULLET)
            out.append("\\begin{itemize}")
            for it in items:
                out.append("\\item " + inline(it))
            out.append("\\end{itemize}")
            continue

        # ---- numbered list
        if NUMBER.match(line):
            items, i = gather_list(lines, i, n, NUMBER)
            out.append("\\begin{enumerate}")
            for it in items:
                out.append("\\item " + inline(it))
            out.append("\\end{enumerate}")
            continue

        # ---- blank line
        if stripped == "":
            out.append("")
            i += 1
            continue

        # ---- ordinary paragraph
        #
        # Gather the whole paragraph before converting, because **bold** and
        # *italic* spans in the source frequently run across a line break and
        # would not be recognised if each line were handled on its own.
        para = []
        while i < n:
            s = lines[i].strip()
            if s == "":
                break
            if s.startswith(("```", "|", ">", "#")):
                break
            if re.match(r"^\s*[-*+]\s+", lines[i]) or re.match(r"^\s*\d+\.\s+", lines[i]):
                break
            if s in ("---", "***", "___"):
                break
            para.append(s)
            i += 1
        out.append(inline(" ".join(para)))

    return "\n".join(out)


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "rebound_rust.md"
    dst = sys.argv[2] if len(sys.argv) > 2 else "rebound_rust.tex"
    md = io.open(src, encoding="utf-8", newline="").read().replace("\r\n", "\n")
    body = convert(md)
    stamp = datetime.date.today().strftime("%B %-d, %Y") if sys.platform != "win32" \
        else datetime.date.today().strftime("%B %d, %Y").replace(" 0", " ")
    tex = PREAMBLE.replace("DATE_PLACEHOLDER", stamp) + body + POSTAMBLE

    # Safety net: any character neither escaped nor mapped would stop
    # pdflatex with an unhelpful "Unicode character not set up for use"
    # error dozens of pages in. Report them here instead, by name.
    import unicodedata
    stray = {}
    for ch in tex:
        if ord(ch) > 127:
            stray[ch] = stray.get(ch, 0) + 1
    if stray:
        print("WARNING: %d unmapped non-ASCII character(s) remain:" % len(stray))
        for ch, count in sorted(stray.items(), key=lambda kv: -kv[1]):
            try:
                name = unicodedata.name(ch)
            except ValueError:
                name = "?"
            print("  U+%04X x%-4d %s  <- add to UNICODE in md2tex.py" % (ord(ch), count, name))

    io.open(dst, "w", encoding="utf-8", newline="\n").write(tex)
    print("wrote %s (%d bytes) from %s" % (dst, len(tex), src))
    return 1 if stray else 0


if __name__ == "__main__":
    sys.exit(main())
