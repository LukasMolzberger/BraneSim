"""LaTeX report generator for experiments."""

from __future__ import annotations

from datetime import datetime
from typing import Iterable

from .base import ReportData, Scalar


def _latex_escape(text: str) -> str:
    """Escape LaTeX special characters."""
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for key, value in replacements.items():
        text = text.replace(key, value)
    return text


def _format_value(value: Scalar) -> str:
    if isinstance(value, float):
        return f"{value:.6e}"
    return str(value)


def _render_table(rows: Iterable[tuple[str, str]]) -> str:
    lines = [
        r"\begin{longtable}{p{0.35\linewidth}p{0.6\linewidth}}",
        r"\toprule",
        r"\textbf{Key} & \textbf{Value} \\",
        r"\midrule",
        r"\endhead",
    ]
    for key, value in rows:
        lines.append(f"{_latex_escape(key)} & {_latex_escape(value)} \\\\")
    lines.extend([r"\bottomrule", r"\end{longtable}"])
    return "\n".join(lines)


def _render_symbol_table(rows: Iterable[tuple[str, str, Scalar, str]]) -> str:
    lines = [
        r"\begin{longtable}{p{0.25\linewidth}p{0.18\linewidth}p{0.2\linewidth}p{0.32\linewidth}}",
        r"\toprule",
        r"\textbf{Implementation} & \textbf{Paper Symbol} & \textbf{Value} & \textbf{Description} \\",
        r"\midrule",
        r"\endhead",
    ]
    for key, symbol, value, desc in rows:
        key_text = _latex_escape(key)
        symbol_text = symbol if symbol else "--"
        value_text = _latex_escape(_format_value(value))
        desc_text = _latex_escape(desc)
        lines.append(f"{key_text} & {symbol_text} & {value_text} & {desc_text} \\\\")
    lines.extend([r"\bottomrule", r"\end{longtable}"])
    return "\n".join(lines)


class LatexReportGenerator:
    """Generate a .tex report with parameters, assumptions, and notes."""

    def generate(self, report: ReportData, output_path: str) -> None:
        created_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        doc = [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=1in]{geometry}",
            r"\usepackage{longtable}",
            r"\usepackage{booktabs}",
            r"\usepackage{hyperref}",
            r"\begin{document}",
            r"\section*{Experiment Report}",
            f"\\textbf{{Title}}: {_latex_escape(report.title)}\\\\",
            f"\\textbf{{Experiment}}: {_latex_escape(report.experiment_name)}\\\\",
            f"\\textbf{{Run}}: {_latex_escape(report.run_name)}\\\\",
            f"\\textbf{{Generated}}: {_latex_escape(created_at)}",
        ]

        if report.summary:
            doc.extend(
                [
                    r"\section*{Summary}",
                    _latex_escape(report.summary),
                ]
            )

        if report.metadata:
            rows = [(k, _format_value(v)) for k, v in report.metadata.items()]
            doc.extend([r"\section*{Metadata}", _render_table(rows)])

        if report.parameters and not report.symbol_map:
            rows = [(k, _format_value(v)) for k, v in report.parameters.items()]
            doc.extend([r"\section*{Parameters}", _render_table(rows)])

        if report.choices:
            doc.append(r"\section*{Experimental Choices}")
            doc.append(r"\begin{itemize}")
            for item in report.choices:
                doc.append(rf"\item {_latex_escape(item)}")
            doc.append(r"\end{itemize}")

        if report.assumptions:
            doc.append(r"\section*{Assumptions}")
            doc.append(r"\begin{itemize}")
            for item in report.assumptions:
                doc.append(rf"\item {_latex_escape(item)}")
            doc.append(r"\end{itemize}")

        if report.derived:
            rows = [(k, _format_value(v)) for k, v in report.derived.items()]
            doc.extend([r"\section*{Derived Measurements}", _render_table(rows)])

        symbol_rows = list(report.symbol_map)
        if symbol_rows:
            doc.extend([r"\section*{Symbol Mapping}", _render_symbol_table(symbol_rows)])

        if report.figures:
            doc.append(r"\section*{Figures}")
            doc.append(r"\begin{itemize}")
            for fig in report.figures:
                item = f"{fig.caption} (file: {fig.path})"
                doc.append(rf"\item {_latex_escape(item)}")
            doc.append(r"\end{itemize}")

        if report.notes:
            doc.append(r"\section*{Notes}")
            doc.append(r"\begin{itemize}")
            for item in report.notes:
                doc.append(rf"\item {_latex_escape(item)}")
            doc.append(r"\end{itemize}")

        doc.append(r"\end{document}")

        with open(output_path, "w") as handle:
            handle.write("\n".join(doc))
