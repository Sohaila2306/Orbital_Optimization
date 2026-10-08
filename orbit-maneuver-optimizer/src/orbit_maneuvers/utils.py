"""Small formatting helpers for the CLI and reports."""

from __future__ import annotations

import math


def format_duration(seconds: float) -> str:
    if not math.isfinite(seconds):
        return "inf"
    if seconds < 3600:
        return f"{seconds / 60:.1f} min"
    if seconds < 48 * 3600:
        total_min = int(round(seconds / 60))
        return f"{total_min // 60}h {total_min % 60:02d}m"
    days = seconds / 86400
    if days < 2 * 365.25:
        return f"{days:.1f} d"
    return f"{days / 365.25:.2f} yr"
