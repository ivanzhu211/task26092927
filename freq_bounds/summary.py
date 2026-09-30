from __future__ import annotations

from dataclasses import dataclass

CHANNELS = ("TV", "Digital", "Social", "Print")
PAIRS = (
    ("TV", "Digital"),
    ("TV", "Social"),
    ("TV", "Print"),
    ("Digital", "Social"),
    ("Digital", "Print"),
    ("Social", "Print"),
)


def _as_nonneg_int(value, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer, got {value!r}")
    if value < 0:
        raise ValueError(f"{name} must be non-negative, got {value}")
    return value


@dataclass(frozen=True)
class Summary:
    """Aggregate reach and pairwise overlap counts for four channels.

    ``n`` includes people exposed to no channel at all.
    ``overlaps``[(a, b)] counts people covered by *both* a and b, including
    those also covered by any third or fourth channel.
    """

    n: int
    reaches: dict[str, int]
    overlaps: dict[tuple[str, str], int]

    @staticmethod
    def from_mapping(data: dict) -> "Summary":
        if not isinstance(data, dict):
            raise ValueError("summary must be a mapping")
        if "n" not in data:
            raise ValueError("summary requires key 'n'")
        raw_reaches = data.get("reaches", {})
        raw_overlaps = data.get("overlaps", {})

        reaches: dict[str, int] = {}
        for channel in CHANNELS:
            if channel not in raw_reaches:
                raise ValueError(f"missing reach for channel {channel!r}")
            reaches[channel] = _as_nonneg_int(raw_reaches[channel], f"reach[{channel!r}]")
        extra = set(raw_reaches) - set(CHANNELS)
        if extra:
            raise ValueError(f"unknown channels in reaches: {sorted(extra)}")

        def pair_key(key):
            if isinstance(key, tuple):
                a, b = key
            elif isinstance(key, str) and "|" in key:
                a, b = key.split("|", 1)
            else:
                raise ValueError(f"invalid overlap key {key!r}")
            if a not in CHANNELS or b not in CHANNELS or a == b:
                raise ValueError(f"invalid overlap channel pair {key!r}")
            order = {c: i for i, c in enumerate(CHANNELS)}
            return (a, b) if order[a] < order[b] else (b, a)

        overlaps: dict[tuple[str, str], int] = {}
        for raw_key, value in raw_overlaps.items():
            pair = pair_key(raw_key)
            overlaps[pair] = _as_nonneg_int(value, f"overlap[{pair}]")
        missing = [pair for pair in PAIRS if pair not in overlaps]
        if missing:
            raise ValueError(f"missing overlap counts for pairs: {missing}")
        if len(overlaps) != len(PAIRS):
            raise ValueError("duplicate overlap pair keys detected")

        summary = Summary(_as_nonneg_int(data["n"], "n"), reaches, overlaps)
        summary._check_local_consistency()
        return summary

    def _check_local_consistency(self) -> None:
        for channel in CHANNELS:
            if self.reaches[channel] > self.n:
                raise ValueError(
                    f"reach of {channel} ({self.reaches[channel]}) exceeds n ({self.n})"
                )
        for a, b in PAIRS:
            count = self.overlaps[(a, b)]
            if count > self.reaches[a] or count > self.reaches[b]:
                raise ValueError(
                    f"overlap of {(a, b)} ({count}) cannot exceed either channel reach"
                )
