"""Quick usage examples for the frequency_bounds library.

Run:  PYTHONPATH=src python3 examples.py
"""

from frequency_bounds import Summary, bounds, is_feasible


def show(title, summary):
    print(title)
    result = bounds(summary)
    if not result.feasible:
        print("  infeasible: no integer population matches this summary")
        print()
        return
    print(f"  deduplicated reach : {result.reach.minimum} .. {result.reach.maximum}")
    for k, metric in enumerate(
        (result.frequency_1, result.frequency_2, result.frequency_3, result.frequency_4),
        start=1,
    ):
        print(f"  frequency = {k}      : {metric.minimum} .. {metric.maximum}")
    print()


# Scenario 1 from the specification.
show(
    "Scenario 1: N=4, TV/Digital/Social reach 2 (pair overlaps 1), Print 0",
    Summary(
        n=4,
        tv=2, digital=2, social=2, print=0,
        tv_digital=1, tv_social=1, tv_print=0,
        digital_social=1, digital_print=0, social_print=0,
    ),
)

# Scenario 2 from the specification.
print("Scenario 2 feasible?", is_feasible(
    Summary(
        n=4,
        tv=2, digital=2, social=2, print=2,
        tv_digital=1, tv_social=1, tv_print=1,
        digital_social=1, digital_print=1, social_print=1,
    )
))
