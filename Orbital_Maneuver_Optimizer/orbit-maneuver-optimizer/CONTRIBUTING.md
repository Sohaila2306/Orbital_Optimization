# Contributing

Thanks for taking a look! Bug reports, corrections to the physics, and new features are all welcome.

## Getting set up

```bash
git clone https://github.com/<your-username>/orbit-maneuver-optimizer.git
cd orbit-maneuver-optimizer
python -m venv .venv && source .venv/bin/activate   # on Windows: .venv\Scripts\activate
pip install -e ".[dev]"
pytest
```

## Before opening a pull request

- Add a test. For physics changes, the best test is a number from a textbook or another trusted source
  (see `tests/test_transfers.py` for the style).
- Keep everything in SI units inside the package; convert to km/degrees only at the edges (CLI, plots, printing).
- If you add a plot, put it in `visualization.py` and make it take an optional `save_path`.
- If you change something that shows up in the README, regenerate the figures with
  `python examples/make_readme_figures.py`.

## Ideas that would be great to have

See the roadmap at the bottom of the README.
